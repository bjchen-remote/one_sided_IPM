# `ipm_cr` 架构与接口契约

本文说明 `main` 分支唯一维护求解器 `ipm_cr` 的模块边界、数据、缓存和失败行为。使用入口与证据索引见 [README](../README.md)，离散推导与实验解释见 [算法报告](../test/ALGORITHM_REPORT.tex)，网格连接构造及误差预算见 [CONNECTION_METHOD](CONNECTION_METHOD.md)。

## 1. 模块边界与公开入口

```text
src/+ipm_cr/
  config / problem                 完整数值选项 / 物理初值与规范
  init / step / rhs / refresh       建立模拟 / 一个接受步 / 完整RHS / 同身份冷恢复
  sample / run                     通用物理采样 / 可恢复轨道
  context                          每个模拟独立的只读上下文
  pconfig / psolver                 独立谱 Poisson 配置与管理
  +mesh                            网格、ghost、差分和局部采样
  +model                           初值、authority、解析背景及尾项
  +gauge                           原 I/P/C、峰值与规范求根
  +poisson                         源/速度映射计划及 core/polar 谱算子
  +time                            RK4 试步、阶段 CFL、峰事件与回滚
  +op                              固定场重构、速度应用及原生内核
  +io                              checkpoint、history、身份、JSON、SHA
  +diag                            独立重叠审计与研究区域比较
```

公开主链为：

```matlab
opts = ipm_cr.config(overrides);
[s,ctx] = ipm_cr.init(opts);
[s,info] = ipm_cr.step(s,ctx,tauEnd);
eval = ipm_cr.rhs(s,ctx);
s = ipm_cr.refresh(s,ctx);
v = ipm_cr.sample(s,ctx,x,y);
[summary,s,ctx] = ipm_cr.run(runOptions);
```

`init` 直接构造状态，`rhs` 直接组装本阶段背景与完整 RHS；没有第二套初始化或 RHS 分发入口。`step` 只返回一个接受步，`tauEnd` 是该步上限；完整时间循环由 `run` 管理。`rhs` 返回当前状态的完整求值，可用第二输出 `[eval,W]` 取得同步、重建后的完整中心化密度。`refresh` 在同一数学身份下调用它重建接受态 RHS；两者都不跨版本转换轨道。`sample` 接受匹配的规范坐标查询数组，不要求区域掩码。区域比较、采样范数和研究云在 `diag.compare`，不侵入通用采样。

初始化、RK 阶段、诊断、history 和 v2 checkpoint 全程使用同一套 `W/eta/rhs` 字段；不再往返转换 `rho/residual/flow`。旧 v1 数据的命名转换仅在 `io.checkpoint('import',folder)` 这一只读边界发生；`'legacy'` 仍返回未经修改的历史 payload。导入不改写历史文件，也不授权新源码继续旧轨道。

`config` 是完整数值选项的唯一默认来源，保持 `geometry/solverOptions/transport/cfl/maxDt/minDt/timeIntegrator/background` 字段。`pconfig` 只定义谱椭圆选项，不与含时 `config` 混用；`problem` 集中物理初值、原始规范和 authority。`data/initial.json` 是唯一初始文件定位与 SHA 清单，实际文件为 `data/reference/gaussian_t0.mat`。

`model.initial` 仅实现原 Gaussian 初值，运算顺序与原 authority 一致；旧求解器的其他算例不混入维护实现。背景选项控制同一物理初值的解析分解，不改变原初值及 I/P/C 参考。

`model.reference(cfg)` 直接读取已核 SHA 的原始 t0 frame，返回原墙算子、初始尺度、I/P/C 参考和规范初猜。原墙节点来自该 frame，原 `Dx`、密度/墙源一致性和 t0 时钟均重新核验；建立参考不构造输运网格，也不建立旧 Cartesian/polar 状态载体。

`context` 是每次模拟独立的 handle，公开 `mesh/wall/opts/key/id/identity` 只读；固定 Poisson 计划和 `psolver` 句柄由它私有持有；墙采样 `wallR/wallMap` 对外只读，统计经 `ctx.stats()` 读取。不把每个局部函数改成类，不持有跨状态密度/源/规范观测；数组仍以 struct 和数值矩阵传递。模型、时间步与诊断只能使用匹配身份的状态和上下文。

`tools/run_explicit_experiment.py` 和 `tools/run_validation.py` 通过 `tools/common.py` 统一冻结/校验规则；`tools/build_ipm_cr_native.m` 从源码构建 MEX。`tests` 只存验证与研究适配，`examples/poisson_demo.m` 是独立椭圆演示，不推进含时系统。

旧 `test/+ipm_cr`、`time_evolution` 和旧运行器在本地保留但不再作为维护依赖或默认 MATLAB 路径；所需旧公共数值函数已收进相应新包，完整旧算法仍在 [`legacy/solvers`](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers)。`test/prototypes/`、`experiments/`、历史结果及冻结 runtime 保持原位，不能自动清理或借当前包替换旧实验身份。

谱 Poisson 与输运是两套网格。谱配点的端点聚集不进入输运 CFL，输运环带不含原点。`mesh.localMap/tensorMap` 保留当前连接与源映射必需的局部插值。

## 2. 网格和数组

| 对象 | 朝向、作用域与不变量 |
|---|---|
| `g.core` | 行 `Y`、列 `X`；均匀张量网格；默认完整域 `[-1.2,1.2]×[0,1.2]`、步长 `.03`；左右及顶部三层为 ghost，物理墙 `Y=0` 为活动点 |
| `g.polar` | 行 `r`、列 `θ`；`θ∈[0,π]`；内侧三层为 ghost，其余均为活动点；外半径默认 `1e10`，墙列 `Y=0` |
| `g.exchange` | 两侧 ghost 的七点张量插值矩阵和常数差分计划；donor 只来自另一侧活动点 |
| `g.geometryKey` | 实际 `x/y/r/θ/drMetric` 及 ghost、插值宽度；用于源、返回和墙查询缓存身份 |
| `psolver.Core/Polar` | 核心 Chebyshev 张量与分片 `log(r)`/角向谱网格，仅供 Poisson |

两块活动重叠区各自独立演化，不能把一侧活动值投影到另一侧。`mesh.fill` 同时读取交换前两侧的活动 donor 并提交 ghost，不循环读取新 ghost。交换插值使用

\[
f(q)=f_*+\sum_j w_j(f_j-f_*),
\]

以局部 donor `f_*` 为锚，常数场的差值严格为零；原几何的只读矩阵兼容路径保留，当前公开恢复仍要求相同数学身份。`diag.overlap` 在独立非节点查询上分别重构两块活动值，记录 `rho/q/F/Q/u` 差异；ghost 相等不是独立活动区一致性的证明。

### 固定远区径向重分配

`mesh.polar` 仅构造基础轴、坐标和径向度量。旧原点拟合、二阶差分、重复一阶矩阵及整张 `R/Theta` 已无活动消费者；实际输运的 `Dr/Dt` 在 `mesh.build` 截取环带后构造。基础轴和度量的算术保持原样。

`opts.geometry.radialRedistribution.enabled` 默认 `false`，关闭时不改变原轴与 `geometryKey`。开启时，`mesh.radial` 在 `uniform_core` 原连续径向轴的旧索引上积分 BetaCDF(8,8) 监测函数，再反解新整数索引。默认保护 `r≤4`，在 `4<r<20` 平滑过渡，外侧索引密度约增至 `factor=2`。固定网格全程不动，不引入 ALE 或含时重映射。

保护前缀直接复制旧节点；只有七点模板完全处于该前缀的度量行复制原 `radialMetric/(baseNr-1)`，其他度量由新整数索引上的七点差分求得。截取环带后设

\[
J_r=\Delta r_{\rm metric}(n_r-1),\qquad h_\xi=1/(n_r-1).
\]

输运的逻辑导数除以 `Jr`，CFL 径向间距为 `drMetric`。物理导数 `Dr` 直接按实际 `r` 构造，不再除 `Jr`。`gridOptions.nr` 为基准请求，`polar.originalNr/polar.nr` 才是重分配后的完整轴/截取轴大小。新增点改变矩阵列步幅，保护检查比较物理 donor 和权重。

`g.radialRedistribution` 保存新旧索引、插入数、节点和度量保护范围及实际节点数。`mesh.build` 对完整 donor 模板和首活动度量的保护范围做断言，保护半径不足时停止。监测函数为 C7，积分映射及正导数逆映射为 C8；具体公式、整数取整及加密族约定见 [连接方法](CONNECTION_METHOD.md#可选固定远区径向重分配)。

## 3. 状态、背景和完整场

| 数据 | 定义及持久化要求 |
|---|---|
| `s.W` | 完整中心化密度 `W`，核心原点严格为零 |
| `s.eta` | 实际积分变量 `η=W-H`；启用解析背景时必须保存 |
| `s.scale.densityOrigin` | 单独保存的密度常数 `a`，完整密度 `R=a+W` |
| `s.scale` | `logC_l/logC_omega/X_shift`、规范时间、物理时间高低部分 |
| `s.originalAuthority/ctx.wall` | 原初值 I/P/C 参考及原墙节点，不随网格重定 |
| `s.rhs` | 对应同一接受态、尺度、方法的完整 RHS、源、速度、规范和诊断；禁止将另一状态的 RHS 作为有效缓存 |
| `ctx.opts` | 几何、谱阶数、重构/CFL/步长、时间格式、背景的数值身份 |

默认表示为

\[
R=a+W,\quad W=\eta+H,\quad H=G-G(0,0),\qquad
G=e^\omega[\rho_0(x,y)+tA(x,y)],\quad
(x,y)=((X-S)e^{-\ell},Ye^{-\ell}).
\]

`model.gaussian` 给出原 Gaussian primitive 的稳定中心化值和解析导数；`model.tail` 给出原初值物理主尾的光滑解析项 `A`；`model.background` 组合二者。`tailCorrection=true` 为默认，初始 `t=0,η=0`。尾项表示有限物理时间的已知解析部分，未预设最终重标度稳态尾部；`tailCorrection=false` 可关闭尾项，`background.kind='none'` 直接积分 `W`。

`model.background` 返回 `H/x/y/xx/xy/originX/timeH/timeX`。`timeH/timeX` 已含 `dt_physical/dτ`；调用方不能再乘一次。物理时间使用 `physicalTime/physicalTimeLow` 两部分累加和求尾项，不能只保存一个合并时钟。背景及所有参数均属于检查点恢复身份。

`mesh.diff` 按当前传入的物理导数矩阵计算 `D_X/D_Y`，每个模板先形成 `f_j-f_i` 再乘系数，避免密度常数进入差分。当前几何不持有物理差分模板缓存；诊断替换实际矩阵时直接使用新矩阵。常数消去不修复输入本身已经丢失的浮点有效位。

## 4. 一次完整 RHS 的契约

`ipm_cr.rhs(s,ctx)` 在一个入口中创建本阶段背景并完成下列绑定计算：

1. 同时同步 `η` ghost，重建完整 `W=η+H`。
2. 计算残差物理导数并同步其源 ghost，组装 `q=D_Xη+G_X`。
3. 将残差源映射到独立谱节点；解析 `G_X` 在谱节点直接求值，严格求解 `-Δψ=q`，返回同一谱解的 `u=(-ψ_Y,ψ_X)`。
4. 在原墙样本上测量完整 `W` 的 I/P/C；查询 `Uref=u1(C,0)`，求解当前重构对应的 `c_l/c_r/c_omega`。
5. 用真实输运与阶段背景组装 `η_τ` 和完整 `F=W_τ`；独立审计实际 `F` 的 I/P/C 变化率，记录完整 CFL。

接受末态额外构造 `Q=D_Xη_τ+H_{τX}`。内部阶段可省略展示用 `ψ` 返回和 `Q`，但不能省略完整源、速度、规范、实际 `F` 审计和 CFL。`s.rhs.F/s.rhs.q/s.rhs.Q` 分别指完整 `W_τ/W_X/F_X`，不能用残差 RHS 或残差源代替。

`op.prepare/apply` 在一个 RHS 内固定当前密度空间重构，供规范迭代多次查询速度及偏导；闭包不跨状态保存场。默认 `linear5` 为五阶最优迎风线性重构；`smooth_weno5` 使用平滑归一化 WENO-Z；二者要求兼容 MEX。`legacy_weno5` 保留当前包内的兼容 MATLAB 算子。内核构建记录源码/二进制 SHA 与严格浮点选项，关闭 fast-math 和收缩。

平滑分裂的默认 `σ=.1 min(dx,dy)`，使用 `hypot(v,σ)`；实际分裂速度进入 CFL，墙法向通量为零。离散模式与速度系数是方法选项，不能为提速静默切换。有限外边界的五次外延仍需独立远场敏感性验证。

## 5. 规范与时间推进

`gauge.wall` 保留原 I/C 和全窗口自然三次样条峰值 P：先对 `q-max(q)` 求同一最大值，再加回锚值，避免近等高带引入峰位置平台。P 的变化率用同一峰的权重。`gauge.rates` 对固定空间重构的实际输运构造 `(c_l,δ)` 解析偏导，暖启动 Newton/线搜索；原 `gauge.seed` 仅在需要时提供新种子，最终仍按当前完整泛函求根。

Newton 目标为 `16eps`，线搜索停滞且残差不超过 `2e-12` 时允许有记录的有界例外；这不等于达到目标。最终真实完整 RHS 的 I/P/C 审计独立进行，失败停止，不能用求根器残差替代。

`time.step` 默认经典 RK4；`time.rk4` 同时供生产推进与含时制造解使用。密度与以下标量使用同一阶段值和 RK 权重：

\[
\ell'=c_l,\quad\omega'=c_\omega,\quad S'=c_lS+c_r,
\quad t'=e^{\omega-\ell},\quad a'=F_0+c_\omega a.
\]

已接受初始 RHS 复用，计算三个内部阶段及新末态完整 RHS；末态留给下一步。物理时间正增量使用补偿累加。峰跨原样条节点时，用同一接受态出发的真实 RK4 子步定位并分步，保留各阶段 CFL，不投影密度或峰值。`timeIntegrator='ssprk3'` 仅为显式兼容选项；RK4 没有 SSP 保证。

试步限制为两域**全部活动点**上的真实径向/角向或 Cartesian 速率加 `abs(c_omega)`；只有 ghost 不进入最大值。预测留 5% 余量，初始、三个内部阶段与新末态都必须满足 `dt*rate≤cfl`。超限回滚缩步；可恢复求解失败重试，耗尽重试或低于最小步长后保存失败证据并停止。CFL 审计不替代带高阶连接的稳定性证明。

## 6. 谱 Poisson 与缓存身份

| 缓存 | 所有者与失效条件 |
|---|---|
| 核心特征分解、环带逐模 LU、trace Schur LU/齐次响应 | `psolver`；谱几何、边界或配置变化时重建 |
| 残差源映射与谱速度返回映射 | `context` 私有持有 Poisson plan；实际网格身份或完整谱配置变化时重建 |
| 原墙查询映射 | `context` 拥有并只读暴露；实际网格或原墙节点变化时重建 |
| ghost 交换计划 | `mesh.build`；由真实活动 donor、坐标和插值宽度决定 |
| 密度重构、背景公共原值、规范速率闭包 | 仅单 RHS/单调用，不能跨状态复用 |

核心为 `A_Y Ψ+Ψ A_X^T=q`；环带在 `t=log(r)` 上求解 `-ψ_tt-ψ_θθ=r²q`，面板匹配值和一阶径向导数。小 Schur 系统连接两域 trace。缓存齐次响应允许“特解+响应”重建；原离散后向误差不满足门限时，回退直接核心/环带求解并重做原审计，不放宽门限。

`psolver.reconfigure` 先完成候选几何/缓存再提交；失败保留旧对象，新几何拒绝旧解。计数器记录构建、复用和直接回退。每次保留原方程分量后向误差、Schur/实际 trace、独立重叠 `ψ/u/法向导数` 审计；有独立连续源时才计算离网格强残差。数组残差源叠加解析背景仍没有完整连续源真值，该强残差明确为不可用 `NaN`。`requireAccuracy=true` 默认严格验收。

默认外边界 `decay` 是无源外域的离散角向 DtN：`ψ_t+sqrt(μ_k)ψ=0`。它没有证明域外源可忽略，`outerSourceTailCertified=false`；`dirichlet` 需提供外弧流函数。有限外半径与源尾需独立审计。

## 7. 只读输出与运行事务

`ipm_cr.sample(s,ctx,x,y)` 是完整背景路径的通用物理采样。接受态必须带有匹配的 `W/eta/rhs/scale/step`，`x/y` 形状匹配且有限；域外查询行为按采样实现检查。背景在查询点解析求值，仅余项与原生速度插值，`F/Q` 含 `timeH/timeX`。完整物理速度和去平移速度分列；它不推进、不对齐、不改活动值。

`diag.compare` 才使用 `near/interfaceBand/wall/remoteWallStrip/hotspot` 这些比较区域，保留原三层研究的采样与诊断算术。区域掩码不成为求解器状态或一般采样的必需参数；有限查询范数不是连续区域上界。

`ipm_cr.run` 在每段开始核验维护源相对路径与 SHA、初始输入及完整数值选项。接受 checkpoint 顶层保留 `state/geometry/explicitOptions/wallOps`，状态使用公开字段；还保存 history、来源和必要的恢复信息。context 的私有固定计划不直接序列化，恢复时按相同数学身份冷重建；禁止省略 `eta`、尺度、尾项参数、时钟低位或原始规范。

双槽写完、读回并核验 SHA 后才原子更新 `latest.json`，指针中的 `path` 仅为槽文件 basename，失败文件不成为恢复指针。源码身份采用相对路径，允许完整冻结目录在不改变输入内容的前提下迁移；这不允许更换代码或方法后续接旧轨道。

Python 管理器冻结当前 `src`、工具、选项、`data/initial.json` 和其引用输入，再监督同一 `ipm_cr.run` 入口。不遍历复制本地保留的旧求解器或原型，冻结和 MATLAB 入口必须使用一致的源码/输入身份。预算、人工停止、资源中止、求解失败和目标达到分别记录；累计 MATLAB 成本与含启停的监督成本分别报告。旧结果和冻结 runtime 保留其原 SHA，不迁移或重写历史实验身份。

## 8. 维护与证据分工

根 `README` 提供当前公开 API、运行示例和验证地图；本文说明模块与数据契约；`CONNECTION_METHOD` 维护连接构造与误差预算；仍位于 `test/ALGORITHM_REPORT.tex` 的算法报告保留推导、实验数字和适用范围，并更新当前源码索引；`CLEANUP` 记录历史退出维护面的闭包、归档与架构迁移。旧 `test/README.md` 只作导航。

迁移前数值维护统一了公共读出并消除重复采样。差分模板缓存候选虽通过逐位回归，中层局部计时约慢 6%，已恢复原差分算法，候选记录留在 [本地历史证据](../test/results/cleanup_20261022)。当前缓存所有权重构与目录迁移须独立验收；它们不继承旧轨道的源码身份，也不声称已得到整机加速。

验证按改动选择：几何/常数/制造解 → 连接和完整 RHS → 短程时间/空间对照 → 固定端点谱加密 → 完整轨道。迁移验收应额外覆盖上下文独立性、暖/冷计划一致性、相对身份、干净副本运行与 checkpoint 冷恢复。名义五阶重构、四阶 RK、很小的线性残差或两条数值轨道相近，均不能单独认证全局精度。既有 τ10 有限采样观测及缺失控制保持原义，当前 τ200 未验证。

## 9. 工程管理原则的出处

参考 PETSc 的 [DM](https://petsc.org/release/manual/dmbase/) 与 [setup/solve](https://petsc.org/release/manualpages/KSP/KSPSetUp/) 分工：网格离散提供求解数据，准备阶段与逐次求解分开。参考 SUNDIALS 的 [SUNContext](https://sundials.readthedocs.io/en/latest/sundials/SUNContext_link.html) 和 [ARKODE 返回状态](https://sundials.readthedocs.io/en/latest/arkode/Usage/User_callable.html)：每个模拟有独立上下文，事件、目标到达与失败有不同含义。参考 AMReX 的 [Geometry/FillBoundary](https://amrex-codes.github.io/amrex/docs_html/Basics.html) 与 [checkpoint](https://amrex-codes.github.io/amrex/docs_html/IO.html)：几何和场分离，ghost 填充不改活动值，应用自行定义完整恢复数据。

这些官方资料支持职责分离与生命周期原则；本项目的包名、context 接口和文件组合是针对 MATLAB 双网格方法的工程选择，不是上述框架规定的实现。此处不引入框架依赖；RK4 算术、CFL、原峰泛函和真实子步事件定位由本项目的数值契约决定。

`io.checkpoint('save',...)` 由持有运行目录 owner 锁的单一写者调用，写前核对传入的上一代指针与磁盘一致；过期游标会在覆盖任何槽之前被拒绝。正常 `run` 自动管理该锁和指针链。


## 发布配置与数据边界（2026-10-24）

`config()` 保留完整方法接口；`config('case')` 返回唯一用户 case 模板，含 initial/domain/method/run/stop/save，`config(case,'run')` 统一解析后供生产 runner 使用。旧 mesh/poisson 和 explicitOptions 仅在此边界归一化，重复输入拒绝；run 不再维护另一份默认值或别名逻辑。唯一默认解析背景来自选中 manifest 的 Gaussian 参数，init 复用它判断严格零 eta 初态。

`problem(initial)` 和 `io.identity('authority',root,initial)` 绑定项目相对 manifest 与 MAT 的 SHA。`make_initial` 为 Gaussian 参数族生成新的固定测量轴 authority 和 I/P/C。context 在准备前后核验源码及初值身份，将该身份纳入只读 key；pack 必须核对创建时身份，不能用导出时的输入把旧状态重新命名。

内部时间步仅为 `time.step(state,ctx,finalTau,controls)`；controls 限四项时间策略。几何、原墙和完整方法由 ctx 提供，不存在缺失 ctx 时临时重建的维护分支。数值 RK、峰事件、全阶段 CFL 及 Poisson 直接回退不变。

`diag.steady` 是可选的接受态观测，默认关闭；完整原生活动节点和全部重叠副本的 F/Q、固定网格全场/速率窗口变化共同生成候选。tracker、原生 anchor、最后观测和计数随 checkpoint 保存；策略改变明确重新累计，恢复候选不得额外推进。它不属于 RHS，也不参与阶段演化。

`pack` 在新的 data 子目录保存不含求解器句柄的原生科学包、方法/初值/源码身份、真实坐标/active/ghost、完整场与尺度时钟。论文重采样独立于原生 checkpoint 和误差验证。使用细节见 USAGE；测试清单唯一来源为 tools/validation.json。
