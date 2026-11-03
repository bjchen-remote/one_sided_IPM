# ipm_cr 架构与接口

本文描述尚未部署到 NAS `main` 的 compute-only 候选。固定 patch 再分布的活动条带转移已修复，同 τ 转移、三层完整 RHS 阶数测试与实际单事件复算已有局部证据；32项 `repair_release_01` 统一检查已通过（282.715秒），源码与测试输入身份均未改变。长期独立残差降至 `1e-9` 与三段修复版复现尚未验证，不能称整套交付通过。需求见 [MAINLINE_REQUIREMENTS](MAINLINE_REQUIREMENTS.md)，操作见 [USAGE](USAGE.md)，证据与旧失败见 [VALIDATION](VALIDATION.md)。

## 1. 主链与模块

```text
main → config / problem → init → context
run → step → time.step / trial → rhs
rhs → mesh.fill / diff → context.velocity → FEM
    → gauge.wall / rates → op.prepare → 完整 F、Q
run → 周期 remesh → 有界再分布 + 活动转移 + 新 context + 一次完整 RHS
run → io.checkpoint / history → pack
```

| 位置 | 职责 |
|---|---|
| `config / problem` | 唯一方法默认、case、初值 authority、物理设置 |
| `context` | 实际几何、方法、原墙和私有 Poisson 计划 |
| `+mesh` | 轴、度量、D7、活动区、三层 ghost、转移与采样 |
| `+poisson/+fem` | 弱式、缓存求解、trace 连接和解析 Green 富集 |
| `+model / +gauge / rhs` | 背景、物理源、原 I/P/C、直接规范求解与完整 RHS |
| `+op / +time` | 中心对流、耗散、初始率区间、样条 knot 时间网格、同阶段场/尺度/时钟 |
| `remesh` | 周期源 profile、固定参考再分布与活动主态转移 |
| `+io / +diag / pack` | 身份、事务存储、被动诊断和原生数据 |

谱与 FD 后端、独立 GLL 布局和旧 namespace 保留为显式参考，不参与默认主链或继承当前资格。

## 2. API 与主态

```matlab
opts = ipm_cr.config(overrides);
[s,ctx] = ipm_cr.init(opts);
[s,info] = ipm_cr.step(s,ctx,tauEnd,timeControls);
f = ipm_cr.rhs(s,ctx);
s = ipm_cr.refresh(s,ctx);
v = ipm_cr.sample(s,ctx,X,Y);
[s,newCtx,event] = ipm_cr.remesh(s,ctx,settings);
[summary,s,ctx] = ipm_cr.run(c);
```

`step` 完成一次显式推进，内部 knot 分段属于时间网格构造；完整循环归 `run`。`timeControls` 对应 CFL、`maxDt` 和时间格式，没有 `minDt`。`refresh` 重建同一数学身份的 RHS；方法改变须建立新运行或显式 seed 导入。

| 状态 | 含义 |
|---|---|
| `s.W / s.eta` | 完整中心化密度 / 可选解析背景余项，`W=eta+H` |
| `s.rhs` | 同态 q、速度、规范率、F/Q 及诊断 |
| `s.scale` | 尺度、平移、密度原点、规范时间与完整物理时钟 |
| `s.originalAuthority / referencePeak` | 固定初始参考，不随网格重设 |
| `s.remesh` | 上次周期时刻、次数和固定 `referenceCore` |
| `s.import` | 显式 seed 来源、SHA、真实 τ 和转移比较 |

`F=W_τ`、`Q=q_τ`。解析背景路径先微分 η 并加回解析背景导数；对完整 W/F 再差分一般不同于定义的 q/Q，该表示差可记录为诊断。物理密度为 `exp(-logC_omega)*(densityOrigin+W)`。

物理时钟保存 `physicalTime + physicalTimeLow + physicalTimeTail(1:2)` 四字展开。阶段增量、背景时间项与 RK 场和规范标量使用同一 tableau；CSV 的前两字不替代完整时钟。

## 3. 共同主要节点与连接

`medium` 从 `data/mesh.json` 读取实际轴，`coarse` 每轴精确取 `1:2:end`，度量步长相应倍增。核心数组行 Y/列 X 为 `129×241 / 65×121`，环带行 r/列 θ 为 `549×193 / 275×97`，外半径 `1e10`。

要求是 Poisson 与对流主要内部节点完全相同；边界低维例外须单独说明。当前实现全部原生节点共用：`x/y/r/theta`、`X/Y` 和 shape 逐位一致，无额外全局 GLL 未知量。Gauss 点仅作单元内积分。

`mesh.build` 支持 preset、显式生成器或完整 `axes={x,y,r,theta,dxMetric,dyMetric,drMetric}`。实际轴重建微分、度量、活动区、donor、交换、采样和 geometry key。三层 ghost 的活动边界按实际节点确定。

物理微分 `mesh.diff` 的活动行只读取同一 patch 的活动节点：内部保留原中心七点模板，人工接口附近采用同一物理坐标上的七点单侧/偏心闭合。core 的墙法向、polar 的两条物理墙和外圆闭合保持原模板。五阶迎风仍使用同步 ghost，活动重叠值不投影。这个选择避免把两片不同的 RHS 截断误差通过 ghost 交换后再作横跨接口的差分；它真实改变了离散源 `q=D_A eta+H_X`，Poisson 与 `Q=D_A eta_tau+H_Xtau` 使用同一个 `D_A`，不是单独更改残差显示。geometry key 记录 `physicalDerivativeClosure=same_patch_active_D7_v1`，旧轨道不继承为同方法资格。光滑场的局部多项式一致性、单侧权重放大、完整反馈及稳定性须分别测试。

两片活动副本独立演化。`mesh.fill` 同时由旧活动 donor 填 ghost，不投影活动重叠区；`mesh.sample` 以活动 core 为唯一 owner，执行局部七点采样，其边缘 stencil 可以读到同步 ghost。remesh 转移另行只读活动 donor；采样 owner 相同不表示每个重构 donor 都活动。局部参考值加差分的算式保留常数消去；MEX 与分块 MATLAB 实现使用同一公式，缓存只含几何与查询。

## 4. FEM 与外域

默认 native C0 Lagrange 元，`degree=5`、`quadOrder=8`、`stride=[]`。几何动态规划分组控制局部插值条件，少量六次元吸收余数，最外完整径向元为六次；不依赖源或测得误差选组，也不增加节点。非空 native `stride` 仅兼容 `stride=degree`。

核心在物理 x/y、环带在 `t=log(r),theta` 装配弱式。原生 raw `psi` 直接取节点系数；速度来自同一 FE 场的解析单元导数。Gauss 源求积和跨 patch trace 仍需传递。固定几何缓存稀疏 Cholesky 因子和小 Schur 系统；每次 patch 线性解固定作一次原弱残差修正，不按误差挑选迭代次数。构建的角向特征系统、trace、填充与 remesh 重建成本须独立量测。

墙取 `psi=0`、`u2=0`，允许非零上侧墙源，奇延拓可有跳跃。含时 FEM 显式使用有限源圆盘外置零及齐次 DtN；该 compact-source 模型以总阶三角点 jets 构造20系数的全局多项式 Green 场，FEM 解余项。完整场在原生节点及迹位置加回解析项。直接组件的已知完整外源 Green Dirichlet 问题使用其声明模型。

弱行、trace、重叠和 Green 误差是不同诊断，不是生产接受门。Green 比较须同源、同尺度、同截断模型；未知物理外源保持 `outerSourceTailCertified=false`。C0 流函数不自动给出连续或 C1 速度，共享面平均也不消除两侧极限的跳跃。

## 5. RK4 时间网格与规范

默认 `linear5` 五阶迎风、平滑速度系数 `.1`、经典 RK4；`central4`＋KO6 `.2` 为显式中心参考选项。实际非均匀度量进入输运与速率，平滑分裂与所选耗散是离散定义。外圆按已选多项式闭合计算，入流速度与数据模型仅记录为 metadata；本题的物理模型适用性在正式实验前独立确认。

每次先使用当前已完成态的 RHS 率选取初始区间

```text
h = min(maxDt, 0.95*cfl/currentRhs.rate, tauEnd-currentTau)
```

`currentRhs.rate` 包含输运、所选稳定化和幅度反应率。原 I/P/C 的自然三次样条峰使方程在原测量墙 knot 之间分片光滑。分段 RK4 对普通相邻 cell 的连续峰端点跨越，用固定4次 secant/Newton 更新构造到该 knot 的积分段，然后继续余段；所有定位探测从同一起态积分。

该构造确定数学分片的时间网格，不按误差或残差接受状态。没有 stage-only guard、残差阈值、括根 fallback、折半、`minDt`、拒绝/重试或旧态择优。ULP 只识别同一坐标，不把峰或状态投影到 knot。非连续的全局峰分支交换仍直接采用原 global maximizer，不按连续 knot 事件定位。

各段内部阶段和新末态照常计算完整 RHS；阶段率、CFL 和规范/Poisson 质量诊断被动记录。末段完整积分剩余 τ，不恢复旧四ulp软上限规则。空间稳定化与四字时钟公式保持原定义。

`compute_only_release_04` 的时间失败保留：端点差在 `τ≤.11` 约 `1e-15`，在 `.12/.34` 跨 knot 后约 `1e-12`；光滑 `.10/.11` 约 `4.01` 阶。固定四次弦截分段随后已冻结，并通过原 triplet 的 temporal-only 复测及 `compute_only_release_05`，W 时间观测阶3.9606。这是对应冻结版本的时间证据；当前 remesh 修复版的统一检查仍待完成。

I/P/C 保持原墙离散定义。gauge 对 `(c_l,delta)` 用精确速度 Jacobian 直接 Newton，最多8次更新，点不再改变即停；没有条件数门、线搜索、备用 seed 求根或非唯一峰拒绝。并列峰记录被选全局最大值和分支导数，残差与条件数被动记录。零能量比值按公式未定义，不填假零。

## 6. 核心形态累积变化触发 remesh

配置为 `enabled / level / threshold / checkEveryTau / minCheckEveryTau / monitorStrength / blend / maxNodeMotion / candidateAxes`，默认见 [USAGE](USAGE.md#形态触发-remesh)。取全活动正 q 峰值25%超水平集中的主峰连墙分量，按横向最大绝对范围和纵向最大范围计量；相对最后实际 remesh 后参考的任一方向累积变化达到20%才插值。首次仅校准，检查间隔仅控制观测时机，不重置累积参考或强制换网格；按旧形态变化率与阈值余量在 `.005… .5τ` 间调度，保存lastShape及下一检查时刻。几何触发与测试专用的精度验收门严格区分。

每次由当前源 profile 更新固定参考核心上的目标密度，作一次有界凸移动，保持点数、patch 外框、原点和环带轴。固定 `referenceCore` 避免反复乘密度。有 η 时转移 η，否则转移 W，再重建 mesh/context 并计算一次完整 RHS。

固定完整轴端点和节点数的再分布保留各片独立表示：在本片实际活动轴上用最近七点张量多项式重构，新暴露的边界条带继续本片端部多项式，不从另一片替换活动值，不读旧 ghost，也不投影活动重叠。`localMap`/`tensorMap` 的端部延拓须显式启用，默认仍禁止外推。改变 patch 框或节点数的导入使用覆盖路由，属于另一种表示转移。

随后直接继续，没有候选精度验收、反向转移、W/q/F/Q、速度、IPC、重叠门或回滚。原参考、尺度、四字时钟和步号由复制与数据流保留，在独立测试中验证。构造或算式本身无法计算时由运行器记录错误，不切换替代路径。

七点多项式次数不是完整 remesh 精度证明。独立 `test_remesh_order` 用三层非零解析 η，在同一新网格/context 上比较转移态与直接解析构造态的完整 W/q/F/Q；舍入受限的误差对不赋阶。当前通过仅限该平滑、被分辨的同 τ 试验；未分辨窄层与长期轨道仍需分别解释。

## 7. 复现与存储

当前窗口为 medium `0→5`、medium `390→395`、精确半密度 coarse `0→200`。T2 必须使用有来源的真实 τ390 态；三段绑定源码、方法、原初值和规范。结果保存后独立与已有数据比较，报告实际网格、参考底限和方法差异，不通过拟合时间、平移或振幅消除差异。

context key 绑定实际网格、原墙、方法、初值和源码；`opts.geometry` 保留初始处方。恢复从 checkpoint 的实际 `geometry` 冷建计划，不恢复因子 handle。`run.seed={file,sha256,tau}` 是新运行导入：原墙/authority/背景须相容，保留尺度和时钟，重建目标 RHS，转移比较被动记录。

owner 锁、SHA、读回与原子指针属于文件身份及事务存储。`accepted_*` 表示已完成并提交的数值态，不表示误差合格。`run` 保存实际 mesh、事件谱系、原生场、完整时钟和停止原因；出错不自动重试。停止只有目标、声明预算/步数、物理时间、用户 STOP 或错误，没有在线稳态窗口。

既有 Python 监督器冻结 src、mesh、authority、选项及 seed 副本，按单线程分段、外部资源预算和累计 deadline 执行；恢复不自动延长预算。数值质量不进入该生产循环。

## 8. 测试与证据边界

独立测试调用同一生产接口，检查公式、制造解、阶数、稳定性、转移、时钟和冷恢复；独立 `test_steady` 不成为实验早停。正式演化只计算与存储，复现比较另行进行。

旧 `b890305` 的31项统一检查、未分段 suite04 的时间失败和 suite05 的默认空间/时间屏保留原身份。修复前周期 remesh 的粗轨道扰动、实际事件与转移误差原样保留；当前修复消除了一个复算事件的十多倍 Q 跳变，不等于长期恢复到旧 `1e-9` 残差。早先多门 T1/T3 的 owner STOP 和9次 roundtrip 拒绝仅描述旧机制。当前不宣称发布、C1 速度、近线性成本或未知外源尾认证；GLL/P5、Q4、谱、FD与旧均匀轨道按其原配置解释。
