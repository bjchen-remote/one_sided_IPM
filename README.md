# ipm_cr：单侧 IPM 混合网格求解器

默认运行入口为 `main.m`，使用非均匀中网格：核心 `129×241`、环带 `549×193`。Poisson 与对流共用全部主节点；粗网格只将四个轴各隔点取样。默认 FEM、`linear5` 五阶迎风、RK4 和核心形态触发 remesh，`central4` 中心格式保留为显式选项。

本版包含独立 patch 高阶 remesh、累计形态触发及人工接口的活动点微分闭合。接口修复的12项相关回归通过；所测解析标量对流接口再次求导的有限对阶为4.80/4.08。完整非线性轨道、全部派生量的连续精度及整套主干验收仍按 [验收记录](docs/VALIDATION.md) 分别报告，源码发布不等于这些结论已全部成立。

本轮复现实验窗口为中网格 `τ=0→5`、中网格 `390→395`、精确半密度粗网格 `0→200`，重点与已有数据作独立对照。需求及最新指令见 [主干契约](docs/MAINLINE_REQUIREMENTS.md)；计算、独立测试和结果比较分别安排。

[使用](docs/USAGE.md) · [架构/API](docs/ARCHITECTURE.md) · [FEM 方法与历史证据](docs/FEM_POISSON.md) · [参考后端](docs/SAME_MESH_POISSON.md) · [验收记录](docs/VALIDATION.md) · [分支与复现](REPOSITORY.md)

## 当前数值方法

| 部分 | 候选实现 |
|---|---|
| 几何 | 核心内部非均匀加密＋不含原点的环带；三层 ghost，活动重叠副本独立 |
| 中 / 粗网格 | 核心数组 `129×241 / 65×121`，环带 `549×193 / 275×97`；粗轴精确取 `1:2:end` |
| Poisson | 原生节点 C0 FEM，`degree=5`、`quadOrder=8`、`stride=[]`；几何分组包含少量六次元，最外完整径向元为六次；缓存稀疏 Cholesky＋小 trace 系统 |
| 外域 | 显式有限源圆盘外取零＋外域 DtN；该 compact-source 模型采用解析 Green 角点富集 |
| 对流 | 默认 `linear5` 五阶迎风、平滑速度系数 `.1`；`central4`＋KO6 `.2` 保留为显式选项 |
| 物理微分 | 内部中心 D7；人工接口活动行采用本 patch 七点单侧/偏心闭合，q 与 Q 同算子 |
| 时间 / gauge | RK4＋原墙样条 knot 固定四次弦截分段；初始区间按当前 RHS 率选取；原 I/P/C 精确速度 Jacobian、最多8次直接 Newton 更新 |
| remesh | 核心主峰连墙区域的方向尺度累积变化达20%才重分布；固定 patch 外框内保留各片独立活动重构，端部七点多项式继续；新 context 和一次完整 RHS |

Poisson 与对流的主要内部节点完全共用。当前实现进一步共用全部原生 `x/y/r/theta`、物理 `X/Y` 和数组 shape；没有独立的全局 GLL 未知量。Gauss 点仅作局部求积，原生源到求积点和跨 patch 迹仍有传递。raw FEM 的节点 `psi` 直接取系数，速度由同一有限元基函数求导；富集时系数表示余项，完整场加回解析 Green 场。

初始积分区间为

```text
h = min(maxDt, 0.95*cfl/currentRhs.rate, tauEnd-currentTau)
```

默认 `cfl=.3`、`maxDt=.2`。原 I/P/C 使用全局自然三次样条峰，连续峰跨原墙 knot 时方程分片光滑。分段 RK4 对普通相邻 cell 的连续端点跨越使用固定4次 secant/Newton 更新，积分到该 knot 后继续余段；探测均从同一起态积分。非连续的全局峰分支交换按原 global maximizer 直接演化。

这是原方程的时间网格构造，精度判断仍由独立测试完成。没有残差阈值验收、stage-only guard、括根 fallback、折半、`minDt`、拒绝重试或旧态择优；ULP 仅识别同一坐标，不投影峰或状态。阶段 CFL、规范残差和 Poisson 质量诊断被动记录；平滑速度项与 KO6 保持所选空间公式。

## MATLAB 与库入口

根目录 `main.m` 顶部可编辑 case。默认 `rho0=[]`，使用 `data/initial.json` 登记的原 Gaussian；默认目标 `τ=200`、本次预算600秒、输出为项目根 `result`。

```matlab
[c,rho0] = main('config');       % 只读，不构建或求解
c.run.tauEnd = 5;
c.run.outputDir = 'result/medium_0_5';
[summary,s,ctx] = main(c);
```

`main(c,rho0)` 也接受向量化初值表达式。表达式按实际网格物化并绑定 SHA，重新计算原墙参考；改变初值或初始网格使用新输入与目录。直接库入口为：

```matlab
addpath('src');
opts = ipm_cr.config();
[s,ctx] = ipm_cr.init(opts);
[s,info] = ipm_cr.step(s,ctx,.003);
f = ipm_cr.rhs(s,ctx);
v = ipm_cr.sample(s,ctx,[0;.3],[0;.1]);
```

完整循环只有 `ipm_cr.run(case)`；`main` 是其路径与输入包装。采样可用同公式的 MEX 或分块 MATLAB 实现；默认 `linear5` 与显式 `smooth_weno5` 需要 WENO MEX。`main` 可构建本机内核，二进制不入库；Windows 实机运行与恢复尚未验收。

remesh 的形态触发与重分布配置由 `ipm_cr.remesh('config')` 给出，见 [使用说明](docs/USAGE.md#形态触发-remesh)。固定节点数与完整轴端点时，新暴露条带继续本 patch 最近七活动点的多项式，不读 ghost 或切换到另一片活动副本；通用采样默认仍不外推。完成转移与一次 fresh RHS 后直接继续，没有误差触发阈值、候选验收或回滚。原参考、尺度和全部时钟在测试中核验。

## 复现、停止与存储

| 实验 | 网格与实际窗口 |
|---|---|
| T1 | medium，从登记原始初值 `0→5` |
| T2 | medium，从有来源的真实 `τ=390` 态到 `395` |
| T3 | 原四轴各减半密度的 coarse，从相同物理初值 `0→200` |

晚期 seed 显式声明 `file/sha256/tau`，建立新运行并记录跨方法/网格导入；不能称为同方法恢复或从零演化。与已有数据的比较在保存结果后独立进行，匹配实际 τ、初值、规范、尺度和时钟，不拟合平移、振幅或时间掩盖差异。

`stop` 只含 `physicalTime`。目标、调用预算、累计步数和用户 `STOP` 管理运行时长；没有 `steady_candidate` 数值早停。既有 Python 监督器负责冻结、单线程分段、资源预算和累计绝对截止，`resume` 不自动增加预算。

| 输出 | 内容 |
|---|---|
| `latest.json`、`accepted_0/1.mat` | 原子提交的已完成态与 SHA；文件名中的 accepted 不代表精度认证 |
| `initial.mat/json` | 本次真实起态及 SHA，滚动检查点不覆盖 |
| checkpoint、`mesh/` | 完整主态、实际网格、原墙参考、方法/输入、四字时钟、形态参考控制器和 remesh 谱系 |
| summary、history、`data/` | 真实终点、停止原因、被动诊断、成本及原生 W/η/q/F/Q/速度 |

`W` 是中心化密度，完整规范密度为 `densityOrigin+W`；解析背景路径积分 η。`F=W_τ`、`Q=q_τ`，不能用余项 RHS 或对完整 F 再差分替代 Q。显示重采样与原生数据、独立误差参考分开保存。

## 独立测试与旧证据

生产轨道不执行误差门。组件、制造解、空间/时间阶、稳定性、转移和冷恢复测试调用同一生产 API，在正式实验之前独立安排；稳态分析保留为独立 `test_steady`，不进入实验停止逻辑。已有静态 Green 证据按其冻结来源引用。

修复前的 coarse 复算曾在同一 remesh 时刻把完整活动点 `Q∞` 从 `.0029181` 推至 `.0365613`；修复后该事件变为 `.0029181→.0029185`。真实 η 往返误差从 `6.65e-3` 降至 `2.82e-5`。这些局部证据与三层非零 η 的完整 RHS 阶数测试分别记录在 [验收记录](docs/VALIDATION.md)，不替代长期复现。七点次数不保证未分辨窄层精度或每一项小截断误差均下降。

`compute_only_release_04` 完成31项 requested checks，但时间资格为 false。本次诊断中新旧端点差在 `τ≤.11` 约 `1e-15`，`.12/.34` 跨自然样条 knot 后升至 `1e-12`；光滑 `.10/.11` 的时间观测阶仍约 `4.01`。该失败保留。数学分段已通过原 triplet 的独立复测及完整 suite05；F 和部分标量低于分辨屏，仍不赋予观测阶。

旧 `b890305` 的31项统一检查耗时 `237.63 s`；当时非零墙角源的完整速度/角点速度观测阶为 `4.2927–4.2960 / 4.0242–4.1302`，RK4 的 W 时间观测阶为 `3.9566`。这些只证明旧冻结配置的相应有限测试；旧 radial 完整空间资格仍为 false，不能转移为当前 compute-only 源码或新复现窗口的通过记录。

旧多门版本 T1 正常停止于约 `τ=1.2878`，T3 正常停止于约 `τ=23.297`、3680步，保存了旧往返门9次拒绝；两者为 owner STOP，不能记为新版本成功或数值失败。旧独立 GLL/P5 的角点阶数失败，以及 Q4、谱、FD、旧均匀网格轨道与 frozen runtime，保留各自范围。尚无本候选整套发布、全局 C1 速度、通用 `O(N)` 成本或未知物理外源尾的认证。
