# `ipm_cr` 混合网格求解器

`test/+ipm_cr` 在均匀 Cartesian 张量核心和不含原点的极坐标环带上推进密度，用独立的谱网格求解 Poisson 方程。当前含时路径为经典显式 RK4、解析背景与残差 `η` 拆分、全部活动点与全部阶段的物理 CFL 限制。原 `time_evolution/+ipm`、`time_evolution/+ipm_c` 提供物理模型、原初值和 I/P/C 参考，本目录不修改它们。

最新已完成的证据是启用固定远区径向重分配的 **τ10 三层采样自收敛**：`W/q/F/Q/u1-Uref/u1/u2` 的观测阶分别为 `4.577/4.714/4.479/4.00569/4.619/4.506/4.078`。`Q` 刚超过四阶，需要更细一级检验渐近趋势。中层缩步和固定末态 Poisson 加密均已完成；最细层缩步、连续全域误差上界、新路径 τ200 及等精度速度比较尚未建立。

## 运行与恢复

以下命令从项目根目录执行，采用最新 τ10 实验的中层方法。径向重分配为显式选项，未设置时仍使用原网格；默认 Poisson 次数为 `48/64/16`，示例为已测 `64/80/20`。

```matlab
addpath('test','test/tests','time_evolution');
if isempty(which('ipm_cr.explicit.wenoPairSmoothMex'))
    build_ipm_cr_native();
end
rr=struct('enabled',true,'factor',2, ...
    'preserveRadius',4,'fullRefinementRadius',20);
geom=struct('coreSpacing',.03, ...
    'grid',struct('nr',385,'ntheta',129,'firstSpacing',.03), ...
    'radialRedistribution',rr);
opts=struct('outputDir',fullfile(pwd,'test/results/my_redistributed_run'), ...
    'tauEnd',10,'wallSeconds',600, ...
    'geometry',geom,'solverOptions',struct('coreDegree',[64,64], ...
    'angularDegree',80,'polarDegree',20,'rOuter',1e10), ...
    'timeIntegrator','rk4','cfl',.3,'maxDt',.2);
summary=run_ipm_cr_explicit(opts);
% 同一源码、同一数值设置、同一目录再次调用，继续最后已接受态。
summary=run_ipm_cr_explicit(opts);
```

`outputDir` 必须是绝对路径。`wallSeconds` 是本次调用的时间预算，到预算后正常写检查点并返回 `wall_budget`；只有 `target_reached` 表示达到目标。若要停止，在输出目录创建 `STOP`，运行器在接受步边界保存退出。后续恢复前删除自己创建的 `STOP`。改变源码、authority 或方法选项时使用新目录，不能绕过身份核验续接旧轨道。

长实验使用同一入口的 [Python 管理器](run_explicit_experiment.py)，先把 MATLAB 示例中的数值选项写成 JSON，再运行：

```sh
python3 test/run_explicit_experiment.py start \
  --output "$PWD/test/results/my_frozen_run" \
  --options /absolute/path/to/explicit_options.json \
  --target 10 --hours 1 --segment-seconds 600
python3 test/run_explicit_experiment.py status \
  --output "$PWD/test/results/my_frozen_run"
```

管理器冻结包、入口、配置及原始 authority，单线程运行，约束 RSS、可用内存、磁盘和总时间，记录每段 MATLAB 与监督进程耗时。`resume` 使用原任务记录中的绝对截止时间，过期任务不会因再次调用自动延长。额外预算应单独记录，保留原任务账本。

## 读取数据与绘图

| 文件 | 用途 |
|---|---|
| `summary.json`、`summary.mat` | 实际状态、终点、累计耗时与方法摘要 |
| `latest.json` | 最后已接受检查点的文件名、序号和 SHA-256 |
| `accepted_0.mat`、`accepted_1.mat` | 交替槽；读取指针指定且 SHA 正确的槽 |
| 检查点内 `state/geometry/explicitOptions/wallOps` | 主态、实际节点、方法及原始规范；缓存可以重建 |
| `history/historyColumns` | 每步实际 CFL、重试、Poisson、I/P/C、峰事件与耗时 |

公共只读接口 [`ipm_cr.explicit.sampleState`](+ipm_cr/+explicit/sampleState.m) 从接受态、实际网格、方法及带区域掩码的公共查询云重建完整场，供精度比较复用。

`state.rho` 是完整中心化密度 `W`，`state.residual` 是积分变量 `η`，完整密度为 `state.scale.densityOrigin + W`。`flow.q=W_X`、`flow.F=W_τ`、`flow.Q=F_X`；`flow.physical.u1/u2` 为绝对 Darcy 速度，`u1-Uref` 另行表示中心速度。绘图需同时保留当前尺度、背景和实际网格，不能将 `η` 当成完整密度。

```matlab
addpath('test/plots');
plot_ipm_cr_mesh('/absolute/path/to/fresh_mesh_figures', ...
    '/absolute/path/to/checkpoint_folder');
plot_ipm_cr_redistribution_convergence([],[],[], ...
    fullfile(pwd,'test/results/my_redistribution_figures'));
```

第一个函数显示已接受态的全局 `rho` 与核心、连接区 `rho_x`；第二个只读指定的新旧三层 JSON，重画网格、自差、公共查询和成本对比，默认输入为下表中的 τ10 实验。两者都要求独立输出目录。展示重采样不替代原生或公共查询精度指标。

## 源码与验证地图

| 位置 | 职责 |
|---|---|
| [`tests/run_ipm_cr_explicit.m`](tests/run_ipm_cr_explicit.m) | 唯一含时运行入口、身份检查、历史、双槽检查点 |
| [`+explicit/initialize.m`](+ipm_cr/+explicit/initialize.m)、[`advance.m`](+ipm_cr/+explicit/advance.m) | 原初值、状态、RK4 和同步标量推进 |
| [`+explicit/geometry.m`](+ipm_cr/+explicit/geometry.m)、[`redistributeRadialAxis.m`](+ipm_cr/+explicit/redistributeRadialAxis.m) | 活动/ghost 网格、固定远区径向重分配和几何计划 |
| [`+explicit/synchronize.m`](+ipm_cr/+explicit/synchronize.m)、[`differentiate.m`](+ipm_cr/+explicit/differentiate.m)、[`transport.m`](+ipm_cr/+explicit/transport.m) | 同时 ghost 交换、常数消去差分、五阶输运 |
| [`+explicit/background.m`](+ipm_cr/+explicit/background.m)、[`evaluateResidual.m`](+ipm_cr/+explicit/evaluateResidual.m)、[`evaluate.m`](+ipm_cr/+explicit/evaluate.m) | 解析背景、完整源/RHS 与规范耦合 |
| [`+explicit/poisson.m`](+ipm_cr/+explicit/poisson.m)、[`Solver.m`](+ipm_cr/Solver.m) | 输运/谱网格映射、固定分解和严格椭圆验收 |
| [`+core`](+ipm_cr/+core)、[`+polar`](+ipm_cr/+polar)、[`+internal`](+ipm_cr/+internal) | 核心/环带谱算子与内部几何工具 |
| [`+carrier/localMap.m`](+ipm_cr/+carrier/localMap.m)、[`tensorMap.m`](+ipm_cr/+carrier/tensorMap.m) | 当前连接与源映射仍使用的局部插值 |

采用与改动对应的分层验证，所有输出使用独立目录。下列入口不会自动组合为一次全局资格认证。

| 目标 | 维护入口 |
|---|---|
| 原生重构内核构建与局部检查 | [`build_ipm_cr_native`](tests/build_ipm_cr_native.m)、[`test_ipm_cr_native`](tests/test_ipm_cr_native.m) |
| 谱 Poisson 制造解与离散线性等价 | [`test_ipm_cr`](tests/test_ipm_cr.m)、[`test_ipm_cr_superposition`](tests/test_ipm_cr_superposition.m) |
| 背景算术、解析尾项、完整规范 | [`test_ipm_cr_background`](tests/test_ipm_cr_background.m)、[`test_ipm_cr_tail_background`](tests/test_ipm_cr_tail_background.m)、[`test_ipm_cr_gauge_acceleration`](tests/test_ipm_cr_gauge_acceleration.m) |
| 连接和径向重分配 | [`test_ipm_cr_connection_operators`](tests/test_ipm_cr_connection_operators.m)、[`test_ipm_cr_radial_redistribution`](tests/test_ipm_cr_radial_redistribution.m) |
| RK4 输运制造解与真实短窗时间加密 | [`test_ipm_cr_fourth_order`](tests/test_ipm_cr_fourth_order.m) |
| 原初值三层完整轨道 | [`test_ipm_cr_ivp_spatial`](tests/test_ipm_cr_ivp_spatial.m) |
| 维护重构等价回归 | [`test_ipm_cr_refactor`](tests/test_ipm_cr_refactor.m)；检验完整 RHS、公共采样与短 RK 轨道的逐字段等价 |
| 固定状态源、Poisson 与成本诊断 | [`inspect_ipm_cr_explicit_initial`](tests/inspect_ipm_cr_explicit_initial.m)、[`inspect_ipm_cr_explicit_endpoint`](tests/inspect_ipm_cr_explicit_endpoint.m)、[`profile_ipm_cr_rhs`](tests/profile_ipm_cr_rhs.m) |

`linear5` 为默认重构，`smoothSpeedCoefficient=.1`。`smooth_weno5` 使用同一内核的平滑 WENO-Z 权重，`legacy_weno5` 保留旧算子诊断路径；选项变化需要单独验证。构建不选择重构模式，要求 MEX 的模式缺少兼容内核时明确报错。`main_cr.m`、`run_validation.py` 仅为独立椭圆入口；SSPRK3 不属于当前默认 RK4 主链。`prototypes/` 是仍在使用的并行稳态加速研发区，跨目录依赖与研究记录保持独立，本轮不清理或改写。

## 结果与文档

| 证据 | 可支持的结论 |
|---|---|
| [重分配静态检查](results/radial_redistribution_static_20261022_v2/radial_redistribution.json) | 核心、halo、donor、保护度量一致及算子精度；一对误差低于诊断阈值，未作完整空间资格 |
| [τ10 三层报告](results/radial_redistribution_tau10_20261022/study/ivp_spatial.json)、[修复分析](results/radial_redistribution_tau10_20261022/fix_result.json) | 同一原初值、全 `15975` 公共查询的采样自收敛和运行成本 |
| [中层缩步对照](results/radial_redistribution_tau10_20261022/time_check_completed.json) | 当前中层 CFL `.3→.18` 的时间贡献条件筛查；未覆盖最细层 |
| [固定 τ10 端点谱加密](results/radial_redistribution_tau10_20261022/study/endpoint_poisson_tau10_v2.json) | 同一最细接受态下 `64/80/20→96/128/24` 的敏感性；未重推进轨道 |
| [收敛与网格图](results/radial_redistribution_tau10_20261022/figures/redistribution_convergence_v3/redistribution_convergence.pdf) | 来自保存 JSON 的已审阅展示、CSV 与来源清单 |

三层接受步数为 `764/1523/3115`，与原网格相同；MATLAB 累计耗时为 `83.56/260.30/1349.08 s`，单步成本增加。该结果没有建立等误差速度优势，也不保证 50 分钟内完成 τ200。原始三层报告的 `fullIVPShortWindowQualified=false` 和最细层时间控制缺失保持原义，后续分项控制不改写冻结报告。

本轮整理保留 `sampleState` 公共读出并消除空间研究中的重复采样。差分模板缓存候选虽然通过逐位等价回归，中层局部计时却由旧实现中位数 `3.96452 ms` 增至 `4.20885 ms`，因此已撤回；候选源码与记录保存在 [本轮维护证据](results/cleanup_20261022)。该微基准不代表完整 RHS 或整机加速。

算法推导和验证边界见 [算法报告](ALGORITHM_REPORT.tex)；模块、数据和缓存契约见 [架构说明](ARCHITECTURE.md)；网格连接构造与误差预算见 [连接方法](CONNECTION_METHOD.md)；历史代码归档与恢复见 [清理记录](CLEANUP.md)。

历史 SSPRK3 [τ200 轨道](results/explicit_tau200_20261021_v1/native/summary.json) 及 [极坐标对照](results/explicit_vs_polar_tau200_20261021/comparison.json) 仅属于各自冻结版本。其 MATLAB 累计耗时为 `43.97/86.99 min`，配置和误差未匹配；它们不能作为当前 RK4/η 路径的精度或等精度加速结论。历史结果、检查点、冻结副本和原初值 authority 均作为证据保留。
