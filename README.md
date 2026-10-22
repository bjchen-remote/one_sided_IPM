# ipm_cr：单侧 IPM 混合网格求解器

`main` 维护 `ipm_cr`：均匀 Cartesian 张量核心、不含原点的极坐标环带、独立谱 Poisson 求解，以及全部活动点、全部阶段受物理 CFL 约束的显式 RK4。解析背景与残差分离用于控制远场表示及舍入误差。

[架构与接口](docs/ARCHITECTURE.md) · [网格连接](docs/CONNECTION_METHOD.md) · [算法报告](test/ALGORITHM_REPORT.tex) · [分支与复现](REPOSITORY.md) · [迁移记录](docs/CLEANUP.md) · [本轮验收](docs/VALIDATION.md)

## 安装与短程计算

需要 MATLAB、Python 3；原生内核需要兼容的 C++ 编译器。MATLAB 在项目根目录执行：

```matlab
restoredefaultpath;
addpath('src','tools');
build_ipm_cr_native();       % 从源码在本机生成内核，二进制不入库
opts = ipm_cr.config();      % 完整数值设置的唯一默认入口
[s,ctx] = ipm_cr.init(opts);
[s,info] = ipm_cr.step(s,ctx,.03);  % 一个接受步，且不越过 tauEnd
v = ipm_cr.sample(s,ctx,[0;.3],[0;.1]);
```

`step` 返回一个已接受步；它不是推进到指定时间的循环。`tauEnd` 是该步的上限；完整轨道交给 `ipm_cr.run`。采样只需要匹配的物理查询 `x/y`，不要求研究区域掩码。`s.W` 是完整中心化密度，`s.eta` 是实际积分余项，`s.rhs` 保存该接受态的完整 RHS 与诊断，`s.scale` 保存尺度及补偿物理时钟。完整密度是 `s.scale.densityOrigin+s.W`，不能将 `eta` 当作完整密度。`ctx` 持有网格、原墙参考、方法与私有固定计划；墙采样映射对外只读，计数经 `ctx.stats()` 读取。

独立椭圆演示不推进含时系统：

```matlab
addpath('examples');
[~,~,diagnostic] = poisson_demo();
assert(diagnostic.accepted);
```

## 完整轨道、停止与恢复

下面使用既有 τ10 实验的中层方法。远区径向重分配必须显式开启；未设置时使用原网格。默认 Poisson 次数为 `48/64/16`，示例使用已测的 `64/80/20`。

```matlab
rr = struct('enabled',true,'factor',2, ...
    'preserveRadius',4,'fullRefinementRadius',20);
g = struct('coreSpacing',.03, ...
    'grid',struct('nr',385,'ntheta',129,'firstSpacing',.03), ...
    'radialRedistribution',rr);
p = struct('coreDegree',[64,64],'angularDegree',80, ...
    'polarDegree',20,'rOuter',1e10,'requireAccuracy',true);
o = struct('outputDir',fullfile(pwd,'output/my_tau10'), ...
    'tauEnd',10,'wallSeconds',600,'geometry',g, ...
    'solverOptions',p,'timeIntegrator','rk4','cfl',.3,'maxDt',.2);
[summary,s,ctx] = ipm_cr.run(o);
% 同一源码、同一数值设置、同一目录再次调用，继续最后已接受态。
[summary,s,ctx] = ipm_cr.run(o);
```

`outputDir` 使用绝对路径。`wallSeconds` 是本次调用预算；`wall_budget` 表示已安全保存后返回，只有 `target_reached` 表示达到目标。在输出目录创建 `STOP` 会在接受步边界保存退出；恢复前删除自己创建的 `STOP`。改变源码、初始输入或数值设置时使用新目录，不能绕过身份检查续接旧轨道。`ipm_cr.refresh(s,ctx)` 只在同一数学身份下重建 RHS，不负责跨版本 checkpoint 转换。

长实验由 Python 管理器冻结输入并监督同一入口。先把数值设置写成 JSON，再执行：

```sh
python3 tools/run_explicit_experiment.py start \
  --output "$PWD/output/my_frozen_run" \
  --options /absolute/path/to/options.json \
  --target 10 --hours 1 --segment-seconds 600
python3 tools/run_explicit_experiment.py status \
  --output "$PWD/output/my_frozen_run"
```

管理器冻结 `src`、所需工具、选项和初始输入，并核对相对路径的源码身份。单线程串行运行，限制 RSS、可用内存、磁盘和总时间；MATLAB 段成本与监督进程总成本分别记录。`resume` 沿用原任务的绝对截止时间；额外预算必须独立记录，不能自动延长原账本。MATLAB 路径按本机设置，部分资源监控仅适用 macOS。

## 输出与验证

| 输出 | 含义 |
|---|---|
| `summary.json`、`summary.mat` | 实际状态、终点、累计耗时和方法摘要 |
| `latest.json` | 最后接受检查点的槽文件名、序号和 SHA-256 |
| `accepted_0.mat`、`accepted_1.mat` | 交替槽；只读取指针指定且 SHA 正确的槽 |
| checkpoint 的 `state/geometry/explicitOptions/wallOps` | 主态、实际网格、数值选项和原始墙规范；固定计划可冷重建 |
| `history/historyColumns` | 每步实际 CFL、重试、Poisson、I/P/C、峰事件与耗时 |

`ipm_cr.rhs(s,ctx)` 求当前状态的完整耦合 RHS，`ipm_cr.sample(s,ctx,x,y)` 只读重建完整场；带区域掩码的三层研究比较由 `ipm_cr.diag.compare` 负责。物理速度和去平移速度分开，`q=W_X`、`F=W_τ`、`Q=F_X` 不与残差导数混用。

只读网格绘图为 `examples/plot_ipm_cr_mesh.m`，通过统一 checkpoint API 读取当前 v2 或显式只读导入历史 v1，不重建 RHS；绘图所需辅助函数随源码提供。

测试从 `tests` 运行，输出写入新的独立目录。原生构建入口在 `tools`，演示与只读绘图在 `examples`；不再把整套求解器放入测试目录。

| 验证对象 | 入口 |
|---|---|
| 公开状态、上下文身份与通用采样 | `test_arch` |
| 断点事务与目录迁移 | `test_runtime` |
| 原生内核 | `test_ipm_cr_native` |
| 谱 Poisson、线性叠加和直接回退 | `test_ipm_cr`、`test_ipm_cr_superposition` |
| 背景、解析尾项、完整规范 | `test_ipm_cr_background`、`test_ipm_cr_tail_background`、`test_ipm_cr_gauge_acceleration` |
| 连接与固定径向重分配 | `test_ipm_cr_connection_operators`、`test_ipm_cr_radial_redistribution` |
| RK4 制造解与真实短窗时间加密 | `test_ipm_cr_fourth_order` |
| 原初值三层轨道 | `test_ipm_cr_ivp_spatial` |
| 维护等价回归 | `test_ipm_cr_refactor`，需另行提供原 baseline MAT |
| 固定端点谱加密与成本诊断 | `inspect_ipm_cr_explicit_endpoint`、`profile_ipm_cr_rhs` |

`tests/run_ipm_cr_explicit.m` 是旧研究签名的适配入口，生产调用使用 `ipm_cr.run`。历史脚本和冻结 runtime 按其原身份运行，不为迁移后包生成新的历史资格声明。

## 源码与数值证据边界

`src/+ipm_cr` 是完整的维护源码；不依赖旧 `+ipm/+ipm_c` 包。`config` 管完整数值选项，`pconfig/psolver` 只管独立椭圆系统，`problem` 集中物理初值和规范。初始文件由 [data/initial.json](data/initial.json) 定位并校验 SHA。完整旧算法在 [`legacy/solvers`](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers)；本地保留旧目录、结果和原型不属于维护源码。

既有固定远区重分配 τ10 三层、全部 `15975` 公共查询上的自收敛，`W/q/F/Q/u1-Uref/u1/u2` 观测阶为 `4.577/4.714/4.479/4.00569/4.619/4.506/4.078`。`Q` 刚超过四阶；中层时间缩步与固定最细端点谱加密已完成，最细层缩步、整条轨道谱加密、连续全域误差上界及当前路径 τ200 尚未建立。原报告保留 `fullIVPShortWindowQualified=false`。

既有三层接受步数为 `764/1523/3115`，MATLAB 累计成本 `83.56/260.30/1349.08 s`；单步成本增加，没有建立等误差速度优势。此次架构迁移不新增四阶资格、加速或 50 分钟完成 τ200 的保证。数学推导、原始数据数字和历史源身份见 [算法报告](test/ALGORITHM_REPORT.tex) 与 [连接方法](docs/CONNECTION_METHOD.md)；本地证据链接在干净克隆中可能不存在。
