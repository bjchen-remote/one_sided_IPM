# 运行 case 与论文数据

## 直接运行根目录 main.m

在 MATLAB 编辑器打开根目录 `main.m`，修改顶部 case 和 `rho0=@(x,y)...` 并点击 **Run**。默认表达式为 `@(x,y)(1+tanh(x)).*exp(-y.^2)`，目标τ=200、调用预算600秒、外半径 `1e10`、输出为项目根 `result`。入口选择 h=.03、基准 nr=385/ntheta=129、远区重分配 factor=2、Poisson次数64/80/20；完整方法仍由 `ipm_cr.config` 解析。达到预算时保存接受态，不保证达到时间目标。

```matlab
[c,rho0] = main('config');       % 无求解、无 MEX 构建，返回 case 和表达式
c.run.tauEnd = .003;
c.run.wallSeconds = 120;
c.run.outputDir = 'result/short run';
rho0 = @(x,y)(1+tanh(x)).*exp(-y.^2);
[summary,s,ctx] = main(c,rho0);
% 保持目录、源码、方法和初值，再次调用即继续最后接受态。
c.run.tauEnd = .006;
[summary,s,ctx] = main(c,rho0);
```

`main()` 使用顶部 case 和函数表达式。`[c,rho0]=main('config')` 只返回二者，不建立网格、上下文或数值轨道。`main(c)` 使用 `c.initial` 文件；`main(c,rho0)` 接受任意向量化表达式，返回 `[summary,state,ctx]`。表达式使用 `.*`、`./`、`.^` 等逐元素运算，须返回与坐标同形状、有限、实数密度。

```matlab
[c,rho0] = main('config');
rho0 = @(x,y)(1+tanh(2*x)).*exp(-y.^2);
c.run.outputDir = 'result/my_initial';
[summary,s,ctx] = main(c,rho0);
```

`tools/set_initial.m` 在实际 core/polar 节点及原固定参考轴上物化表达式的真实数值，保存独立 SHA 输入并重算 I/P/C。不持久化函数句柄，不进行 Gaussian 拟合。该路径使用 `background.kind='none'`，直接积分完整中心化 W，状态没有 η，不能宣称初始 η=0。物化输入绑定输运网格，换网格须用原表达式重建。改变初值、网格或方法时换输出目录；身份冲突会拒绝且保留旧数据。当前入口只接函数表达式，不提供外部离散数组接口。

`main` 按自身文件位置找到 `src/tools/data`。相对输出路径以项目根目录定位，绝对路径保留；调用不通过 `cd` 改变工作目录。默认算法需要本平台 MEX，缺失时自动从 `wenoMex.cpp` 构建。若 MATLAB 尚未选择 C++ 编译器，先按错误提示执行 `mex -setup C++`。构建失败会停止，不切换到另一种输运算法。直接运行和续算不需要 Python。

Windows 示例：解压至 `C:\work\one_sided_IPM`，打开其根目录 `main.m` 后运行。命令窗口从其他目录调用时，可先执行 `addpath('C:\work\one_sided_IPM')`。保留所有 MATLAB 包目录名中的 `+`，包括 `src\+ipm_cr`。输出可使用含空格的盘符绝对路径或 UNC 路径；生成文件名不含冒号。`c.initial` 仍使用 `data/.../initial.json` 这样的项目相对路径。现有验收尚未包含 Windows 实机的编译、运行和恢复。

## 高级 case 接口

配置的默认来源为 `ipm_cr.config`。`ipm_cr.config()` 保留数值方法接口；`ipm_cr.config('case')` 生成供用户编辑的 case，`ipm_cr.run` 接受该结构或 JSON 文件。高级接口由调用者准备 MATLAB 路径和 MEX，其相对输出路径仍以当前工作目录定位。`examples/case.json` 是明确选择中层网格和谱次数的 τ10 示例，不是第二份默认配置。

```matlab
addpath('src','tools');
build_ipm_cr_native();
c = ipm_cr.config('case');
c.run.tauEnd = 200;               % 规范时间 τ 的绝对终点
c.run.wallSeconds = 12*3600;      % 本次调用的墙钟预算
c.run.outputDir = 'result';       % 默认；相对当前 MATLAB 工作目录
c.domain.rmax = 1e10;             % 默认；同时解析输运与谱外半径
c.save.label = 'resolution_h003';
[summary,s,ctx] = ipm_cr.run(c);
% 或：ipm_cr.run('examples/case.json')
```

case 分为六块：

| 块 | 用户设置 | 含义 |
|---|---|---|
| `initial` | 默认 `data/initial.json` | 项目相对 manifest；SHA 绑定实际初值和固定原墙 I/P/C |
| `domain` | `rmax=1e10` | 输运和 Poisson 的共同外半径；重复设置不一致时拒绝 |
| `method` | 数值覆盖值，默认空 struct | `geometry/solverOptions/transport/background/cfl/maxDt/minDt/timeIntegrator` |
| `run` | `tauEnd=200`、`wallSeconds=600`、`outputDir='result'` | 演化终点、调用预算、输出和断点频率 |
| `stop` | 物理时间阈值、默认关闭的稳态观测 | 在完整接受态检查；与方法身份分开 |
| `save` | `data=true`、`history=true`、`label=''` | 原生科学包、CSV 和可辨认实验标签 |

`method` 中可选核心步长/尺寸、极坐标点数/径向映射/固定远区重分配、谱次数、重构模式、解析背景与尾项，以及时间控制。未覆盖的模块参数沿各构造器的验证规则解析。不存在通用“论文精度”开关；改变数值方法后需对应的时间、网格、Poisson 和远场敏感性实验。退役的 `originRings/originFitDegree` 明确报错。

## 停止与恢复

`run.maxSteps` 默认 2,000,000，表示自 t0 起的累计接受步数上限。`tauEnd` 是 RK 步的准确上限；物理时间由高低两部分补偿时钟保存，`stop.physicalTime` 非空时在首次越过阈值的接受端点停机，并报告越界量，不伪装成精确物理时间事件。`wallSeconds` 也在接受步边界检查，初始化、完整单步和最终持久化可能越过该软预算；严格总截止由外部监督器管理。

| `summary.status` | 解释 |
|---|---|
| `target_reached` | 达到有限 τ 目标 |
| `physical_time` | 接受态首次到达物理时间阈值 |
| `steady_candidate` | 预登记的完整残差与连续窗口条件满足 |
| `wall_budget`、`maximum_steps` | 预算或累计步数用尽，已保存接受态 |
| `stopped` | 在输出目录发现用户 `STOP` 文件 |
| `failed` | 验收、身份、执行或输出失败；已发布接受指针不被失败态替换 |

同一源码、方法和初值下，再次调用原 case 继续 `latest.json` 指定的接受态；`run.resume=false` 可明确拒绝已有轨道。用户创建的 `STOP` 需由用户明确移除后才能恢复。不要将不同版本源码直接接到旧检查点上。改变停止策略时，保存 `stopPolicyReset` 记录并重新累计窗口，不改初值参考或已计算轨道。

```matlab
c.stop.steady.enabled = true;
c.stop.steady.residualTol = 1e-6;
c.stop.steady.changeTol = 1e-5;
c.stop.steady.windowTau = 10;
c.stop.steady.windows = 3;
c.stop.steady.sampleEveryTau = 1;
c.stop.steady.minTau = 200;
```

观测保留 near（默认 r≤3）、body（r≤200）和 global 的全部原生活动节点，包括两块独立演化的重叠副本。完整 `F=W_τ` 以 `max(1,||W||∞)` 归一化，完整 `Q=q_τ` 以 `max(1,||q||∞)` 归一化；两项在三个区域均须通过。每个非重叠窗口至少跨 `windowTau`，比较固定网格完整 W/q/F/Q、未去平移的规范速度、去平移速度和全部五个规范速率；中间任一观测残差失败也使窗口失败。只在新闭合的通过窗口报告候选，不凭重复观察累计次数。观测频率之间的连续时间行为没有因此得到证明。

候选不是连续误差、全域精度或网格收敛证书。`eta` 的 RHS、接近机器精度的 I/P/C 残差、有限时间目标和近场通过均不替代这些条件。CFL、Poisson 原验收及直接回退照常执行。

## 初值

主入口的任意表达式通过 `set_initial` 物化，采用无背景 W 路径。独立 Gaussian 工具仍可使用匹配的解析背景，此时才有严格零 η 初态。

独立的 `make_initial` 工具用于 Gaussian primitive 参数族，不限制 `main(c,rho0)` 的表达式。它生成独立 authority，保留固定原测量轴、零时钟、单位尺度，并重算 I0/P0/C0：

```matlab
addpath('tools');
[initial,receipt] = make_initial('wide_x',struct('widthX',1.0,'skewWeight',.15));
c = ipm_cr.config('case');
c.initial = initial;
```

该独立工具可调 `amplitude/centerX/widthX/widthY/skewWeight`，拒绝覆盖已有 case。`data/cases/` 默认不入 Git；共享问题时须同时保留 JSON、MAT 和 SHA。manifest/MAT 路径均相对项目根，禁止外逃。任意表达式使用上述 `main(c,rho0)`；不接受外部离散数组或 η₀ 数组。新初值不继承旧轨道或旧精度结论。

## 数据与复现

默认目录 `result` 包含接受指针、双槽 MAT、运行日志、完整摘要、解析后的 `case.json`、逐步 `history.csv`，以及 `data/step.../` 原生科学包。科学包含 `fields.mat`、`metadata.json` 和可选 `history.csv`，分别记录 SHA；无 LU 或 context 句柄。每次导出写新目录，保留此前数据。

```matlab
d = load(summary.dataPackage.folder + "/fields.mat");
W = d.data.fields.core.W;
rhoCanonical = d.data.scale.densityOrigin + W;
xPhysical = (d.data.mesh.core.X-d.data.scale.X_shift)*exp(-d.data.scale.logC_l);
rhoPhysical = rhoCanonical*exp(-d.data.scale.logC_omega);
```

原生核心/环带分别保存 X/Y、x/y 或 r/θ、active/ghost、W/η/q/F/Q/u1/u2 和中心化速度。重叠两份值用于审计，物理积分不能重复计数。坐标、密度、速度、规范时间和物理时钟的定义均写入包；绘图重采样不是原生误差证据。科学包用于分析；可靠续算仍从校验过的原始 checkpoint 开始。

在 macOS/Linux 上，严格预算、内存监测与冻结执行使用同一 Python 3 管理器：

```sh
python3 tools/run_explicit_experiment.py start --case examples/case.json
python3 tools/run_explicit_experiment.py status --output "$PWD/result"
```

未给命令行覆盖时从 case 的 `run` 读取输出目录、目标和总预算；命令行 `--output/--target/--hours` 是显式覆盖。冻结目录内 `native/` 是实际运行和数据输出位置。管理器分段调用，逐段校验源、选项和初值，恢复沿用原绝对截止时间。`--prepare-only` 只准备，不证明 MATLAB 或数值资格；不要对已有运行目录重复 `start`。该监督器目前没有 Windows 实现，不能用 MATLAB 入口的路径修复推导其跨平台可用性。

## 验证入口

维护清单唯一来源为 `tools/validation.json`。macOS/Linux 的 Python 启动器用 `--list` 查看检查及科学范围，用 `--suite quick/components/full` 执行维护检查；`research` 和 `history` 显式选择。MATLAB 测试函数可从 `tests` 调用。缺少外部历史数据记录 skipped，不能当成历史比较通过。默认检查不依赖本地 results；MEX 从源码构建并验证 receipt。独立制造解和误差公式仍分别实现。
