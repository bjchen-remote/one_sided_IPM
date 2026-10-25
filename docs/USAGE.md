# 运行 case 与论文数据

配置的默认来源为 `ipm_cr.config`。`ipm_cr.config()` 保留数值方法接口；`ipm_cr.config('case')` 生成供用户编辑的 case，`ipm_cr.run` 接受该结构或 JSON 文件。`examples/case.json` 是明确选择中层网格和谱次数的 τ10 示例，不是第二份默认配置。

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

物理初值与 `method.background` 是两件事：后者只选择同一密度的解析表示。默认使用与初值完全匹配的解析背景，初始 η 严格为零。

第一版支持 Gaussian primitive 参数族，新参数生成独立 authority，保留固定原测量轴、零时钟、单位尺度，并重新定义该问题的 I0/P0/C0：

```matlab
addpath('tools');
[initial,receipt] = make_initial('wide_x',struct('widthX',1.0,'skewWeight',.15));
c = ipm_cr.config('case');
c.initial = initial;
```

可调 `amplitude/centerX/widthX/widthY/skewWeight`。生成器拒绝覆盖已有 case。`data/cases/` 默认不入 Git；有意共享某个问题时须同时保留其 JSON、MAT 和 SHA。manifest 与 MAT 路径均相对于项目根，禁止外逃。任意函数句柄或随网格变化的 η₀ 数组尚未开放；它们需要额外的初值、尾部与固定参考契约。新初值不能继承旧轨道或旧精度结论。

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

严格预算、内存监测与冻结执行使用同一 Python 管理器：

```sh
python3 tools/run_explicit_experiment.py start --case examples/case.json
python3 tools/run_explicit_experiment.py status --output "$PWD/result"
```

未给命令行覆盖时从 case 的 `run` 读取输出目录、目标和总预算；命令行 `--output/--target/--hours` 是显式覆盖。冻结目录内 `native/` 是实际运行和数据输出位置。管理器分段调用，逐段校验源、选项和初值，恢复沿用原绝对截止时间。`--prepare-only` 只准备，不证明 MATLAB 或数值资格；不要对已有运行目录重复 `start`。

## 验证入口

维护清单唯一来源为 `tools/validation.json`。`--list` 查看检查及科学范围，`--suite quick/components/full` 执行相应维护检查；`research` 和 `history` 显式选择。缺少外部历史数据记录 skipped，不能当成历史比较通过。默认检查不依赖本地 results；MEX 从源码构建并验证 receipt。独立制造解和误差公式仍分别实现。
