# 运行、恢复与原生数据

本文描述当前公开接口，默认使用非均匀中网格。核心形态触发、高阶转移和接口物理微分闭合分别登记其测试范围；完整长期残差复现与连续精度不从源码发布推导为已通过。需求见 [MAINLINE_REQUIREMENTS](MAINLINE_REQUIREMENTS.md)，方法见 [ARCHITECTURE](ARCHITECTURE.md)，当前与历史证据见 [VALIDATION](VALIDATION.md)。

## 入口与当前方法

```matlab
[c,rho0] = main('config');       % 只读；rho0=[] 使用登记原 Gaussian
c.run.outputDir = 'result/medium_0_5';
c.run.tauEnd = 5;
c.run.wallSeconds = 600;
[summary,s,ctx] = main(c);
```

默认是非均匀 medium、原生节点 FEM（P5/Gauss8、少量P6、最外径向P6）、linear5 五阶迎风、平滑速度系数 `.1`、RK4、CFL `.3` 和形态触发核心 remesh。默认目标200、本次预算600秒；示例的输出目录与目标须按具体实验设置。

Poisson 与对流主要内部节点完全共用；当前实现全部 `x/y/r/theta`、`X/Y` 和 shape 共用。Gauss 点只积分；raw `psi` 直接取系数，速度由同一 FE 场求导；富集时加回解析 Green 场。没有独立全局 GLL 解网格。C0、基函数次数或缓存本身不证明 C1 速度、高阶实测精度或 `O(N)`。

默认五阶迎风；若要研究中心参考格式，可在运行前显式设置 `c.method.transport=struct('reconstruction','central4','centralDissipation',.2,'smoothSpeedCoefficient',.1)`，并使用新的输出目录。主节点和中网格默认不因此改变。接口的 q/Poisson 源与 Q 使用同一套活动点物理微分闭合，详见 [架构](ARCHITECTURE.md)。

`config` 是唯一默认来源。`main` 的相对输出以项目根解析，直接 `ipm_cr.run` 以 MATLAB 当前目录解析；库调用先 `addpath('src')`。`main` 可构建本机 MEX；采样的 MEX/分块 MATLAB 路径使用同一公式，默认五阶迎风与显式 WENO 模式需要其内核。Windows 实机运行与恢复尚未验证。

```matlab
c = ipm_cr.config('case');
c.method.geometry = struct('preset','coarse');
c.run.outputDir = 'result/coarse_0_200';
c.run.tauEnd = 200;
[summary,s,ctx] = ipm_cr.run(c);
```

`coarse` 从 medium 的实际四轴精确取 `1:2:end` 并重建度量和计划：核心数组行Y/列X `129×241→65×121`，环带行r/列θ `549×193→275×97`。预设固定 `R=1e10`；自定义域/几何应显式使用生成器或完整 axes。native 非空 `stride` 仅兼容 `stride=degree`，不能生成另一套 Poisson 主体节点。

## Case 与初值

| 块 | 内容 |
|---|---|
| `initial` | 默认 `data/initial.json`，绑定 MAT SHA 和原墙参考 |
| `domain` | 共同外半径 |
| `method` | geometry、Poisson、transport、background、九参数 remesh、cfl/maxDt/timeIntegrator |
| `run` | 输出、绝对 τ 目标、本次预算、checkpoint 周期、累计 maxSteps、resume、seed |
| `stop` | 只有 `physicalTime`；空值表示不用该停止条件 |
| `save` | 原生科学包、history 和 label |

`main(c,rho0)` 接受向量化密度表达式，例如 `@(x,y)(1+tanh(x)).*exp(-y.^2)`。表达式在实际节点和原参考轴物化，重算原 I/P/C 并绑定 SHA；积分完整 W，不使用 Gaussian 余项 η。改变初值或初始网格使用新输入和新目录。

参考后端须显式选择；谱、FD、独立 GLL 组件与旧图形/案例只解释其原配置。生产 FEM 声明有限源圆盘外置零模型；外圆输运使用已选多项式闭合。入流速度和缺少物理数据的状态只记录诊断，模型是否适合本题由独立测试与实验设计事前确认。

## 步长、停止与同方法恢复

```text
h = min(maxDt, 0.95*cfl/currentRhs.rate, tauEnd-currentTau)
```

先按当前 RHS 率选初始区间 `h`。原 I/P/C 的连续样条峰跨原墙 knot 时，方程分片光滑；分段 RK4 以固定4次 secant/Newton 更新定位普通相邻 cell 的连续端点跨越，积分到 knot 后继续余段。全部探测从同一起态积分，非连续的全局峰分支交换直接按原 global maximizer 演化。

这属于原方程的时间网格构造，精度仍在独立测试中判断。没有残差阈值验收、stage-only guard、括根 fallback、折半、`minDt`、拒绝重试或旧态择优。ULP 仅识别同一坐标，不投影峰或状态。各段的完整 RK 阶段与末态 RHS 照常计算，CFL、规范和 Poisson 质量诊断被动保存；入流容差开关、在线稳态精度早停与旧四ulp末段规则保持删除。

`wallSeconds` 在完整步边界检查；初始化、单步和落盘可能越过软调用预算，严格总截止由外部监督器执行。创建输出目录 `STOP` 后在完整态停止；移除自己的 STOP 后，同一源码、初值、方法和目录可继续。实际网格从 checkpoint 冷建，不从初始 preset 替换。

| `summary.status` | 含义 |
|---|---|
| `target_reached` | 已完成态到达 τ 目标 |
| `wall_budget / maximum_steps` | 本次预算或累计步数上限 |
| `stopped / physical_time` | 用户 STOP 或跨过物理时间阈值 |
| `failed` | 计算、输入身份或存储错误；没有自动数值重试 |

文件名 `accepted_*` 表示已提交记录，不是精度接受门。源码、方法或初始 authority 改变后建立新运行；跨方法/网格晚期导入使用下一节的 seed。

## 形态触发 remesh

`ipm_cr.remesh('config')` 返回形态触发与节点重分布设置：

| 参数 | 默认 | 数学作用 |
|---|---:|---|
| `enabled` | `true` | 启用形态驱动 remesh |
| `level` | `.25` | 正 q 峰值的超水平比例 |
| `threshold` | `.20` | 相对最后实际 remesh 的方向尺度累积变化 |
| `checkEveryTau` | `.5` | 最大检查间隔；到时只检测，不强制插值 |
| `minCheckEveryTau` | `.005` | 原算法的最小观测间隔，不是积分最小步长 |
| `monitorStrength` | `2` | 当前源 profile 对固定参考目标密度的权重 |
| `blend` | `.25` | 朝目标轴移动的凸组合系数 |
| `maxNodeMotion` | `.3` | 以当前局部间距计的移动上限 |
| `candidateAxes` | `[]` | 空值使用源 profile；非空显式指定完整目标 axes |

首次完成步只校准25%正源主峰连墙分量的横纵最大范围。以后与最后一次实际换网格后保存的范围比较，任一方向累积相对变化达到20%才重分布；普通检查不能重置该参考。固定形状经过很长时间也不会重复插值。分辨率指标/profile只构造目标密度，几何阈值不属于数值精度验收门。固定参考核心防止重复累乘密度；保持点数、patch 外框、原点和环带轴。

检测调度也沿用旧算法：由最近两次形态观测估计方向变化率，用距20%阈值的剩余量安排下一次检查，间隔取 `max(.005,min(.5,.8*remaining/speed))`；形状不变时取最大间隔。该预测只安排观测，实际插值仍只由已观测的累积变化触发，不改变PDE积分步长。`lastShape / lastCheckTau / nextCheckTau` 随检查点保存。

执行顺序为有界移动 → 活动七点高阶转移 → `mesh.build`/新 context → 一次 fresh 完整 RHS → 直接继续。有 η 时转移 η，否则转移 W；原参考、尺度和四字时钟沿数据流保留。没有候选验收、反向转移检查、场量误差门或回滚；转移精度、时钟和冷恢复由独立测试验证。实际计算错误由运行器记录。

节点数和完整轴端点不变时，每个 patch 只用本片活动 donor 的最近七点多项式；新暴露活动条带继续本片端部多项式，避免把另一片独立表示注入条带。通用 `localMap`/`tensorMap` 默认禁止外推，该延拓由 remesh 显式启用；改变 patch 框或节点数的 seed 导入仍按实际覆盖转移。七点次数不保证未分辨剖面的精度或无过冲。

## 本轮复现窗口与晚期 seed

| 编号 | 设置 |
|---|---|
| T1 | medium，原 Gaussian `0→5` |
| T2 | medium，有来源的真实 τ390 态 `390→395` |
| T3 | 精确半密度 coarse，相同物理初值 `0→200` |

结果保存后独立对照已有数据，匹配实际 τ、原始问题、规范、尺度和时钟；不同方法、网格、参考底限及历史差异须显式报告。到达目标不等于复现质量或连续精度通过。

```matlab
c = ipm_cr.config('case');
c.run.outputDir = 'result/medium_390_395';
c.run.seed = struct('file','/absolute/accepted_tau390.mat', ...
    'sha256','<实际文件的64位小写SHA256>','tau',390);
c.run.tauEnd = 395;
[summary,s,ctx] = ipm_cr.run(c);
```

seed 是包含完整 checkpoint 的真实态；核对文件 SHA、实际 τ、原墙、初值 authority 与背景表示。保留尺度、步号和全部时钟，构建目标 RHS，W/q/F/Q 转移比较被动记录。该导入建立新运行，不重标时间，也不称同方法 resume 或从零演化。

## 监督器与数据

```sh
python3 tools/run_explicit_experiment.py start --case /absolute/case.json --output /absolute/run
python3 tools/run_explicit_experiment.py status --output /absolute/run
python3 tools/run_explicit_experiment.py resume --output /absolute/run
```

`start --prepare-only` 只冻结不计算。既有监督器冻结 src、mesh.json、authority、选项及必要工具；seed 复制到 `input_seed.mat` 并保留原来源及 SHA，段脚本使用冻结副本。每段与恢复核对文件身份，单线程执行，外部控制资源预算、STOP 和累计绝对 deadline；`resume` 不自动延长预算。它不根据数值误差决定是否继续。

| 文件/接口 | 内容 |
|---|---|
| `latest.json`、`accepted_0/1.mat` | SHA 绑定的原子已完成态记录 |
| `initial.mat/json` | 真实起态和 SHA，滚动记录不覆盖 |
| checkpoint、`mesh/` | 主态、真实 geometry、原墙、方法/输入、周期控制器、history、事件谱系 |
| summary、`history.csv` | 真实终点、停止原因、被动 CFL/gauge/Poisson/峰诊断和成本 |
| `data/.../fields.mat` 与 metadata | 原生 W/η/q/F/Q、速度、网格、尺度、字段定义和 SHA |
| `ipm_cr.sample(s,ctx,X,Y)` | 同坐标只读采样，独立于显示重采样 |

完整规范密度为 `densityOrigin+W`，物理密度再乘 `exp(-logC_omega)`。物理时钟高、低和两个尾字分开保留，CSV 不替代完整时钟。完整 F/Q、物理速度与去平移速度分别解释。

## 测试与旧结果

组件、制造解、阶数、稳定性、转移和冷恢复在正式实验之外调用同一生产 API；独立 `test_steady` 不进入运行停止逻辑。已有静态 Green 审计按冻结证据引用，新的计算与复现计划另行记录。

`test_remesh_order` 对三层非零解析 η，比较同一新网格上的转移态与直接解析态完整 W/q/F/Q，记录全活动点与固定公共查询误差、位置和可分辨阶。它与真实 η 往返及单事件复算共同验证本次修复，不替代修复版 `0→200` 或长期残差复现；统一32项结果留待 `repair_release_01` 完成。

`compute_only_release_04` 已完成31项 requested checks，但时间资格为 false。本次新旧端点差在 `τ≤.11` 约 `1e-15`，`.12/.34` 跨原自然样条 knot 后约 `1e-12`；光滑 `.10/.11` 仍约 `4.01` 阶。该历史失败保留；固定四次弦截分段已通过原 triplet 复测及完整 suite05，W 为3.9606阶。

旧 `b890305` 的31项统一检查和有限空间/时间筛查只属于该冻结版本。旧 T1/T3 为 owner STOP，分别停于约τ1.2878和23.297，后者3680步且保留旧往返门9次拒绝。它们不代表新窗口完成。旧 GLL/P5、Q4、谱、FD、旧均匀轨道、示例和 frozen runtime 保留原来源，不自动认证当前精度、C1 速度、低复杂度或未知物理源尾。
