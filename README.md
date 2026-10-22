# 单侧 IPM 程序

源码跟踪范围、原始初值、构建与复现边界见 [REPOSITORY.md](REPOSITORY.md)。张量核心与极坐标环带的混合网格实现位于 [test/](test/README.md)，其 [算法报告](test/ALGORITHM_REPORT.tex)和[架构说明](test/ARCHITECTURE.md)随源码维护。AI 脚手架、本地实验输出和编译产物不纳入 Git。

极坐标数据与直接绘图示例放在 [plotting/](plotting/README.md)。在 MATLAB 中 `addpath('plotting'); data=load_polar_data();` 后，可直接使用 `imagesc(data.x,data.y,data.rho)` 和 `plot(data.t,data.c_r)`。

2026-10-21，极坐标入口 `time_evolution/main_c.m` 使用完整 `rho`、三阶段 SSPRK3 与直接 CFL，保留自动换网格及 Poisson 分解缓存。Poisson 径向采用四阶（k=4、degree=3）B 样条 Galerkin，角向采用 FFT；网格外边界为 `r_end=1e10`，近均匀核心首间距 `.03`。精度默认配置321×257，长时间探索配置385×128。新 `ipm_c.data/cartesian/plot` 接口提供逐步轨道、实际网格及场切片，支持指定方形区域的有限元样条重采样。求解器讲解分为[整体架构](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ARCHITECTURE_REPORT.md)、[算法细节](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ALGORITHM_DETAILS.md)和[验证记录](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_VALIDATION_NOTES.md)；使用方法见 [IPM_C.md](time_evolution/IPM_C.md)；此前有限差分 Poisson 的 ABBA 计时与精度记录归属历史版本。

2026-10-14 新增原生 MATLAB 数据输出：维护程序的数值 source20 保持原算法，保存结果时自动产出 `data.mat` 与不可变原生帧目录。可用 `ipm.output.loadData`、`plot(data.tau,data.cr)`、`ipm.output.selectTime` 和 `mesh(data.x,data.y,data.rho_xt)`。协议与真实验证见 [DATA_OUTPUT.md](time_evolution/DATA_OUTPUT.md)。当前 source24 的独立输出副本位于 `runs/native_data_runtime_20261014/source24`，原冻结数值源码未修改；验证与协调状态另行记录。

Cartesian 时间演化维护源在 `time_evolution/+ipm`，默认1e10盒、513×321节点，时间入口为两阶段二阶 `sdirk2_alexander_practical_v2`。保持原点附近壁面R_X积分，以平方梯度重心定位，采用二维level-set面积20%触发的光滑CDF重网格。旧五步RK、三阶时间测试（3.068）和规范19版记录属于历史版本；具体源码版本、实测范围和失败记录见 `time_evolution/ARCHITECTURE.md` 与 `VALIDATION_STATUS.json`。完整域1e-5及真实保持0.5τ目标尚未达到；外部Newton/PTC加速仍是独立实验。部署历史保存在 `releases/one_sided_ipm_wall_integral_20261007`，下文solver历史数据不代表当前版本测试。

Cartesian 独立包使用 `time_evolution/main.m` 或完整发布目录的脚本；极坐标使用 `time_evolution/main_c.m`。根目录的 `main.m` 与 `solver/` 是历史原版，下面使用方式及数值结果均属于该历史路径。运行生成的结果放在 `runs/`。

## 历史原版使用方式

在 MATLAB 中将当前文件夹切换到本目录。打开 `main.m`，在文件开头修改网格节点数、计算区域、初值和 C 的三个规范条件。默认值是 401 × 321 个节点，即 400 × 320 个网格单元。`C_s=c_l/(c_l-c_omega)` 由求解结果计算，不在输入中指定。

```matlab
main('check')                            % 检查网格、初值和参数，不演化
main('smoke')                            % 演化一步，检查完整调用链
main                                   % 按 main.m 中的设置正式运行
main('run','/absolute/path/to/new_run') % 指定新的输出目录
```

`check` 只验证初始网格、初值和参数组合；完整演化是否成功，以正式运行的 `run.log` 为准。
运行时每个记录时刻只输出一行：规范时间 `tau`、物理时间 `t`、三个规范速率 `c_l`、`c_omega`、`c_r`、导出量 `C_s`、全域/核心区的 `max|∂τ R_X|`、物理质量漂移、X/Y 最大相邻网格间距比以及已执行的重网格次数 `remesh`。审核的漂移列改为峰值/曲率相对初值的漂移。
`drift(mass)` 比较的是当前物理覆盖范围内的质量与初值；物理覆盖范围随缩放和平移变化，且外边界开放，因此它本身不是封闭域的守恒误差。

新运行默认开启形态自适应重网格。求解器按壁面峰值和不同高度截面的分辨率判断是否需要重新分布 X/Y 网格，完成插值并校验后才接受新网格；`remesh` 增加表示实际发生了重网格。旧检查点续算沿用其原有的重网格设置。

网格对照在 `tau=0.7` 给出：321 × 257 到 401 × 321 时，三个规范速率和全域范数的差异约为 0.7%～4%；401 × 321 到 481 × 385 时约为 0.05%～0.8%。但 `tau=0` 的 `c_l` 在 321、401、481、561、641 个 X 节点上分别约为 0.785、0.877、1.089、1.878、0.636，初始规范速率尚未呈现网格收敛。默认网格是运行设置，不是精度保证；长时间结果还需逐段作细网格对照。

正式运行会打开一幅实时图：上半部显示边界中心附近的 `R_X=∂_X R`，下半部显示壁面 `Y=0` 的同一剖面。视窗随核心宽度调整，每次终端记录时更新，也会跟随重网格。关闭图窗不会中断求解。

未指定输出目录时，程序自动创建 `runs/run_时间戳/`。每次运行都要求新目录，不覆盖旧结果。运行结果包括 `settings_used.mat`、`result_or_target.mat`、`monitor_full_absolute_rate.csv`、`run.log` 和 `checkpoints/checkpoint_*.mat`。

使用已经生成的检查点继续运行或独立审核：

```matlab
main('resume','/absolute/path/to/checkpoint.mat','/absolute/path/to/new_segment',14)
main('audit','/absolute/path/to/checkpoint.mat','/absolute/path/to/new_audit')
```

`resume` 的最后一个参数是新的规范时间终点。续算沿用检查点中的方程、网格、初值和 C 条件；后来对 `main.m` 的修改不会改变旧轨迹。`audit` 针对固定网格、单侧各向同性的检查点。

最新 τ=20 数据的奇偶分解、局部稳态候选、网格/区域复核及解析构造见 [稳态试算结果](runs/steady_20260928_tau20/RESULTS.md)。文末包含二维倾斜结构、全域应变、内外展开，以及由峰脊和左右半高宽构成的反向几何坐标；分析代码统一在 `solver/steady_analysis/`。候选以 MATLAB 文件 `steady_candidates.mat` 保存在结果目录，当前尚未通过半平面稳态的收敛验证，也尚未得到极限行波的精确表达式。

## 程序架构

当前维护版本：

```text
time_evolution/
  ├─ main.m                   Cartesian 开发函数入口
  ├─ main_c.m                 极坐标 CFL + SSPRK3 入口
  ├─ +ipm_c/                 极坐标、直接 Poisson、显式推进与自动径向换网格
  ├─ server/main_server.m     服务器脚本模板
  ├─ +ipm/                   完整时间演化包
  │    ├─ solve.m            正式收敛观察、循环、保存
  │    ├─ +workflow/         配置映射、新运行和续算
  │    ├─ +config/           配置校验
  │    ├─ +mesh/             网格、离散算子和 Schur 因子
  │    ├─ +field/            初值、Green、Poisson、速度和 WENO
  │    ├─ +evolve/           规范方程、共享SDIRK、步长和重网格调度
  │    ├─ +remesh/           二维几何、候选、转移及验收
  │    ├─ +diagnostics/      观测和质量判据
  │    └─ +output/           日志、检查点和结果
  ├─ tests/                  历史组件回归和实用规范阶数／重网格／长测
  └─ ARCHITECTURE.md         按源码核对的架构及未验证边界
releases/one_sided_ipm_wall_integral_20261007/
  └─ 完整部署副本：+ipm、main.m、tests、说明及校验清单
```

Cartesian 服务器新运行：`main → workflow.run_fresh → solve → initialize/advance`。该维护入口为半平面、各向同性、WENO5与 `sdirk2_alexander_practical_v2`；其此前 SSPRK3 的记录属于历史版本。极坐标 `main_c` 默认采用独立的 CFL + SSPRK3 显式路径；该路径本轮组件、原生推进、实际换网格事务与冷启动恢复已通过。

历史原版仍在根 `main.m` 和 `solver/`，其独立审核和稳态分析在 `solver/audit_method/`、`solver/steady_analysis/`。旧版曲率轨迹和历史数据不能作为当前积分规范版本的测试结果。

## 仿射自适应 r 研究支线

`experiments/affine_tensor` 独立实现 `X=xi+s*r,Y=r` 的平行四边形张量网格，重点优化固定节点数下的 r 布点，同时保留完整 CFL 与原壁面 I/P/C 规范。从原解析初值 tau=0 出发，目标为 tau=200；当前案例由 [CURRENT_RUN.json](experiments/affine_tensor/CURRENT_RUN.json) 指向，执行、误差检查和恢复说明见 [RUN_ZH.txt](experiments/affine_tensor/RUN_ZH.txt)。进度以当前段报告和账本为准，组件验证不代表已完成长轨迹或空间收敛认证。
