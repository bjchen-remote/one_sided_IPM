# 归档：极坐标共享 SDIRK2 架构

本文件保存 2026-10-20 切换到 CFL + SSPRK3 之前的说明及当轮 SDIRK2 实测记录。以下内容属于历史架构，不能作为当前默认显式求解器的验证结果。当前入口、配置与状态见 [IPM_C.md](IPM_C.md)。

以下旧命令依赖当时的源码快照；当前 `ipm_c.config()` 和 `main_c` 已切换默认路径，不能用当前入口直接复现旧 SDIRK2 方法。原报告与失败记录留在原运行目录中。

---

# 原共享 SDIRK2 说明

`main_c.m` 默认直接推进完整 `rho`，时间算法与当前 `time_evolution/main.m` 一致：
`sdirk2_alexander_practical_v2`。它复用 `ipm.evolve.stepRk`、`advance`、
`selectTimestep` 和共享 `sdirkContract`，通过空间回调提供极坐标离散。

**新共享架构已完成本轮有界回归。** 旧 Cartesian 算术逐位一致；唯一原点
制造解的密度收敛阶为 2.005/2.002，物理钟为 1.974/1.987；真实极坐标完整 RHS
推进、实际不同径向网格的提交/回滚、两槽检查点和独立进程恢复均通过。
演化空间精度与长轨迹仍未认证。
此前 primary-v 显式基线和 expRK/Radau 实验的通过项、失败项及耗时均保留，
但不能作为这一新默认架构的性能、空间精度或长轨迹资格。

## 方程与状态

坐标为 `X=r*cos(theta)`、`Y=r*sin(theta)`，`0<=theta<=pi`。沿用原半平面方程、
固定 Cartesian 墙节点上的 I/P/C 观测与初始参考：

```text
q = D_X rho,  -Delta psi = q,  u=(-psi_Y,psi_X)
rho_tau = -(u1+cl*X+cr)*rho_X -(u2+cl*Y)*rho_Y + cw*rho
ell'=cl, omega'=cw, S'=cl*S+cr, t'=exp(omega-ell)
```

默认状态保存完整 `rho` 及 `+ipm` 的尺度/物理钟状态。唯一原点压缩后共有
`1+(nr-1)*ntheta` 个自由度；两条物理墙射线接成一条 signed wall，原点更新一次。
默认路径不采用旧 `rho=b+v` 主余项、Gaussian primary-w 或额外指数输运系统。

I/P/C 初始参考保持 `I0=0.8340921004503201`、`P0=0.6369025096521013`、
`C0=0.06283125189791006`。原 Gaussian 参数为 `A=.7,a=0,sx=.75,sy=1,w=.2`；
WENO 设置仍为 `epsilon=1e-12,alphaMode=local,inflowDensityDegree=1`。

## 共享时间架构

二阶段 SDIRK 使用 `gamma=1-1/sqrt(2)`，每阶段求完整非线性方程
`R=rho_stage-rho_known-gamma*h*F(rho_stage)=0`。原算法和控制量保留：

- 同时检查完整密度残差及 `D_X R`；两阶段、系数、容差与尺度消元使用共享合同。
- 原 ILUT 只作方向预条件。v2 的廉价缺陷修正必须使实际完整残差 merit 至少减半；
  否则使用包含 Darcy 与 gauge 的完整中央差分 Jv-GMRES，以及原 Newton 回溯规则。
- `selectTimestep` 使用原最大步长与 canonical/物理终点约束。隐式路径保留全域 CFL
  诊断，阶段失败按原规则减半重试并逐步恢复；默认最大 canonical 步长为 `.2`。
- 接受态 RHS 缓存、fresh 终点求值和接受事务沿用共享路径。缓存绑定实际几何、算子
  与 remesh 版本；它不序列化为恢复所依赖的物理状态。

默认没有旧显式基线的逐步 whole/halves 对照，也没有实验 expRK/Radau 的 M1/M2
强迫 action。这里的 `.2` 是配置上限；实际接受步长、数值误差和成本均须实测。

## 空间回调与模块

`ops.spatialOperators` 提供完整 RHS、`derivativeX`、唯一原点 `pack/unpack`、
极坐标一阶迎风 ILUT 预条件器及空间重网格接口。共享时间核无需伪造 Cartesian
`x/y/Dx`。空间部分负责完整 source、Poisson/Green、Darcy、WENO、joined wall、
原固定墙 I/P/C 与完整 `Q=D_X F`；这些离散改变仍须独立空间精度验证。

- `+ipm_c/+spatial`：极坐标空间后端及共享时间核适配。
- `+ipm/+evolve`：共享 SDIRK 合同、阶段求解、步长选择和接受/重试流程。
- `+ipm_c/+workflow/solveShared.m`：默认有界运行、输出及恢复组织。
- `+ipm_c/+output`：共享架构的检查点和输出；与旧 primary-v 检查点分开识别。
- `+ipm_c/+workflow/solveExplicit.m` 与原 `+evolve` 显式模块：legacy 兼容路径。

默认配置由 `ipm_c.config()` 返回；`ipm_c.config('explicit_legacy')` 才返回旧显式配置。
旧 primary-v 检查点不能静默转换为完整 rho 的新轨迹。已有 legacy 恢复仍按其原配置
与源码清单检查，历史目录和失败报告不覆盖。

## 执行与验证状态

```matlab
addpath('/Users/cbjin/Documents/one_sided_IPM/time_evolution');
main_c('check');                    % 默认共享架构的初始完整 RHS 检查
main_c('test');                     % 新共享组件验证入口
main_c('pilot','/absolute/new_run',1e-4);
main_c('resume','/absolute/existing_shared_run',2e-4);
```

本轮验证目录为 `runs/ipm_c_ipm_architecture_20261020_01/validation_02`，
完整结论与源码清单见 `summary.json`、`components.json`、`native/native_report.json`、
`native/restart_report.json` 和 `source_inventory.json`。65×65 网格推进两步到
`tau=2e-4`；独立 MATLAB 进程恢复后，密度、尺度和时钟与连续推进逐位相同。
重复发布同一接受态时，两槽检查点保持原指针安全，`data.mat` 与保留帧的 Q 一致。

默认321×257网格（82,241个独立自由度）单线程实测：请求 H=.2 后按原规则
减半三次，以 H=.025 接受，CFL=11.378；含失败尝试的 advance 耗时28.834秒，
接受的两阶段共79次完整 RHS 求值。该结果表明无显式 CFL 截断，但隐式阶段
收敛及失败重试仍限制吞吐；`.2` 不等于每一步都能接受。详见 `default_step.json`。
固定 H=.025 再连续推进两步到 tau=.05，分别耗时15.823/13.503秒，
计时脚本与回执见 `shared_fixed_step_cost_v2.m`、`fixed_step_cost.json`；
第一次补充计时的记录脚本因空结构赋值失败，原失败日志保留，未改求解器。
这些是单线程短程计时，不能直接推断 tau=10 的成本或相同精度下的加速比。
汇总状态见 [IPM_C_SHARED_VALIDATION.json](IPM_C_SHARED_VALIDATION.json)。
长期目标仍为完整 `Q` 低于 `1e-5` 后真实演化持续 `.5tau`；当前尚未取得新架构的
长程认证，也未证明600秒推进 canonical `Delta tau=10`。

## legacy 证据与独立实验

切换前说明完整保存在 [IPM_C_LEGACY_EXPLICIT.md](IPM_C_LEGACY_EXPLICIT.md)，
其中的执行命令、high/low 钟、五阶段显式步与 whole/halves 门都属于旧路径。

旧 `runs/ipm_c_native_qualification_20261020_03` 中，641×513 与1281×1025解析 t0
轨迹到达 `tau=1e-4` 并恢复；原直接 F 余项比较失败保留，独立因式化补充比较通过。
其中完整 F/Q 差为 `3.42e-8/4.06e-8`，7/9点 F/Q 差为 `8.20e-10/4.21e-9`。
该旧显式长验证后来在粗格 `tau=.0010883543291254748` 的接受检查点收束。

更早极坐标轨迹的空间失败仍见
`experiments/polar_tensor/results/polar_tau200_campaign_20261018_01/mission.json`：
`tau=.15019879090138588` 的完整张量加密 Q 差约 `2.11e-7`，未过 `1e-7` 门；
最后通过完整张量同场空间检查的端点为 `tau=.07072849807506046`。

[experiments/ipm_c_speed](../experiments/ipm_c_speed/README.md) 保留独立 expRK/Radau、
冻结输运及 primary-w 效率实验；它们不接入当前默认求解器。原局部时间通过、CDF空间
失败、物理钟失败及吞吐记录仍在
[RESULTS.json](design_ipm_c_10tau/RESULTS.json) 和原运行目录中，不能转作新共享架构结论。
