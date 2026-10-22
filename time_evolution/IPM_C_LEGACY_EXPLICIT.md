# legacy：旧 primary-v 显式基线档案

本文保存共享 SDIRK 架构切换前的说明与当时证据；所有“当前”、入口和执行命令均属于旧版本快照。现默认接口见 [IPM_C.md](IPM_C.md)。不要将下列旧驱动命令用于新默认配置；旧表示需显式选择 `ipm_c.config('explicit_legacy')`，旧检查点不得静默转换为完整 rho 新轨迹。

---

# +ipm_c：极坐标求解器与验证状态

本包求解原半平面 IPM 真实含时初值问题，使用现有固定 I/P/C 规范。
入口为 `main_c.m`。它是可运行的显式基线；新包尚未完成长时间精度认证。
原 `+ipm` 与极坐标实验的生产轨迹均未修改。

本轮已实际通过 10 项 MATLAB 组件/时间/检查点测试及 3 项执行器进程故障测试。
641×513 与 1281×1025 两条独立解析 t0 轨迹均到达 τ=1e-4（分别3、11个接受步），
分别通过新 MATLAB 进程恢复。原直接 F 余项比较在已知远壁 Gaussian 层失败，
失败记录保留；独立补充比较按旧审计的精确因式化通过全部原门。
全域 F/Q 差为 3.42e-8 / 4.06e-8，7/9 点 F/Q 差为 8.20e-10 / 4.21e-9。
数据位于 `runs/ipm_c_native_qualification_20261020_03`。

已启动最长12小时的串行长验证，下一共同端点为 τ=.00490287607550439；
实际状态以该目录的 `mission.json`、各段执行回执、当前检查点与真实进程为准。
τ200 尚未达到；显式基线尚未解决长期效率。

用户随后要求600秒内推进10单位仿真时间、径向沿用`+ipm`布点思路且角向均匀。
据此，本次慢速长验证已在粗格 τ=.0010883543291254748 的完整接受检查点收束，
资源转向 `experiments/ipm_c_speed` 的独立网格与求解器性能候选。
新目标尚未达成，计时/精度合同及可行性记录在 `design_ipm_c_10tau`。

## 为什么另建包

截至 2026-10-20，原极坐标主轨迹停在 τ=0.15019879090138588。
该状态的完整张量加密 Q 差约 2.11e-7，超过 1e-7 门限；
7/9 点重构不确定性约 1.76e-7，超过 1e-8 门限。
最后通过完整张量同场空间检查的是 τ=0.07072849807506046。
独立径向加密轨迹仅到 τ=0.0010150299189975354。
原记录见 `experiments/polar_tensor/results/polar_tau200_campaign_20261018_01/mission.json`。

另一障碍是近原点 Q 的敏感性：某 H=.003 试算的密度差仅 3.78e-12，
但二倍 Q 差达到 4.79e-5。CFL、线性残差和 Poisson 残差通过不能代替 Q 门。
仅改包结构也不能消除这个数学和精度问题。

## 数值合同

坐标为 X=r cosθ、Y=r sinθ，0≤θ≤π；原点为同一个物理点。

```
q = rho_X,  -Delta psi = q,  u=(-psi_Y,psi_X)
rho_tau = -(u1+cl*X+cr)*rho_X -(u2+cl*Y)*rho_Y + cw*rho
ell'=cl, omega'=cw, S'=cl*S+cr, t'=exp(omega-ell)
```

保存并直接推进主余项 v，rho=b+v；b 是原 Gaussian 初值经实际尺度/位移变换后的
解析背景。source、F、Q 和 I/P/C 审计都包含完整背景。不得从两个已舍入的总场相减
来恢复待推进的 v 或 v_tau。物理钟携带 high/low 两分量。

I 为原固定墙节点上 [-1,1] 的截断带符号积分，P 为完整墙区间的自然样条全局最大值，
C 为固定窗口 q² 重心。初始参考跨网格共享：
I0=0.8340921004503201，P0=0.6369025096521013，C0=0.06283125189791006。
Gaussian 参数为 A=.7、a=0、sx=.75、sy=1、w=.2。
原始 MAT 的输运设置是 epsilon=1e-12、alphaMode=local、inflowDensityDegree=1。

每阶段和接受端点均重新求解完整物理 RHS 与规范；Poisson 物理模态残差及完整
I/P/C 残差门均为 1e-5。时间基线为原 HO expRK4s5 的零快块极限，经典四阶、五阶段。
每整步必须通过两半步对照和比较后的独立新鲜 RHS 复核，只提交整步。

| 时间误差量（观测差乘 2，无 Richardson 除数） | 门限 |
|---|---:|
| 全节点 rho 差 / max(1,abs(rho)) | 1e-7 |
| 全节点 source 差 | 1e-8 |
| 全节点 Q 差 | 1e-7 |
| 两个对数尺度 | 1e-9 |
| 位移差 / max(1,abs(S)) | 1e-9 |
| 物理钟 high+low 差 | 1e-12 + 1e-10 abs(t) |

完整 Q 的稳态目标保留 1e-6 并持续 .5τ；这不是启动接受条件。
当前短程研究报告始终标记长期认证和稳态完成为 false。

## 模块

- `+mesh`：sinh 径向网格、唯一原点、物理导数与几何。
- `+field`：Poisson/Green、Darcy 速度、WENO 输运、解析背景和完整 RHS。
- `+gauge`：原固定墙观测、非线性规范求解、初始参考权威。
- `+evolve`：五阶段步、物理钟、完整误差比较、实际 RHS 调用计数。
- `+output`：源码和二进制清单、流式 SHA、两槽原子保存/真实回读/恢复。
- `solve`：有界段、拒步、接受事务、标量日志；不在状态中保存整套旧 flow 和临时缓存。

10 个数学/初值模块从现有极坐标基准只作命名空间迁移，逐文件来源见
`ipm_c_module_provenance.json`。共享 `+ipm` 的 FD、WENO、初值和三个墙观测函数，
恢复清单同时绑定这些依赖及 MEX/C/CPP。这仍要求整个 `time_evolution` 目录在路径上。

## 执行

```matlab
addpath('/Users/cbjin/Documents/one_sided_IPM/time_evolution');
main_c('check');
main_c('test');
main_c('pilot','/absolute/path/to/new_run',1e-4);
main_c('resume','/absolute/path/to/existing_run',2e-4);
```

从任意当前目录执行两条独立原生轨迹的首个共同端点验证：

```sh
python3 /Users/cbjin/Documents/one_sided_IPM/time_evolution/run_ipm_c.py qualify \
  --output /absolute/path/to/new_qualification
```

驱动串行运行组件测试、641×513 与 1281×1025 原解析 t0 轨迹至精确 τ=1e-4、
分别在新 MATLAB 进程中恢复，最后比较完整细网格字段及 7/9 点重构不确定性。
每段为 300 秒软预算，完整验证批可跨过该边界；外层每进程 480 秒硬预算、单线程、
采样 RSS 上限 2.25 GiB，失败保留全部结果。
只在正常段预算耗尽时同配置承接，最多三段；数值失败不自动反复重试。
不会更改另一个聊天的 STOP、checkpoint 或自动化。

首个共同端点实际通过后，可串行执行下一段较长验证；指定秒数约束总计算预算：

```sh
python3 /Users/cbjin/Documents/one_sided_IPM/time_evolution/run_ipm_c.py continue \
  --output /absolute/path/to/existing_qualification \
  --target .00490287607550439 --seconds 43200
```

它保留首端点快照，再推进两条原生轨迹并重新对照；检查资源不满足时直接结束，
不建立等待器。根目录 `STOP` 可让它在当前有界段结束后停止。

## 长时间模拟的验收顺序

1. 全部组件、时间阶数、high/low 钟、拒绝门、检查点与真实分进程恢复通过。
2. 两条独立原生轨迹到同一 τ=1e-4，完整全域字段与重构不确定性过门。
3. 推进到 τ=0.00490287607550439 的共同端点，保留演化空间检查和物理 high/low 钟。
4. 依次扩展到 .15、1、10、200，每窗都重新验证演化后的解析能力和误差。

首版 `solve` 对 τ>.005 返回 `SpatialQualificationRequired`，防止把起步局部结果当成长程资格。
后续正式长程控制器需要接入独立空间窗证据及失败停止，然后才开放对应窗口。
这一限制不是用户审批步骤；缺失的是数值证据。

角向/径向耦合快块留作独立实现和验证：apply、shift solve、forcing action、完整时间步
分别验收，加入完整 residual 和 factor/solve/apply/RHS 分项计时。
当前显式基线未解决长期效率，不能据此承诺 τ200 已可经济完成。
具体快块接口、求解算法及验收门见 [IPM_C_FAST_PLAN.md](IPM_C_FAST_PLAN.md)。

最终实测结果在本次 `runs/ipm_c_*` 目录和 `IPM_C_VALIDATION.json` 中更新。

## 600 秒推进 10 个时间单位的新要求

独立效率候选位于 [experiments/ipm_c_speed](../experiments/ipm_c_speed/README.md)，
采用 `+ipm` 风格的径向正 CDF 和均匀角轴；未改本包冻结数值模块。
完整耦合快块的局部时间检验已在原 t0、321×257、H=.0125 上通过，
紧凑 CDF 空间检验与整体吞吐仍不通过，不能据此开放本包的长程限制。
最新逐项结果和原报告 SHA 见 [RESULTS.json](design_ipm_c_10tau/RESULTS.json)。
