# 归档设计：完整冻结输运快块（当前不采用）

**历史未采用设计，停止作为当前实施计划。** 当前 `+ipm_c` 保留完整 rho、CFL + SSPRK3、直接 Poisson 及自动径向换网格；面积20%触发、径向正 CDF、高阶转移和完整审计后的提交/回滚见 [IPM_C.md](IPM_C.md)。以下快块、隐式迭代与拟议 API 不在当前默认包中实现，不能据此解释当前时间步。时间步失败减半重试已从当前显式路径移除，换网格候选回滚则保留。

以下保留当时的设计，不能作为当前实现或资格声明。文中的“显式基线”“冻结检查点”及拟议 API 均指当时的历史路径；当前包不保留这些替代时间算法的运行入口。

## 要解决的具体问题

旧角向快块没有同时处理径向、两条墙射线和原点的输运耦合。历史 seq8 的 `H=.003` 虽通过线性残差、Poisson、gauge 与径向诊断 cap，二倍 whole/halves 全域 Q 差仍为 `4.7860765e-5`，超过 `1e-7`；最大值在第一正径圈。该圈从余项速率到 Q 的观测放大约 `2.42e4`。seq22 的真实角向 WENO 权重已接近理想权重，而完整冻结输运方向导数与旧 B 的差明显更大，因此下一候选应保留径向/墙/原点耦合。这个证据不证明单一根因，也不保证新快块能使用更大 H。

完整张量空间门在旧 τ≈`.15019879` 仍失败。快块改进不能代替原生空间轨迹比较；降低线性求解成本也不能消除原点附近 Q 的空间敏感性。

## 数学合同：完整方程留在每个 stage

以主变量和增广变量 `y=(v, logC_l, logC_omega, X_shift, physicalClockIncrement)` 写完整系统 `y'=f(y)`。`f` 必须调用现有完整 `ipm_c.evolve.evaluate`，包括解析背景、完整 source、Poisson/Darcy、固定 I/P/C gauge、原 WENO 输运、尺度/位移及物理时钟速率。不能用快块替代 stage 的完整 RHS。

每个宏步仅在已接受起点冻结完整输运速度、`c_omega`、实际 WENO 权重、local-alpha 选择及原边界闭合，构造齐次线性算子 `B=diag(B_v,0)`。`B_v` 包括径向、角向、joined signed wall 和唯一原点的完整冻结输运，以及一次 `c_omega*v`。若某种边界闭合包含非齐次常数，该项留在余项，不能塞进齐次 B。

恒等分解是 `g(y)=f(y)-B*y`。实现可使用等价且较少取消的 affine stage defect：

```text
D_i = f(y_i) - f(y_n) - B*(y_i-y_n)
f(y_i) = f(y_n) + B*(y_i-y_n) + D_i
```

每个 stage 必须审计这个重组恒等式。Darcy 对状态的依赖、gauge 求解、WENO 权重变化、背景强迫、尺度与时钟非线性均保留在 `g`/`D_i` 中。因此分解不改连续 PDE 和原规范；冻结 B 影响时间离散及其误差，仍须实际验证。B 不是完整物理 Jacobian，也没有在这里获得刷新 B、变步长情形下的统一刚性四阶保证。

## 原点与墙的算子布局

- 使用 `N=1+(nr-1)*ntheta` 个自由度，`r=0` 只有一个标量；复数 pack/unpack 不得使用共轭转置。
- 内部正径节点使用径向与角向输运之和；物理墙与原点行由 joined signed-wall 方程覆盖。墙闭合不是再加一份输运。
- 输运部分装配完成后仅加一次 reaction。常量应满足 `B_v*1=c_omega*1`，不能靠事后改对角或修行和使测试通过。
- shifted solve 必须保留墙/内点对原点的完整反馈。原点行含输运耦合时，不能把原点独立写成 `b0/(1-sigma*c_omega)`。
- 所有源项、第一圈、远域与原点均保留；不启用拟合 jet、平滑、投影或局部掩码来取得通过。

## 拟议 API 与线性求解合同

```matlab
[B,freezeInfo] = ipm_c.fast.freeze(state,context,options);
[x,work,cache] = ipm_c.fast.solveShift(B,sigma,b,cache,options);
[action,work] = ipm_c.fast.forcingAction(B,coefficients,fraction,h,M,options);
[trial,obs,work] = ipm_c.evolve.stepFast(state,h,context,options);
```

`freeze` 只接受完整 fresh RHS 已审计的起点，返回 `apply/pack/unpack`、B 身份与冻结配置。`stepFast` 是未提交试算；接受、fresh 审计和 checkpoint 发布继续由统一控制器负责，显式路径保留作对照与回退。

`solveShift` 用 matrix-free GMRES 解原系统 `(I-sigma*B_v)x=b`。径向、角向 band solves 及 joined-wall 小块用于方向预条件；ADI 乘积只能充当预条件器，不能把乘积系统的解当作原 shifted 解。预条件器与缓存必须绑定当前 B；不能跨速度、权重、gauge 或网格变化复用旧因子。

每次返回前，以原自由度、原单位、未平衡的 `B.apply` 独立重算：

```text
||b - (x - sigma*B.apply(x))||_2 / max(||b||_2,realmin) <= 1e-11
```

内部 GMRES 残差、预条件残差和变换坐标残差均不能替代这个门。设置明确的 Krylov 维数、迭代数、墙钟与内存上限；到限或真残差失败就拒绝该试算。记录全部失败尝试、迭代、B products、分解/求解/残差时间及实际峰值 RSS。缓存预算包括新建因子时的暂存，不只计算最终驻留因子；不得构造全网格稠密矩阵。

`forcingAction` 才负责实际 HO5 系数所需的强迫线性 IVP/Radau 微步；shifted 解通过不等于 forcing action 准确。保留所有实/共轭 block、原 coupled residual 与 imaginary-leak 门，M 与 2M 成本及失败均计入。使用实值 B 的共轭因子复用属于另一个实现优化，不能在本阶段默认采用。

## 必须依次通过的资格

| 阶段 | 实际输入和决定性通过门 | 通过后允许的下一步 |
|---|---|---|
| 1. 小格制造 | 常量、Cartesian X/Y/二次多项式、光滑 Gaussian、实/复一般场；pack/unpack 精确；实/复线性和独立装配对照 ≤`1e-12`；常量 reaction 一次；原点/墙逐项检查；已知 shifted 解误差 ≤`1e-9`，原坐标真残差 ≤`1e-11`。 | 进入实际系数探针，尚无物理 H 资格。 |
| 2. 实际系数 | 从新包实际接受态冻结原速度、权重和边界标志；先小规模明确矩阵对照，再完整原网格 matrix-free shifted solves。覆盖实际 Radau 的全部三类 shift、half/full fraction 和 M/2M；保持同一 `1e-12` apply 一致性与 `1e-11` 真残差门；记录真实内存/迭代/全部耗时。人工截断小域只作诊断，不能签全域资格。 | 进入同起点、同终点物理试算。 |
| 3. 时间方法 | 原完整 RHS 每 stage 执行；split defect 重组通过；whole、two halves、M/2M 均通过下述原门；另以足够小步显式 N4/N8 参考比较同端点。经典非线性增广制造仍应呈四阶；这不代签刚性阶数。完整 fresh、保存回读、冷进程 restart 均通过。 | 仅允许已测父态/H/资源范围内短窗口。 |
| 4. 独立原生空间轨迹 | 两分辨率均从同一解析 t0、同一固定 Cartesian I/P/C 独立演化，精确落在同一 canonical τ；时间误差分别控制。用当前 `validatePair` 的完整细节点、7/9 重构和变换/时钟门比较，不以晚态迁移重启代替细轨迹。 | 逐窗口延长；一个端点通过不能批准 τ200 长程。 |

时间资格保留二倍完整差、Richardson 除数 1：rho 相对 `1e-7`，source `1e-8`，全域 Q `1e-7`，log-scale/相对 shift `1e-9`，物理 high+low clock `1e-12+1e-10*abs(t)`。原 stage/endpoint transport cap、Poisson `1e-5`、gauge `1e-5` 及线性内门全部保留。附加同端点完整 F 与 rates 比较各 ≤`1e-7`；任何新增快块稳定性规则须单独验证，不能直接删掉当前 cap。

空间资格保留全细节点 rho 相对 `1e-7`、source `1e-8`、F/Q/rates `1e-7`，以及原尺度/shift/物理高低钟门；7/9 参考宽度 F/Q 一致性均为 `1e-8`。记录每个全域最大值的位置和前四正径圈诊断，圈内结果不替代全域门。先完成新包当前短共同端点，再到已有精确端点 `.00490287607550439`；更晚窗口须重新取得独立空间和时间证据。

所有性能结论以同父态、同 H、同全部接受门的完整试算和实际接受 τ/总墙钟比较；必须包含拒步、失败、完整 RHS、checkpoint 与 restart 成本。单次 solve、局部 RHS 或分解计时不能写成到 τ200 的总加速比。

相关历史证据：[第一圈诊断](../experiments/polar_time_audit/FIRST_RING_DIAGNOSIS_20261018.md)、[实际系数局部失败](../experiments/polar_coupled_fast/results/actual_coeff_shift_seq000074_20261018_02/probe/report.json)、[原空间门失败与运行账本](../experiments/polar_tensor/CAMPAIGN_TAU200.md)。
