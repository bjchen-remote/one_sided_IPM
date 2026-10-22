# 极坐标 IPM：直接求解、轨道数据与绘图

阅读求解器时，先看[整体架构](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ARCHITECTURE_REPORT.md)，再按需查[算法细节](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ALGORITHM_DETAILS.md)。实验、计时和失败记录集中在[验证记录](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_VALIDATION_NOTES.md)；本页保留使用与数据接口。

当前生产采用 v4。完整重标度密度为 \(R=R_\circ+W\)，其中 \(R_\circ(\tau)=R(0,\tau)\)、\(W=R-R_\circ\)。`state.rho` 保存 `W`；标量 \(R_\circ\) 对应代码 `scale.densityOrigin` 和局部变量 `a`。输出字段 `rho` 仍表示完整 `R`。这是一种线性存储分解，不加入空间解析背景。半圆域为 `r∈[0,1e10], θ∈[0,π]`。

时间推进使用直接 CFL 与 SSPRK3，保留自动径向换网格和 Poisson 的 `decomposition(...,'lu')` 缓存。没有隐式时间阶段、GMRES/ILU 或失败时间步重试。候选网格审计与回滚发生在已接受端点。

## 每个模块求什么

|模块|已知量|求出的量|
|---|---|---|
|密度微分|中心数组 `W` 和实际网格|完整密度导数 `q=D_X W=R_X`|
|Poisson|`q` 和流函数边界值|解 `−ΔΨ=q`，得到流函数及 Darcy 速度|
|规范|当前密度梯度与速度|三个缩放／平移速率 `cl,cw,cr`|
|输运|`W`、速度和规范速率|`W_tau`、\(R_{\circ,\tau}\)|
|时间推进|上述右端与 CFL 步长|下一时刻的 \(W,R_\circ\)、尺度与时钟|
|换网格|接受态及完整源特征|同一时刻的新网格表示|
|输出与绘图|接受态及其实际网格|轨道、切片和 Cartesian 重采样结构体|

规范求根内部会反复调用输运，复用本次 Darcy 场；表中的“规范”和“输运”是职责划分，具体耦合关系见架构总览。

## 方程、状态与保常计算

计算坐标为 `X=r cosθ, Y=r sinθ`。先求

\[
-\Delta\Psi=R_X,\qquad U=(-\Psi_Y,\Psi_X),
\qquad V_X=U_1+c_lX+c_r,\quad V_Y=U_2+c_lY.
\]

令 `T_h(W;V)` 表示无反应 WENO 输运。内存中的实际演化系统为

\[
G=T_h(W;V)+c_\omega W,\qquad g_0=G(0),
\]
\[
\boxed{W_\tau=G-g_0,\qquad R_{\circ,\tau}=g_0+c_\omega R_\circ.}
\]

完整密度右端与诊断仍是

\[
F=W_\tau+R_{\circ,\tau}=R_\tau,\qquad Q=D_XW_\tau=D_XF.
\]

\(W,R_\circ\) 使用同一组 SSPRK3 阶段，因此在精确算术中与推进完整 `R` 线性等价。尺度满足 `(log Cl)'=cl, (log Cω)'=cw, S'=cl S+cr`。物理钟满足 `t'=Cω/Cl`，以高位 `physicalTime` 和低位 `physicalTimeLow` 补偿累计；不改变规范时间 `τ`、CFL 或步长。`τ=200` 不等于物理时间 `t=200`。

空间差分和原墙 `Dx` 保留原模板，但按每行局部差值 `Σ_j D_ij(f_j−f_i)` 计算，使常平台在乘加前精确消去。原点保持唯一零中心值，导数用前四圈的六次 Cartesian 拟合恢复。

当前径向 WENO 在原点另一侧补三个负半径样本，由对侧上半圆实际数据与原点拟合的奇 Y 部分构造，详见[原点延拓](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ALGORITHM_DETAILS.md#transport)。这项修复在冻结06完成之后采用；下文06的全程统计属于旧原点闭合。

专用 `weno5UnitPairMex` 同时预备单位速度 `±1`、单位度量的左右导数，再按实际速度和网格度量组装。核保留原权重归一化、幅度尺度、epsilon 和入流外推规则；左右 face 使用共同参考，ghost 小偏移分开保存。无 MEX 时保留共享参考路径。同一 RHS 的规范探测复用固定 `W` 的 WENO 导数和墙观测闭包；空间入口、规范入口与两次独立审计仍共做4次完整状态观测。

原 I/P/C 定义和参考不变：`I0=.8340921004503201, P0=.6369025096521013, C0=.06283125189791006`。规范的小型求根不属于隐式时间推进。

## 四阶 B 样条 Poisson 与自动网格

Poisson 角向用 FFT 实现正弦分解，径向未知模态用 **四阶 `k=4`、三次 `degree=3` B 样条 Galerkin** 展开，得到各自的 `A_m z_m=b_m`，其中 `A_m=K+m²M`。网格不变时缓存直接 LU；源和镜像 Green 外弧值每次更新。

Poisson 源在径向使用局部 `C¹` 三次 Hermite 重构，节点斜率来自缓存五点模板，单区间最多依赖6个节点。流函数余项的 B 样条为径向 `C²`；完整流函数加上墙源 lift 后一般仅为 `C¹`。独立全域缩放强残差默认门槛 `.02`，用于检测粗空间缺陷，不是精度证书。

默认网格为 `321×257`；冻结06实验为 `385×128`，首间距请求下限 `.06`。外半径始终为 `1e10`。完整源在 `.1/.2/.5` 峰值水平上的面积继续触发换网格；核心半径使用

\[
r_c=\max(r_{c,\min},\min(r_{.1},2r_{.5})),\qquad
h_0=\max(h_{\rm requested},r_c/[.6(N_r-1)]).
\]

其中 `r_.1,r_.5` 是完整水平集的最外径。近均匀核心最多占60%的逻辑轴，其余保留给完整源尾部，源不裁剪。若当前首间距大于目标间距，独立的 `resolutionRequested` 也能请求细化。

换网格只转移 `W`，保持 \(R_\circ\)、尺度与时钟不变；五次模板以局部差值求值。完整质量仍以 \(R_\circ+W\) 计算。质量／墙积分相对跳变门均为 `1e-3`；另要求固定墙梯度峰相对跳变不超过 `peakTransferTolerance=1e-4`。同新坐标的速度／Q 缩放差门分别为 `.05/.1`，并检查完整 Poisson 和规范。候选失败沿用几何阻尼或回滚，已接受时间步保留。成功提交且 `output.verbose=true` 时在既有日志打印前后 I/P、峰跳和强残差。

## data 结构体与版本

```matlab
addpath('/Users/cbjin/Documents/one_sided_IPM/time_evolution');
d = ipm_c.data('/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_fresh_v4_20261021_06/native');
plot(d.tau,d.c);
old = ipm_c.data(d,10,'nearest');
```

`d.tau,d.t,d.step,d.dt` 始终是全轨道向量，`d.selected` 保存当前切片的标量时间、尺度、速率及网格。每接受步的20列标量保留完整 `rho/F` 范数、导数、速度、质量、规范速率、尺度和残差。`t` 列为物理钟高位，切片与检查点另保存低位。

schema4 帧的 `rho_centered,rho_origin` 是权威 \(W,R_\circ\)；`rho_tau_centered,rho_origin_tau` 是对应变化率。兼容 `rho,rho_tau/F` 是完整 `R,F`，`Q/rho_xt` 是 canonical `D_XF`，不是物理混合时间导数。另保存 `rho_X,rho_Y,psi,u1,u2,vx,vy`；`u1/u2` 是 Darcy 速度，`vx/vy` 含缩放和平移，`transportOrigin` 保存 `g0`。

实际网格在 `d.grid` 及 `d.r,d.theta,d.X,d.Y`。`data.mat` 保存轨道与切片目录，大场位于 `data_frames`，相同网格在 `grids` 只存一次。相对引用支持整体搬移读取。新检查点要求严格 v4 kind，旧检查点只能经显式迁移。

已停止的04试算保留父03全部历史标量、旧场帧和实际网格。混合 schema3/schema4 由 `meta.mixedFrameSchemas` 声明；读取器按每帧实际版本解释，选旧帧时清除中心状态字段，保留父检查点和迁移来源。正式06从初态记录 schema4；旧帧不会因读取器兼容而获得新的存储精度。

## 指定方形区域的绘图接口

```matlab
opts = struct('fields',{{'rho','rho_X','psi','u1','u2'}}, ...
              'coordinates','native','method','cubic');
cartesianSlice = ipm_c.cartesian(d,[-2 2 0 4],401,opts);
graphics = ipm_c.plot(cartesianSlice,'rho_X');
exportgraphics(graphics.figure,'rho_X.pdf','ContentType','vector');
```

区域可为任意 `[xmin xmax ymin ymax]`，分辨率可用 `[nx ny]`。默认 P3 Hermite/Bernstein 有限元使用实际三角网，原点为单顶点；为 `C⁰` 连续且可能超调。`method='linear'` 是 P1，域外均为 NaN。schema4 的密度先插 `W` 再加 \(R_\circ\)，变化率先插 `W_tau` 再加 \(R_{\circ,\tau}\)。

`coordinates='physical'` 转换坐标和密度、导数、Darcy 速度及流函数单位。`rho_tau/F,Q/rho_xt,vx/vy` 保留 canonical 定义。显式请求 `rho_t` 时，schema4 使用稳定公式

\[
\rho_t\approx\frac{C_l}{C_\omega^2}
[W_\tau+g_0-c_\omega W+(c_lX+c_r)R_X+c_lYR_Y],
\]

避免先构造大完整量再相减；旧 schema3 仍用原完整场公式。绘图插值与 Poisson 的径向 B 样条是两个独立算子。

## 验收、效率与冻结06运行记录

原点延拓修复前的[真实短轨及恢复验收](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/qualified_audit_20261021_022938_410/result.json)已通过：实际10步原版／当前轨道对照、主态检查点逐位读回、独立冷恢复、混合选帧与物理绘图。相同规范热种子下重建 RHS 和速率差为零；冷恢复下一步主态归一化差约 `3.06e-16`，最大派生差约 `2.65e-12`。

同一真实10步短测约 `.485 → .329秒/步`，约1.47倍；只计推进，不含恢复及输出，短测可能受 JIT/缓存热身影响。独立 [WENO预备测试](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/validate_transport_unit_candidate_results.json)约7.37倍，[固定状态RHS缓存测试](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/measurement_cache_validation_20261021_022108252/summary.json)约1.24倍；不能相乘当作整步加速。常数 `1e100` 平台及仿射核测试通过，与原输运最大差约 `3.12e-14`。

04试算已正常停止于 `τ=33.431256616977`；其相对初始实际墙峰约 `.0849464%` 的累计偏差主要来自更早第22→23次换网格，不能归因于04的 v4 表示迁移。05仅准备，未执行。新的[峰跳回归](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/peak_gate_validation.json)在原 `τ≈8.29638` 接受态拒绝 `.00080200768` 的完整提案，沿用 `1/16` 几何阻尼后以 `6.92742364e-5` 峰跳通过全部门槛；初始网格回归通过。

正式06已从同一 Gaussian 初值 `τ=0` 算到 **`τ=200`**，采用冻结 v4、原规范参考及峰转移门。任务与数值报告均为 `target_reached`；[最终只读审计](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/final_audit_ipm_c_tau200_fresh_v4_20261021_06.json)为 `passed`、`requestedTau200Complete=true`，0条问题。实际结果为：

|项目|结果|
|---|---|
|轨道与输出|30568个接受步、30569行标量、411帧、26份实际网格、31次换网格|
|最终物理时间|`t=5.319851681889767`，另存低位 `2.80049424504877e-16`|
|保存切片最大 P/I 相对漂移|`3.30423721717e-9 / 3.29872651417e-9`|
|保存切片最大质心恢复律误差|`4.54737511535e-10`，比较 \(C_0+(C_{\rm init}-C_0)e^{-.5\tau}\)|
|保存切片全域缩放强残差|`.00193862030102–.00978460158484`，门槛 `.02`|
|成功换网格最大墙峰相对跳变|`2.48575776201e-9`，门槛 `1e-4`|
|最终 \(\lVert\mathbf Q\rVert_\infty\)|`7.43935891907e-6`|
|正式任务墙钟|7772.588015秒，即2小时9分32.588秒；含准备、启动、恢复、输出及分段空隙|

P/I 漂移相对初始实际值 `.6369024213358077/.8340921023536607`；C 恢复律使用初始实际质心 `.0628312757084418` 和原固定参考。统计覆盖全部411个保存切片，不是每个 RK 阶段的全程上界。逐步 Poisson backward error／规范残差最大值另为 `3.99591326346e-16/9.91210510480e-11`。完整场有限、中心原点精确为零、SHA及目录一致性等审计均通过。

第7段保留退出码 `−15`：归档许可初始化栈、空日志和前后相同检查点 SHA 证明该段没有数值推进，第8段从该检查点恢复。后续分段预算由300秒改为900秒；数值架构、CFL 和原北京时间09:13:35截止不变。13段累计墙钟 `7727.657199秒`（含此启动故障），约 `.252802秒/接受步`，不是纯求解核耗时。正式任务在2026年10月21日北京时间04:55:31结束。

数据：[实验汇总](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/campaign_result.json)、[物理钟检查](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/physical_clock_repeat_audit.json)、[06清单](/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_fresh_v4_20261021_06/mission.json)、[终点报告](/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_fresh_v4_20261021_06/native/report.json)、[逐步标量](/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_fresh_v4_20261021_06/native/monitor.csv)、[data目录](/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_fresh_v4_20261021_06/native/data.mat)、[第7段恢复证据](/Users/cbjin/Documents/one_sided_IPM/output/ipm_c_repair_20261021_01/segment7_launch_recovery.json)。

这是冻结06、旧原点闭合的385×128粗网格含时实验；完成目标和通过硬门不构成空间收敛、长期解精度或一般流场全局四阶证明，也不能直接认证当前原点延拓的全程轨道。更早异常轨道在 `τ≈4.93555` 停止，其关键观测、原数据和后续修复证据均保留在[验证记录](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_VALIDATION_NOTES.md)。
