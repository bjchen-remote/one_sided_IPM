# IPM：Cartesian 与极坐标时间演化

Cartesian `main.m` 是完整程序，不依赖本机 runs 目录。计算盒为 `[-1e10,1e10]×[0,1e10]`；默认513×321网格。R_X 初值为原点附近快速衰减的 Gaussian 混合，R 是其解析原函数，左侧密度为零。

`main.m` 中可修改 Gaussian 的幅度、中心、X/Y宽度及偏斜权重；权重为0时是单个 Gaussian。程序保持自身初值的壁面 `∫[-1,1] R_X dX`，使用完整壁面峰值和重心规范；二维 R_X level-set 面积相对上次成功换网格变化20%时 remesh。

Cartesian 时间后端为两阶段二阶 L 稳定 SDIRK，默认最大步长0.2。目标仍是完整二维 `max(abs(R_X_tau),[],'all')<1e-5`，随后真实演化持续0.5τ。**该长期目标尚未达到；程序中的停止判定没有放宽。** 稳态加速器仍是包外实验。

极坐标入口 `main_c.m` 推进完整 `rho`，采用三阶段 SSPRK3 和直接 CFL，保留自动径向换网格。Poisson 使用四阶 k=4（三次）径向 B 样条 Galerkin、角向 FFT 与缓存直接 LU。外边界 `r_end=1e10`，`uniform_core` 首间距 `.03`，原点为六次局部 Cartesian 拟合。默认321×257，探索实验385×128。当前 B 样条版本的验证须与此前有限差分结果区分。

`ipm_c.data` 读取逐步标量轨道并选择实际网格上的场切片，`ipm_c.cartesian` 在指定区域重采样，`ipm_c.plot` 绘图。执行器 `run_polar_experiment.py` 冻结源码、分段检查点续算，目标为 canonical `τ=200`。完整接口、数据定义及精度范围见 [IPM_C.md](IPM_C.md)。

Cartesian 路径的本机实际验证包括：Gaussian启动、旧检查点兼容、二阶时间精化、默认步长下的纯MATLAB推进/重放、35%面积触发及数据输出。服务器R2022b尚未执行，由用户对接。具体范围和失败记录见 `VALIDATION_STATUS.json`、`INSTALLATION26.json` 与 `validation_evidence/`。

Cartesian 保存的运行自动生成 `data.mat` 和 `data_frames/`；含 `residual_inf`、复合输运速度 `vx/vy/vmag`。读取、选时和打包见 `DATA_OUTPUT.md`。

```matlab
data = ipm.output.loadData('/absolute/run/data.mat');
data = ipm.output.selectTime(data,3.5);
mesh(data.x,data.y,data.rho_xt)
```

Cartesian 本地维护目录使用 `main('check')`、`main('smoke',newFolder)`、`main('run',newFolder)`；发布目录的 main 是服务器脚本。恢复保留检查点原初值、规范和数值配置；恢复旧初值不会把它改成 Gaussian。新初值须新开运行。

直接复制整个 `+ipm` 和服务器 `main.m`，无需zip；服务器命令见 `README_SERVER_ZH.md`。
