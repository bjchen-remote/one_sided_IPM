# CHESS 的 MATLAB 作业

旧服务器发布包 `one_sided_ipm_wall_integral_20261007` 从 CHESS 生成的作业目录启动 MATLAB 脚本，从固定的绝对路径加载求解器。这个发布包沿用同一提交流程，每个作业把实际源码、输入、日志和检查点留在自己的目录。

## 部署与提交

1. 将完整的 `ipm_cr_chess_configurable` 目录上传到 `/data/user/hd58131/ipm/`。如果使用其他位置，把 `server/main_chess.m` 顶部的 `releaseRoot` 改成实际绝对路径，也可设置 `IPM_CHESS_RELEASE_ROOT`。
2. 在 CHESS UI 创建 **1 节点、1 MATLAB 进程**的作业，选择所需的 **CPU 数**。在 `server/main_chess.m` 顶部设置参数，将该脚本提交给 CHESS。作业时限留出 MEX 编译、初始化和最后保存的时间。
3. 在作业目录查看 `native/run.log`、`native/latest.json` 和 `native/summary.json`。最终实际网格点数和 MATLAB 线程设置也记录在 `chess_run.json`。

当前发布入口默认 `k=8`、`threads=16`、`targetTau=400`、`wallSeconds=64*3600`。对应四方向点数为 X1=513、X2=257、r=513、theta=361，包含 ghost 与角向端点。`wallSeconds` 是每次进程运行的计算墙钟预算，独立于 CHESS 调度时限。提交前按实际任务修改；使用下方四方向手动计数时先将 `k=[]`。

发布版已集成核心区缓存稀疏 Cholesky 求解对象。原网格、离散方程和固定残差修正保持不变，普通 RK4 步仍为四次完整 Poisson。本机 k=8、4 线程重复短测的完整步中位数由 2.17 s 降为 1.17 s，初始化由约 21–24 s 降为 9.6–10.1 s。制造解、重网格和冷恢复测试通过；服务器 R2022b、16 CPU 的实际耗时需上机测量。优化会改变浮点运算顺序，晚期极小 Q 残差不保证逐位相同。

完整包包含 `src/`、`data/`、`examples/`、`tools/`、`server/` 和 `main.m`，版本见 `SOURCE_COMMIT`，逐文件校验见 `SHA256SUMS`。采样 MEX 在目标平台编译；本机二进制、实验结果及检查点不放入发布包。

上传新包后用新的 CHESS 作业目录启动，才会冻结并使用新版。已有作业从其自身 `runtime/` 恢复，上传发布包不会改变它使用的实现。沿用旧版本恢复已有作业时，`IPM_CHESS_JOB_DIR` 指向原作业绝对路径，脚本网格设置须与该目录的 `case.json` 一致；原检查点的来源校验继续生效。

例如，这组参数从非均匀中网格规则生成自定义阶梯网格。

```matlab
caseName = 'examples/bspline4.json';
k = [];        % Clear the default k before setting counts separately.
nx = 225;       % X1，含 ghost，使用奇数保留原点
ny = 125;       % X2，含 ghost
nr = 465;       % 实际环带径向点数，含内侧 ghost
ntheta = 161;   % 角向点数，含 0 和 pi
threads = 16;   % 与 CHESS UI 选定的 CPU 数一致
targetTau = 76;
wallSeconds = 8*3600;
```

每个网格数可单独填写。留空 `[]` 时由所选算例规则决定。显式设置的计数在自适应重网格后继续生效。设置 `nr` 时，径向网格保留内边界等距带，在外侧用平滑过渡达到指定点数。初次建网格和重网格都会检查相邻网格间距比。超过目标时发出警告，作业继续使用有效的网格。坐标、度量和 ghost 插值覆盖仍需满足构造条件。

使用中网格的其余默认设置时，可以这样配置 256 个径向点。

```matlab
nr = 256;
radialRatioLimit = 1.2;
```

`radialRatioLimit` 设置径向相邻间距比的告警阈值，留空使用算例预设，也可设置 `IPM_CHESS_RADIAL_RATIO_LIMIT` 或在算例中填写 `method.geometry.radialRatioLimit`。预检会打印实际最大径向网格比。

`threads=[]` 从 `SLURM_CPUS_PER_TASK` 或 `PBS_NCPUS` 读取配置。UI 未传出这些变量时，直接填写 `threads`。在同一个作业目录中续算可调整线程数，网格和算例由该目录的 `case.json` 固定。改网格时，在 UI 新建作业目录。

### 用一个网格级别设置点数

可以填写 `k` 并将四个手动计数留空。四向点数按
`[2^(k+1)+1, 2^k+1, 2^(k+1)+1, 90*2^(k-6)+1]` 计算。命令方式也可设置 `IPM_CHESS_K`，算例 JSON 可填写 `method.geometry.k`。

| k | nx | ny | nr | ntheta | 径向网格比告警阈值 |
|---|---:|---:|---:|---:|---:|
| 6 | 129 | 65 | 129 | 91 | 1.3 |
| 7 | 257 | 129 | 257 | 181 | 1.15 |
| 8 | 513 | 257 | 513 | 361 | 1.1 |

`caseName` 选择几何与演化预设，`k` 设置四向点数和表中的径向比告警阈值。显式填写 `radialRatioLimit` 可覆盖该阈值。直接选择粗网格预设时，径向告警阈值也为 1.3。

`k=7` 的 X1、X2 和 r 三向告警阈值均为 1.15，重网格时保留该级别和阈值。它只调整告警标准，原有核心节点分配保持不变。显式填写 `radialRatioLimit` 时只覆盖 r 方向。

### UI 的 MATLAB 命令方式

如果 UI 表单接收 MATLAB 命令，先在 CHESS 生成的目录里保存 `pwd`，再通过绝对路径调用提交脚本。因 MATLAB 的 `run` 会切换当前目录，这一步用来定位作业输出。

```matlab
setenv('IPM_CHESS_JOB_DIR',pwd)
% This example uses manual counts. Set k=[] in main_chess.m first.
setenv('IPM_CHESS_NX','225')
setenv('IPM_CHESS_NY','125')
setenv('IPM_CHESS_NR','465')
setenv('IPM_CHESS_NTHETA','161')
setenv('IPM_CHESS_THREADS','16')
run('/data/user/hd58131/ipm/ipm_cr_chess_configurable/server/main_chess.m')
```

这些环境变量优先于脚本顶部对应设置。直接提交脚本文件时，默认使用它在 CHESS 中的路径确定作业目录。

### 首次上机检查

先在命令前加 `setenv('IPM_PREFLIGHT_ONLY','1')`。它构造实际网格，打印四向点数，并检查 CHESS 分配的 CPU。随后在新作业用 `setenv('IPM_SMOKE_ONLY','1')` 按所选网格跑到 `0.002τ`，现场编译服务器平台的 MEX。核对 `native/summary.json` 中的 `status`、`threads`、`coreSize` 和 `polarSize`。

## 已分配计算节点上的终端入口

终端入口以新目录运行监督器。可在自己的算例 JSON 中填写 `method.geometry.nodeCounts`，并在启动命令指定 CPU。

```bash
export MATLAB_BIN=/absolute/path/to/matlab
export IPM_RUNS_DIR=/absolute/path/to/persistent/ipm_runs
bash server/launch.sh start --case /absolute/path/to/case.json --threads 16 --target 76 --hours 8
bash server/launch.sh resume /absolute/path/to/ipm_runs/ipm_cr_...
```

服务器包对应 FP64 B 样条求解器。实际 CHESS 节点上的 MATLAB 版本、MEX 工具链与线程配置以服务器短测为准。
