# 主分支、旧分支与复现边界

`main` 以 `ipm_cr` 为唯一维护求解器。旧算法保存在 [`legacy/solvers`](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers)，迁移点为 `2010a4e2e0714eb092ee0a6207e834b52844ced8`。该旧分支保留当时全部 595 个跟踪文件，也含迁移前的混合网格源码，方便复现跨算法对照。

## 旧算法归属

| 旧分支中的路径 | 算法或用途 |
|---|---|
| 根 `main.m`、`solver/` | 历史 Cartesian 原版及历史稳态分析 |
| `time_evolution/main.m`、完整 `+ipm` | Cartesian、实用规范、SDIRK2、重网格与旧验证 |
| `time_evolution/main_c.m`、完整 `+ipm_c` | 纯极坐标、B 样条 Poisson、SSPRK3 与径向自适应 |
| `time_evolution/server`、旧运行管理器和旧测试 | 上述旧求解器的运行、部署与验证 |
| 根 `plotting/` | 纯极坐标保存数据的读取和示例 |
| `test/tests/inspect_ipm_cr_explicit_initial.m` | 包含完整旧纯极坐标初始化的初态对照 |
| `test/tests/compare_ipm_cr_explicit_tau200.m` | 历史 hybrid/polar τ200 轨道比较 |
| `test/plots/plot_ipm_cr_tau10_radial_diagnostic.m` | 绑定历史冻结 runtime 的诊断绘图 |

本地旧文件保持原位，但在 `main` 中不跟踪，避免影响并行研究。直接在这个共享工作区切换旧分支可能遇到本地忽略文件冲突；要运行旧算法，优先克隆旧分支到独立目录：

```sh
git clone --branch legacy/solvers --single-branch https://github.com/bjchen-remote/one_sided_IPM.git one_sided_IPM_legacy
```

## 主分支的公共依赖

`test/SHARED_DEPENDENCIES.json` 明确列出 14 个仍需调用的公共函数。保留其原 namespace 和相对路径，不改函数内容：有限差分矩阵、原初值、原 I/P/C 泛函、径向网格、原始 authority、规范测量与求根种子，以及 SHA-256。其中 MATLAB WENO5 只供仍支持的显式兼容选项。

`time_evolution/+ipm`、`+ipm_c` 在主分支中是公共支持子集，不是可独立运行的旧求解器。请勿把旧分支的完整包加入当前 MATLAB 路径。Python 冻结器只复制清单内文件，防止本地忽略的旧源码混入新实验。

## 构建与初始输入

主入口和构建见 [README](README.md)，含时用法、测试地图和精度边界见 [test/README](test/README.md)。MEX 二进制不入库，C++ 源码和构建入口保留。Python 管理器的 MATLAB 路径需按本机设置；部分资源监控与直接 CLT 构建仅适用 macOS。

默认初始化读取下列原始 τ0 authority，约 4.35 MiB，是唯一随 Git 保存的运行数据：

```text
runs/local_validation_20261013/gaussian26_fresh_warm_01/data_frames/frame_step000000000000_tau0000000000000000.mat
```

保持原路径可避免改变数值定义及输入身份。其 SHA-256 为 `fa80e7d2b4fedb6b28186dcaee844975569fac41da06d26b35a2519d7521db22`。

## 本地数据与历史证据

AI 脚手架、检查点、轨道、日志、锁、编译产物和大型输出不入库；`experiments/`、`test/prototypes/` 的并行研发保留本地，不作为本次分支整理的删除目标。

文档中指向 `results/`、`output/`、历史 `runs/` 或 `releases/` 的链接表示本地历史证据，克隆源码后不一定存在。旧轨道绘图、历史审核及 `test_ipm_cr_refactor` 所需 baseline MAT 必须另行提供。历史报告保持原路径和原实验结论；它们不表示主分支仍维护整套旧算法。恢复旧实验时应使用其冻结源码或相应旧分支。

本次整理不改变数值方程、RK4、CFL、连接、Poisson 或兼容算子的数学实现，也不新增全局精度或长时间性能认证。
