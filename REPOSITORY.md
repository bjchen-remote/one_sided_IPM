# 主分支、旧分支与复现边界

`main` 以 `src/+ipm_cr` 为唯一维护求解器。完整旧算法保存在 [`legacy/solvers`](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers)，其历史分离点为 `2010a4e2e0714eb092ee0a6207e834b52844ced8`。旧分支保留当时全部 595 个跟踪文件，也包含迁移前的混合网格源码，供跨算法对照和历史复现；后续主线目录迁移不修改该分支。

## 主线目录

| 目录 | 职责 |
|---|---|
| `src/+ipm_cr` | 求解器公开入口、context 与职责包；不依赖旧 namespace |
| `data/initial.json`、`data/reference` | 统一初始文件定位、原 SHA、物理/规范设置和原始 τ0 数据 |
| `tools` | 原生构建、输入冻结、Python 管理器与公共工具 |
| `tests` | 组件验证、制造解、完整轨道研究及旧研究调用的薄适配 |
| `examples` | 独立椭圆演示和只读科研绘图 |
| `docs` | 架构、连接方法与清理迁移记录 |
| `test/ALGORITHM_REPORT.tex` | 原位保留的算法报告，当前源码索引与历史数据分开维护 |

所需原公共函数已按职责进入 `mesh/model/gauge/op/io`；旧 `test/SHARED_DEPENDENCIES.json` 只描述迁移前的 14 项共享支持，不是当前运行依赖清单。请只将 `src` 及所需 `tools/tests/examples` 加入当前 MATLAB 路径，避免将本地保留的旧包混入新实验。

## 旧算法归属

| 旧分支中的路径 | 算法或用途 |
|---|---|
| 根 `main.m`、`solver/` | 历史 Cartesian 原版及历史稳态分析 |
| `time_evolution/main.m`、完整 `+ipm` | Cartesian、实用规范、SDIRK2、重网格与旧验证 |
| `time_evolution/main_c.m`、完整 `+ipm_c` | 纯极坐标、B 样条 Poisson、SSPRK3 与径向自适应 |
| `time_evolution/server`、旧管理器和旧测试 | 上述求解器的运行、部署与验证 |
| 根 `plotting/` | 纯极坐标保存数据的读取和示例 |
| `test/tests/inspect_ipm_cr_explicit_initial.m` | 依赖完整旧纯极坐标初始化的初态对照 |
| `test/tests/compare_ipm_cr_explicit_tau200.m` | 历史 hybrid/polar τ200 轨道比较 |
| `test/plots/plot_ipm_cr_tau10_radial_diagnostic.m` | 绑定历史冻结 runtime 的诊断绘图 |

本地旧文件保持原位，在 `main` 中取消跟踪；它们不是自动删除目标。需要旧算法时，在独立目录克隆：

```sh
git clone --branch legacy/solvers --single-branch https://github.com/bjchen-remote/one_sided_IPM.git one_sided_IPM_legacy
```

共享工作区中的原型和旧目录可能存在并行研究依赖，不应通过切分支、批量清理忽略文件来切换运行环境。

## 构建、初始输入与身份

用法见 [README](README.md)，模块/状态契约见 [架构](docs/ARCHITECTURE.md)，迁移范围与历史归档见 [清理记录](docs/CLEANUP.md)。MEX 二进制不入库，保留 C++ 源码及构建入口。`config` 集中完整数值默认，`pconfig` 只定义谱系统，初始输入路径和原始物理设置集中在 `data/initial.json`，不再在多处硬编码历史 `runs` 路径。

原 τ0 authority 内容保存在 `data/reference/gaussian_t0.mat`，约 4.35 MiB，是唯一随 Git 保存的运行数据；其 SHA-256 保持为：

```text
fa80e7d2b4fedb6b28186dcaee844975569fac41da06d26b35a2519d7521db22
```

历史输入路径 `runs/local_validation_20261013/gaussian26_fresh_warm_01/data_frames/frame_step000000000000_tau0000000000000000.mat` 在本地与冻结证据中保持原位。新包通过清单定位输入，不能将旧路径身份误称为新路径身份。

新运行按相对路径和内容 SHA 记录源码/输入，冻结副本只复制声明的当前实现与输入。移动完整冻结目录可保持内容身份；变更源码、方法、初值或缺少时钟低位不属于合法恢复。context 的私有计划在冷恢复时重建，不以反序列化旧 handle 绕过身份核验。

## 本地数据与历史证据

AI 脚手架、checkpoint、轨道、日志、锁、编译产物和大型输出不入库。`experiments/`、`test/prototypes/` 的并行研发保留本地，不是本次删除目标；历史结果、旧输入和 frozen runtime 均不重写。

文档指向 `test/results/`、`output/`、历史 `runs/` 或 `releases/` 的链接表示本地历史证据，干净克隆中不一定存在。`test_ipm_cr_refactor` 的旧 baseline MAT 及历史轨道需另行提供。恢复旧实验使用其冻结源码和原记录身份；当前 API/目录迁移不替代新长程精度试验，也不能给旧轨道改标新来源。

本次实施以保持方程、RK4 算术、CFL、连接和 Poisson 数值语义为约束；实际等价性由对应回归输出说明。尚未完成的精度控制与 τ200/50分钟目标保持未完成，不因架构改善而视为完成。
