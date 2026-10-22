# 活动代码清理记录

## 2026-10-22：Git 主线与旧算法分离

`main` 以 `ipm_cr` 为唯一维护求解器；完整旧算法保留于 [legacy/solvers](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers)，对应迁移前提交 `2010a4e2e0714eb092ee0a6207e834b52844ced8`。根历史原版、旧 Cartesian/纯极坐标入口、各自测试与绘图，以及三个依赖旧运行环境的对照入口不再由主分支跟踪。本地文件保持原位，不作物理删除。

主线保留 [SHARED_DEPENDENCIES.json](SHARED_DEPENDENCIES.json) 声明的 14 个公共函数，原路径和函数内容均不变；冻结管理器与含时入口按清单记录实际共享输入，避免将本地保留的完整旧包混入新实验。当前数值方法和兼容算子不因分支整理改变。下面较早的“保留”记录描述的是各次整理当时的范围；当前 Git 归属以本节及根 [REPOSITORY.md](../REPOSITORY.md) 为准。

本次从暂存区导出 82 个主线文件，在无完整旧求解器的干净 MATLAB 路径上重建原生内核，完成局部数学检查、独立椭圆演示、两种径向网格短 RK4 轨道的逐位回归，以及双槽检查点续算到 τ=.006。均通过；数值阶段合计约 9.23 秒，仅为维护验证耗时，不是长轨道性能测量。两个 Python 管理器的 prepare-only 打包检查也通过。详细本地证据见 `results/branch_split_20261022/clean_validation/validation.json`。

本文件记录从维护面移出的依赖闭包、归档核验和恢复方式。旧结果、检查点、原始 authority、冻结 runtime 以及 `time_evolution/+ipm/+ipm_c` 保持原身份。代码优化与等价性验证另外记录，不把代码清理计为新的数值精度证据。

## 2026-10-22：当前 RK4 维护面整理

本轮从 `test/tests` 移出 **7 个历史专属入口，共 100,270 字节**。删除前建立维护源快照，要求候选 SHA 与快照一致；随后按原路径归档，逐成员核对 SHA 后才删除。

| 归档文件 | 退出当前维护面的依据 |
|---|---|
| `test_ipm_cr_connection_orbit.m` | 实验时间格式固定为 SSPRK3，不是当前 RK4 链；其连接穿越设计保留在冻结证据，当前 RK4 输运制造解由 `test_ipm_cr_fourth_order` 维护 |
| `inspect_ipm_cr_connection_orbit_refinement.m` | 只读取上述历史三层与 exact-ghost 对照输出，无当前维护代码调用 |
| `diagnose_ipm_cr_far_poisson_chain.m` | 严格定位旧 τ10 冻结 runtime 的固定尺度制造解 A/B/C 分解，不是通用运行入口 |
| `diagnose_ipm_cr_far_transport.m` | 使用旧 τ10 接受态尺度和其冻结内核的静态一维径向诊断 |
| `diagnose_ipm_cr_moving_gaussian.m` | 使用旧 τ10 冻结内核的独立一维移动 Gaussian 诊断，不推进 IPM 系统 |
| `diagnose_ipm_cr_tau10_operators.m` | 读取旧 τ10 原始接受态，比较宽模板的固定状态诊断 |
| `diagnose_ipm_cr_tau10_radial_refinement.m` | 固定旧 τ10 中层方法、只允许 `nr/firstSpacing` 变化的单次原初值控制实验 |

这些入口没有被当前含时入口、维护测试、绘图和运行管理器调用。删除闭包前的扫描覆盖 `test` 与 `time_evolution` 的非结果源码；历史 `results/**` 对原入口的引用仍是历史复现依赖，不能改写。项目内 skill 导航的旧穿越入口指向同步更新为当前 RK4 与历史归档的区分。

恢复归档：[releases/ipm_cr_diagnostics_20261022.tar.gz](../releases/ipm_cr_diagnostics_20261022.tar.gz)。

```text
SHA-256 a67daa4fb4f17956568e247585fe8546c5ec38c4eed25dcc9a534297aa06212d
```

[逐文件清单](results/cleanup_20261022/archived_components.json) 记录原路径、字节数与 SHA；[归档核验](results/cleanup_20261022/archive_verification.json) 记录删除前的成员核验。需要复现时先解到独立空目录，校验 SHA 和相应实验的 frozen runtime，再按原路径调用，不能让当前方法替代原依赖。

### 保留的隔离研发与独立工具

`test/prototypes/` 是并行稳态加速研发区，仍有活跃任务及跨目录依赖，本轮**未移动或改写**。`steady_residual_optimization` 依赖 `streamline_fixedpoint` 的投影和候选构造；`residual_stall_diagnostic` 又读取其已验证父运行。目录独立于当前 RK4 主运行链，并不代表无用或失败可删。审计过程中捕捉到新文件变化后，原候选归档在任何删除前中止，最终只归档上述 7 个不变文件。

保留 `diagnose_ipm_cr_poisson_backward` 的通用原方程后向误差审计、`inspect_ipm_cr_explicit_initial/endpoint`、`compare_ipm_cr_explicit_tau200` 及历史只读绘图工具，它们有独立用途。保留 `+carrier/localMap.m/tensorMap.m`，因为当前 ghost 和源映射仍调用它们。保留显式选择的 `legacy_weno5/ssprk3` 兼容路径；本轮清理不更改算法选项或删除原数学回退。

### 测后撤回的缓存候选

本轮尝试为物理差分保存固定稀疏行/列计划，逐字段等价回归通过。但中层差分微基准的中位数由旧实现 `3.96452 ms` 增至缓存候选 `4.20885 ms`，约慢 6%。因此恢复 `geometry.m/differentiate.m` 原实现，移出候选 `internal.derivativePlan`；候选源和测试记录另存于 [本轮维护证据](results/cleanup_20261022)，未作为生产优化采用。当前主链不持有 `derivativePlans`。

保留的代码整理是公共 `sampleState` 与空间研究重复采样的消除，以及维护回归入口 `test_ipm_cr_refactor`。等价回归覆盖其记录的完整 RHS、采样及短程 RK 轨道，不代替新的长程精度实验；本轮没有完整 RHS 或整机提速承诺。

### 当前目录分工

```text
test/
  +ipm_cr/
    Solver.m, config.m          独立谱 Poisson 管理
    +core/, +polar/             核心/环带谱算子
    +internal/                  内部谱与固定几何工具
    +explicit/                  RK4/η 主链与公共只读采样
    +carrier/                   仍在使用的局部张量插值
  tests/                        维护回归、制造解、完整轨道及通用诊断
  plots/                        保存数据的只读绘图
  prototypes/                   并行研发隔离区；不属于自动清理目标
  main_cr.m, run_validation.py   独立椭圆入口
  run_explicit_experiment.py     冻结含时运行管理器
  results/                      实验数据、检查点、冻结运行和清理证据
  README.md                     使用、证据索引和维护地图
  ARCHITECTURE.md                模块、数据、缓存与事务契约
  CONNECTION_METHOD.md          网格连接和误差预算
  ALGORITHM_REPORT.tex          连续方程、离散推导与算法报告
```

## 2026-10-22：早期资格入口归档

此前三个 SSPRK3 专属资格入口已由当前 `test_ipm_cr_fourth_order`、`test_ipm_cr_ivp_spatial` 维护链取代，移出活动目录：`qualify_ipm_cr_explicit.m`、`refine_ipm_cr_explicit_window.m`、`test_ipm_cr_explicit_transport.m`。历史结果和冻结运行保留。

恢复归档：[`releases/ipm_cr_legacy_qualification_20261022.tar.gz`](../releases/ipm_cr_legacy_qualification_20261022.tar.gz)。三个原文件与归档逐项 SHA 核验后才删除，归档含路径清单。

```text
SHA-256 ead5d437238669f7b74ca32a6f6af5d19ad647854260a62e01db62506ba8d6d0
```

## 2026-10-21：退役分支归档

当日移出 20 个退出当前求解链的源码或说明文件，共 128,696 字节。删除前已按原路径归档并逐成员核对 SHA；原维护源码、结果、检查点和冻结副本保持不变。

| 范围 | 文件数 | 依据 |
|---|---:|---|
| `+ipm_cr/+evolve/{initialize,nativeFlow,transport,interpolate}.m` | 4 | 旧原生谱输运分支，当前 `+explicit` 与连接验证不再调用 |
| `+ipm_cr/+gauge/measure.m` | 1 | 仅供上述旧分支使用 |
| `+ipm_cr/+carrier/{poisson,evaluate}.m` | 2 | 旧极坐标输运替换 Poisson 后端包装，不属于双网格演化 |
| `abandoned_time_integrators_20261021/` | 6 | 已撤出的隐式、半拉格朗日草稿、说明和试验 |
| `run_tau200_attempt.py` | 1 | 旧准入实验组织入口 |
| `tests/{test_ipm_cr_carrier,preflight_ipm_cr_tau200,benchmark_ipm_cr_carrier,pilot_ipm_cr_tau200,refine_ipm_cr_tau200_endpoint,optimize_ipm_cr_endpoint}.m` | 6 | 仅服务上述旧包装与准入实验 |

当日静态闭包检查确认，剩余 39 个活动 `.m/.py/.sh` 文件对被移出模块和旧入口无引用；`results/**` 的 379 个已有文件前后路径、大小及修改时间相同。当日未运行 MATLAB，未形成新的精度结论。

恢复归档：[`releases/ipm_cr_legacy_20261021.tar.gz`](../releases/ipm_cr_legacy_20261021.tar.gz)，39,632 字节，含 20 个原路径文件及 `ipm_cr_legacy_manifest.json`。

```text
SHA-256 384c69629b50ba0b7b3603d26f1fc6e493c5668d0c4e82d84c42bba43eb96c8b
```

历史结果里的 `execute.m` 可能引用当时的活动入口。复现应恢复原依赖并核对记录的源码身份，不能借用当前求解器后宣称原实验得到复现。
