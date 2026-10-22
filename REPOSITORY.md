# 源码仓库与复现边界

本仓库管理单侧 IPM 的维护源码、测试、原生内核源码和算法文档。Git 边界由根目录 `.gitignore` 明确约束；本地未跟踪文件保持原位，不因初始化仓库而删除。

## 代码入口

| 路径 | 用途 |
|---|---|
| `time_evolution/+ipm`、`time_evolution/main.m` | Cartesian 时间演化维护版 |
| `time_evolution/+ipm_c`、`time_evolution/main_c.m` | 极坐标维护版 |
| `test/+ipm_cr`、`test/tests/run_ipm_cr_explicit.m` | 张量核心与极坐标环带的混合网格显式求解器 |
| `test/ALGORITHM_REPORT.tex` | 混合网格算法报告，中文独立 LaTeX 源文件 |
| `test/ARCHITECTURE.md`、`test/CONNECTION_METHOD.md` | 模块契约、网格连接与误差预算 |
| `main.m`、`solver/` | 历史原版与历史分析代码，不能代替当前维护版的验证 |
| `plotting/` | 已保存轨道的读取和绘图示例 |

运行需求为 MATLAB 和 Python 3。Python 管理器的解释器与 MATLAB 路径需按本机设置；部分资源监控及直接 CLT 构建接口限 macOS。原生 MEX 二进制不入库，C++ 源码及构建入口入库。混合路径需要兼容的 MEX 时，可在 MATLAB 中执行：

```matlab
addpath('test','test/tests','time_evolution');
build_ipm_cr_native; % 在当前平台构建，并将日志写入独立结果目录
```

求解器使用方法和测试入口分别见 [维护版说明](time_evolution/README.md)、[极坐标说明](time_evolution/IPM_C.md)与[混合网格说明](test/README.md)。请勿同时将历史 `solver/+ipm` 和维护 `time_evolution/+ipm` 放在 MATLAB 路径前端，以免同名包相互遮蔽。

## 必须随源码保存的初始输入

默认极坐标与混合初始化读取：

```text
runs/local_validation_20261013/gaussian26_fresh_warm_01/data_frames/frame_step000000000000_tau0000000000000000.mat
```

该约 4.35 MiB 的文件是原始 τ0 数值 authority，作为唯一运行数据例外随 Git 保存。保持原路径可避免改变数值源码与历史输入身份。其他检查点、轨道和 `data.mat` 不入库。

## 不纳入 Git 的内容

`skills/`、`AGENTS.md`、`agent.m`、`.codex/` 等 AI 脚手架不属于求解器。当前维护源没有对它们的运行依赖。冻结源码、运行结果、日志、锁、MEX 编译物、本机缓存和大型实验输出均忽略；`experiments/`、`test/prototypes/` 的活跃研究代码本轮保持独立，本地保留。

`test/run_explicit_experiment.py`、`test/run_validation.py` 与维护目录的 `run_*.py` 是数值实验管理器，提供资源限制、源码冻结、日志和续算能力，因此保留。

## 验证证据的范围

文档中指向 `results/`、`output/`、历史 `runs/` 或 `releases/` 的链接表示本地历史证据，它们不随本源码仓库上传；克隆后这些链接可能缺少目标文件。读取旧轨道、重画历史对比图或运行需要外部检查点的审核时，必须另行提供对应输入。`test_ipm_cr_refactor` 还要求清理前的 baseline MAT，不能作为克隆后无需额外输入的测试。

算法报告区分已有 τ10 采样精度证据、固定端点诊断和仍未完成的控制。源码管理不新增 τ200 精度认证或运行时间保证。数值实验会创建独立输出目录，输出文件不会自动成为版本管理对象。
