# 历史 `test` 目录导航

维护求解器已迁至 [`src/+ipm_cr`](../src/+ipm_cr)，统一用法见根 [README](../README.md)。当前测试在 [`tests`](../tests)，构建及运行管理在 [`tools`](../tools)，独立演示在 [`examples`](../examples)。

- [架构与接口](../docs/ARCHITECTURE.md)
- [网格连接与验证边界](../docs/CONNECTION_METHOD.md)
- [迁移、清理与恢复记录](../docs/CLEANUP.md)
- [算法报告](ALGORITHM_REPORT.tex)：保留编辑器中的原文件位置，源码/API索引随当前实现更新。

本目录中的 `results/`、`prototypes/`、旧源码及旧工具是本地历史证据或隔离研究，不是当前默认 MATLAB 路径。冻结 runtime、检查点和原型保持原位；请按各自记录的源身份复现，不将其替换为新包后继承原资格结论。完整旧求解器保存在 [`legacy/solvers`](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers)，分支关系见 [REPOSITORY.md](../REPOSITORY.md)。
