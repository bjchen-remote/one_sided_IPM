# ipm_cr 的公共数值支持

在 `main` 分支中，本目录仅保留 [SHARED_DEPENDENCIES.json](../test/SHARED_DEPENDENCIES.json) 列出的 14 个原路径公共函数，供 `test/+ipm_cr` 使用；函数内容保持原样。

完整 Cartesian 和纯极坐标求解器、各自的入口、测试、运行管理器与算法文档见 [legacy/solvers 分支](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers/time_evolution)。根历史原版也在该分支中。

主求解器、构建与运行说明见 [项目 README](../README.md) 和 [混合网格说明](../test/README.md)。本目录的公共支持子集不提供独立旧算法运行入口。保留的 `ipm_c.config` 用于原始初值与规范配置，其旧时间推进配置不用于当前 RK4 主链。
