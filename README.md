# ipm_cr：单侧 IPM 混合网格求解器

主分支 `main` 维护 `ipm_cr`：均匀 Cartesian 张量核心、不含原点的极坐标环带、独立谱 Poisson 求解，以及完整物理 CFL 约束的显式 RK4。解析背景与残差分离用于控制远场和舍入误差。

[使用与验证](test/README.md) · [算法报告](test/ALGORITHM_REPORT.tex) · [架构](test/ARCHITECTURE.md) · [网格连接](test/CONNECTION_METHOD.md) · [仓库与分支说明](REPOSITORY.md)

## 分支与目录

- `main`：当前混合网格求解器、维护测试、构建工具、运行管理器和算法文档。
- [`legacy/solvers`](https://github.com/bjchen-remote/one_sided_IPM/tree/legacy/solvers)：完整旧 Cartesian 求解器、纯极坐标求解器、历史原版、旧绘图与对照入口。该分支保留迁移前的完整快照，含当时的混合网格代码，便于复现。
- `test/+ipm_cr`：当前求解器；保留目录层级以兼容既有相对路径和原始初值身份。
- `time_evolution`：主分支中只保留混合求解器实际调用的 14 个公共数值函数；完整旧求解器在旧分支。

## 构建与运行

需要 MATLAB、Python 3；原生内核需要兼容的 C++ 编译器。先将 MATLAB 当前目录设为项目根，再执行：

```matlab
restoredefaultpath;
addpath('test','test/tests','time_evolution');
build_ipm_cr_native;             % 在本机构建，不上传 MEX 二进制
[~,~,diagnostic]=main_cr;        % 独立椭圆演示，不推进含时系统
assert(diagnostic.accepted);
```

含时运行、固定远区径向重分配、停止与续算示例见 [test/README.md](test/README.md)。运行管理器按 [公共依赖清单](test/SHARED_DEPENDENCIES.json) 冻结实际输入，不复制本地保留的整套旧算法。

原始 τ0 authority 是唯一随 Git 保存的运行数据。检查点、轨道、日志、编译产物、AI 脚手架与本地研发原型均不入库。源码管理与分支整理不新增 τ200 精度或速度保证；当前证据及未完成控制见算法报告。
