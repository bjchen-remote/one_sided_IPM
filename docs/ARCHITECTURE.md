# 单侧 IPM 的统一 B 样条求解器

`main.m` 读入算例，`run.m` 管理接受态和保存，`context.m` 为实际网格准备共享样条与 Poisson 因子。演化、网格更新与独立验证均使用同一套空间离散。

| 论文章节 / 职责 | 文件 | 数量 |
|---|---|---:|
| 入口与演化 | `main.m`、`config`、`init`、`context`、`run`、`rhs`、`step`、`remesh`、`sample`、`refresh` | 10 |
| 初值、背景和初始规范 | `+model/initial`、`background`、`reference` | 3 |
| 网格与转移 | `+mesh/build`、`axes`、`staircase`、`exchange`、`interpolate`、`region`、`transfer` | 7 |
| 共享离散与椭圆问题 | `+space/basis`、`poisson`、`boundary`、`velocity`、`transport` | 5 |
| 墙规范 | `+gauge/observe`、`peak`、`rates` | 3 |
| 时间细节 | `+time/clock`、`knot` | 2 |
| 存储 | `+io/checkpoint`、`identity`、`history`、`export` | 4 |
| 被动诊断 | `diagnostics` | 1 |
| 本机采样 | `+mesh/sampleMex.cpp`、`tools/build_ipm_cr_native.m` | 2 |
| 外部监督 | `tools/run_explicit_experiment.py`、`tools/common.py` | 2 |
| 合计 | 计入正式运行、采样构建和监督器 | **39** |

`src/+ipm_cr` 下有 34 个 MATLAB 文件和一个 C++ 文件。表中的 MATLAB 文件均相对于该目录，入口和 `tools/` 除外。测试、离线初值制作工具、论文和配置数据独立维护。

## 完整 RHS

1. `mesh.exchange` 由另一片的活动 donor 同步三层 ghost。内侧为阶梯裁边的非均匀笛卡尔张量网格，外侧为不含原点的极坐标张量网格。活动重叠点各自推进。ghost 每轴用七点重构。
2. `space.basis` 在原生 `x/y/log(r)/theta` 上构造共用的开结点五次 B 样条。由原生点值求系数，再求 $q=D_1\eta+H_{X_1}$。采用完整密度表示时直接求 $q=D_1W$。
3. `space.poisson` 用同步后的完整 $q$ 组装载荷。阶梯内区缓存稀疏质量算子与对称 Nitsche 迹载荷，外区载荷为 $M_{t,r^2}C M_\theta^{\mathsf T}$。六点 Gauss 点用于单元积分。每次换网格重建载荷、因子及 Schur 系统。
4. `space.boundary` 从角向质量及刚度得到广义特征模态，为圆盘外零源模型施加齐次 DtN。墙上 $\psi=0$，上侧壁源可以非零。`space.velocity` 从同一个流函数恢复速度，并在核心样条支撑之外读取环带解。原生流函数和速度的评价位置与输运节点完全相同。内区人工迹使用 Nitsche 面求积点，外区人工迹使用原生角点。
5. `gauge.observe` 和 `gauge.rates` 在原墙网格计算 I/P/C 和尺度速率。`space.transport` 对每条轴求 $M^{-1}(C\pm P)$ 的两支弱导数，按局部速度符号取支。$P$ 惩罚五阶导数的结点跳跃，缩放系数为 `1/65536`，并乘较小相邻结点跨度的十次方。
6. `rhs` 组装 $F=W_\tau$ 和 $Q=q_\tau$。有背景时推进 $\eta=W-H$，并用共用强导数和显式背景导数组装完整残差。Poisson 在生产阶段只读取原生节点上的完整源。

## 时间与重网格

`step` 固定使用经典 RK4，复用起点 RHS，并计算可供下一步复用的完整末端 RHS。初始步长为 `min(maxDt, .95*cfl/rate, tauEnd-tau)`，默认 `cfl=.4`。原墙样条的连续峰跨相邻结点时，`time.knot` 用最多四次割线更新确定分段。物理时钟的高字、低字和两个尾字随场使用同一 RK 表。

`mesh.region` 测量 25% 和 80% 主峰连墙超水平集。每层比较横向最大绝对范围，以及穿主峰的纵剖面从墙到首个等值线交点的高度。任一层任一方向相对上次实际换网格的锚点变化达 25% 即触发。首次观测建立锚点，后续根据形态变化率调整观测间隔。粗网格默认从 $\tau=50$ 起固定网格。

触发后以 25% 层完整轮廓的横、纵量相对初始参考缩放坐标轴，重建内缘匹配的非均匀环带。角轴和最外半径保持固定。内外网格交界处分配近似相当的物理步长，活动重叠逐角记录。六点张量转移优先使用同一片旧活动 donor，再重建 context 并求一次完整 RHS。原墙参考、尺度、四字时钟和步数随主态传递。

## 持久化与复现

`checkpoint` 记录完整场、实际网格、尺度、来源、方法及换网格谱系。写入后读回确认并校验 SHA，再以原子指针发布。冷恢复由实际网格重建 context，源身份随实现变化。旧轨道用各自的冻结运行时和 Git 版本读取。变更实现时可以在新目录显式导入真实 `file/sha256/tau` 的接受态，生成独立的来源记录。

[使用](USAGE.md) · [维护](MAINTENANCE.md) · [验证](VALIDATION.md) · [主干契约](MAINLINE_REQUIREMENTS.md)
