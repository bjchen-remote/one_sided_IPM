# 极坐标数据的 MATLAB 绘图示例

把本文件夹加入 MATLAB 路径，即可得到一个普通 `data` 结构体：

```matlab
addpath('/Users/cbjin/Documents/one_sided_IPM/plotting');
load('/Users/cbjin/Documents/one_sided_IPM/plotting/data_polar.mat', 'data');
```

这里打包的 `data_polar.mat` 来自修正版算例07在 `tau=200` 的保存态：包含完整原生极坐标场，以及用现有三次有限元插值到 `[-2,2] × [0,4]` 的 `201 × 201` Cartesian 局部场，并保留完整的 `t/tau/c_l/c_omega/c_r` 等标量轨道。需要其他时间、区域或分辨率时，继续使用下方的 `load_polar_data(source, tau, box, resolution)`；不传参数时读取修正版算例的最新保存切片。读取和绘图只使用已保存的数据。

## 直接画图

```matlab
% 局部密度图：选取最接近 tau=50 的实际保存切片。
data = load_polar_data([], 50, [-2 2 0 4], 301);
figure;
imagesc(data.x, data.y, data.rho);
axis xy equal tight;
xlabel('X'); ylabel('Y'); colorbar;

% c_r 关于物理时间 t 的整条轨道。
figure;
plot(data.t, data.c_r, 'LineWidth', 1.4);
xlabel('t'); ylabel('c_r');

% 同一条轨道，改用重标度时间 tau 作横轴。
figure;
plot(data.tau, data.c_r, 'LineWidth', 1.4);
xlabel('\tau'); ylabel('c_r');
```

也可以直接运行本目录的六个脚本：

| 文件 | 内容 |
| --- | --- |
| `example_local_rho.m` | 使用 `data.x/data.y/data.rho` 画局部密度图 |
| `example_local_rho_x.m` | 使用 `data.rho_x` 画局部导数场 |
| `example_c_r_vs_t.m` | 并排画 `c_r(t)` 与 `c_r(tau)` |
| `example_velocity.m` | 画相对 Darcy 速度箭头 |
| `example_polar_radius.m` | 固定半径，画随角度变化的剖面 |
| `example_polar_ray.m` | 固定角度，画随半径变化的剖面 |

脚本开头的 `tau`、`box`、`resolution` 可以直接修改。示例只开图；需要保存时自行使用 `exportgraphics(gcf,'my_plot.png')`。

## data 中的主要字段

| 字段 | 含义 |
| --- | --- |
| `x`, `y` | Cartesian 横、纵坐标轴向量；分别为 `1 × Nx`、`Ny × 1` |
| `X`, `Y` | 对应的二维坐标网格 |
| `rho` | 当前切片的完整重标度密度 `R(X,Y,tau)`，尺寸 `Ny × Nx` |
| `rho_x`, `rho_y` | 当前切片的 `R_X`、`R_Y`；大写别名为 `rho_X`、`rho_Y` |
| `W`, `rho_centered` | 当前切片的 `R-R(0,tau)`，便于观察后期空间形状 |
| `rho_origin` | 当前切片的标量 `R(0,tau)`，满足 `rho = W + rho_origin` |
| `u1`, `u2`, `psi` | 当前切片的重标度 Darcy 速度和流函数 |
| `u1_relative`, `u2_relative` | 相对保存态参考点的 Darcy 速度；`u1_relative = u1-reference_velocity` |
| `vx`, `vy` | 包含膨胀、平移项的重标度输运速度 |
| `Q`, `rho_tau` | 当前切片的 `D_X(R_tau)` 和完整 `R_tau` |
| `mask` | 插值域内的有效点；域外场值保留 `NaN` |
| `t`, `tau` | 所有已保存接受步的物理时间、重标度时间向量 |
| `c_l`, `c_omega`, `c_r` | 与 `t/tau` 对齐的完整系数时序 |
| `c` | 三列分别为 `c_l/c_omega/c_r` |
| `Q_inf`, `Qinf` | 每个接受步的原生网格 `Q` 无穷范数 |
| `selected` | 当前空间切片的实际 `tau/t`、系数、帧号、网格标识等 |
| `fieldTau`, `frames` | 全轨道已保存空间切片的时间列表与原文件目录 |
| `polar_grid` | 当前切片实际使用的原极坐标网格，含 `r/theta/meshId` |
| `polar` | 当前切片的原生极坐标场，数组尺寸为 `Nr × Ntheta`，未经 Cartesian 插值 |
| `grid`, `metadata` | 当前 Cartesian 网格和有限元插值信息 |
| `meta` | 数据源路径、请求与实际选取时间、变量定义 |

场数组对应 `data.selected.tau`；切换空间切片会保留完整的 `t/tau/c_r` 等标量轨道。默认按最近保存时间选帧，可以查看实际取值：

```matlab
data.selected.tau
data.selected.t
data.selected.c_r
data.fieldTau
```

坐标 `x/y` 表示重标度的 `X/Y`；`rho/rho_x/rho_y` 也保留重标度量的定义。物理时间 `t` 在演化后期逼近有限极限，`c_r(t)` 曲线会挤在右端；查看长时间系数收敛时可使用 `data.tau`。

## 选择算例、区域与时间

```matlab
% 修正版07；也可以传 native 目录或 native/data.mat 的路径。
run = '/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_originfix_20261021_07';
data = load_polar_data(run, 100, [-1 1 0 2], [301 201]);

% 原算例06的数据保留，可独立读取。
oldRun = '/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_fresh_v4_20261021_06';
oldData = load_polar_data(oldRun, 100, [-1 1 0 2], [301 201]);

% 把当前结构体另存为自己的 MAT 文件。
save('my_data_snapshot.mat', 'data', '-v7.3');
```

`resolution` 为 `[Nx Ny]`，传一个整数表示两个方向使用相同点数。`box` 为 `[xmin xmax ymin ymax]`。原极坐标场片及每次换网格的信息继续保存在各算例自己的 `native/data_frames` 与 `native/grids` 中。

## 原生极坐标剖面

新增的数据包 `data_polar.mat` 包含原有 Cartesian 数据和 `data.polar`；之前的 `data.mat` 原样保留。调用新版 `load_polar_data` 也会同时得到两种坐标下的数据。`data.polar` 保存选定时刻的完整原网格，与 Cartesian 绘图窗口 `box` 无关。

```matlab
data = load_polar_data([], 50);
p = data.polar;

% 固定半径 r≈1，沿半圆画 R(theta)。选最近的实际径向节点。
[~, i] = min(abs(p.r - 1));
figure; plot(p.theta, p.rho(i,:));
xlabel('\theta'); ylabel('R');
p.r(i)                       % 实际选取的半径

% 固定角度 theta≈pi/4，沿该射线画 R(r)。
[~, j] = min(abs(p.theta - pi/4));
use = p.r >= 0 & p.r <= 4;
figure; plot(p.r(use), p.rho(use,j));
xlabel('r'); ylabel('R');
p.theta(j)                   % 实际选取的角度，单位为弧度
```

`p.r` 是 `Nr × 1`，`p.theta` 是 `1 × Ntheta`，角域为 `[0,pi]`。`p.rho/p.rho_x/p.rho_y/p.W/p.Q/p.psi/p.u1/p.u2/p.vx/p.vy` 均为 `Nr × Ntheta`；`p.rho` 是完整密度，`p.W` 是减去原点值后的密度。固定半径取一行，固定角度取一列；示例按实际网格节点取值，不对剖面额外平滑或插值。若更换为其他场，需同时修改纵轴标签。

另有 `p.rho_r`、`p.rho_theta`（分别为 `dR/dr`、`dR/dtheta`），以及 Darcy 速度分量 `p.ur/p.utheta` 和输运速度分量 `p.vr/p.vtheta`；`utheta/vtheta` 表示切向线速度。`p.selected` 和 `p.grid` 保留当前剖面的实际时间与网格信息。全域半径达到 `1e10`，观察核心区时建议像示例一样指定径向范围。
