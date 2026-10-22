# 极坐标 IPM：算法细节与公式推导

本文展开 `time_evolution/+ipm_c` 的具体数值方法。先读[整体架构](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ARCHITECTURE_REPORT.md)，再按需要查这里；运行结果及其源码版本归属见[验证记录](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_VALIDATION_NOTES.md)，使用命令见[接口说明](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C.md)。

源码核对日期：2026年10月21日。当前为 v4 中心状态表示，径向 WENO 已采用负半径延拓；Poisson 仍采用原四阶 B 样条空间，高模态原点斜率约束的独立候选尚未接入。此前冻结06轨道使用旧原点闭合，其数据不等于当前源码的全程验证。

|查阅内容|入口|
|---|---|
|物理系统、尺度变换和极坐标方程|[1. 连续方程](#equations)|
|中心状态、初值、节点和唯一原点|[2. 状态与网格](#state-grid)|
|按步骤理解源、边界、角模、样条系数与流场|[3. Poisson 求解流程](#poisson)|
|墙观测与两变量规范根|[4. 规范](#gauge)|
|上风导数、有效 epsilon 和原点延拓|[5. WENO 输运](#transport)|
|步长、SSPRK3、物理钟和 RHS 缓存|[6. 时间推进](#time)|
|触发、转移、验收、提交与回滚|[7. 换网格](#remesh)|
|代数误差与独立采样强残差|[8. 空间审计](#audits)|
|字段、物理量恢复和绘图插值|[9. 数据与绘图](#data)|
|检查点恢复和工作量|[10. 恢复与复杂度](#implementation)|

## 记号约定

全篇使用以下约定，代码字段名只在接口部分列出：

|对象|唯一记号|
|---|---|
|物理坐标与场|\((x,y,t;\rho,\psi,\mathbf u)\)|
|重标度坐标与连续场|\((X,Y,\tau;R,\Psi,\mathbf U)\)|
|连续源、密度右端|\(q=\partial_XR\)、\(F(R)\)|
|完整密度与实际存储|\(\mathbf R=R_\circ\mathbf1+\mathbf W\)，\(\mathbf W_{0j}=0\)|
|节点源、完整右端、诊断|\(\mathbf q,\mathbf F_h,\mathbf Q\)|
|角向模式、径向模态函数|\(m\)、\(\phi_m(r)\)|
|径向 B 样条、展开系数|\(N_a(r)\)、\(z_{am}\)|
|内部模态线性系统|\(\mathbf A_m\mathbf z_m=\mathbf b_m\)|
|逻辑径向坐标、首间距、区间宽|\(\xi,h_0,\Delta r_i\)|
|时间编号和时间步|\(n,\Delta\tau_n\)|

下标 \(h\) 只标识空间离散或重构，不表示时间步。\(\mathbf R,\mathbf q,\mathbf F_h,\mathbf Q\) 是节点数组；\(\mathbf u,\mathbf U,\mathbf V\) 是连续速度向量；\(\Psi_h\) 是重构函数。在当前状态处可将 \(\mathbf F_h(\mathbf R)\) 简记为 \(\mathbf F_h\)。模态编号 \(m\)、极坐标节点编号 \(i,j\)、基函数编号 \(a,b\)、时间编号 \(n\) 分工固定。

<a id="equations"></a>

## 1. 连续方程与坐标变换

### 物理方程与边界

物理半平面为 `y≥0`。采用如下符号约定：

\[
\rho_t+u_1\rho_x+u_2\rho_y=0,
\qquad
-\Delta_{x,y}\psi=\rho_x,
\qquad
\mathbf u=(-\psi_y,\psi_x).
\]

因此 `∇·u=0`，密度满足不可压缩 Darcy 流的输运方程。Poisson 右端是密度的水平导数，不能把 `ρ` 本身直接作为右端。

物理墙满足

\[
\psi(x,0,t)=0 \quad\Longrightarrow\quad u_2(x,0,t)=0.
\]

这里没有给密度施加零值边界条件。密度可以在墙上非零，密度水平导数也可以非零。外部计算边界由有限半圆弧截断，流函数的外弧值通过半平面镜像 Green 核近似计算。

### 动态重标度与两种时间

定义正尺度 `C_l,C_ω` 和水平位移 `S`：

\[
X=C_lx+S,\qquad Y=C_ly,
\]
\[
R(X,Y,\tau)=C_\omega\rho(x,y,t),
\quad
\Psi=C_lC_\omega\psi,
\quad
\mathbf U=C_\omega\mathbf u,
\quad
\frac{dt}{d\tau}=\frac{C_\omega}{C_l}.
\]

规范速率定义为

\[
c_l=\frac{d\log C_l}{d\tau},\qquad
c_\omega=\frac{d\log C_\omega}{d\tau},\qquad
S_\tau=c_lS+c_r.
\]

固定 `X,Y` 时，物理坐标满足

\[
x_\tau=-\frac{c_lX+c_r}{C_l},\qquad
y_\tau=-\frac{c_lY}{C_l}.
\]

把链式法则代入密度方程，得到当前代码实际求解的系统：

\[
-\Delta_{X,Y}\Psi=R_X,
\qquad
\mathbf U=(-\Psi_Y,\Psi_X),
\]
\[
\boxed{R_\tau=-(U_1+c_lX+c_r)R_X-(U_2+c_lY)R_Y+c_\omega R.}
\]

尺度和物理钟同步积分：

\[
(\log C_l)_\tau=c_l,\quad
(\log C_\omega)_\tau=c_\omega,\quad
S_\tau=c_lS+c_r,\quad
t_\tau=e^{\log C_\omega-\log C_l}.
\]

`τ` 是规范时间，`t` 是物理时间，机器耗时则是墙钟秒数。三者必须分开。实验目标 `τ=200` 不表示物理时间 `t=200`。

固定网格上，离散密度右端记作 \(\mathbf F_h(\mathbf R)\)，具体定义见节点状态一节。尺度和时钟不作为 Poisson/输运的额外物理输入，仍通过上述 ODE 随同推进，用于恢复物理量。

来源：[空间右端](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+spatial/evaluate.m:1)、[共享时间内核](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm/+evolve/stepRk.m:414)。

### 极坐标方程

在重标度平面引入

\[
X=r\cos\theta,\quad Y=r\sin\theta,
\quad 0\le r\le R_b,\quad 0\le\theta\le\pi,
\quad R_b=10^{10}.
\]

水平和垂直导数为

\[
\partial_X=\cos\theta\,\partial_r-\frac{\sin\theta}{r}\partial_\theta,
\qquad
\partial_Y=\sin\theta\,\partial_r+\frac{\cos\theta}{r}\partial_\theta.
\]

因此 Poisson 源是

\[
q=R_X=\cos\theta R_r-\frac{\sin\theta}{r}R_\theta,
\]

流函数满足

\[
\boxed{-\left(\Psi_{rr}+\frac1r\Psi_r+\frac1{r^2}\Psi_{\theta\theta}\right)=q.}
\]

Darcy 速度的径向和切向分量为

\[
U_r=-\frac{\Psi_\theta}{r},\qquad U_\theta=\Psi_r.
\]

这里 `U_θ` 是切向线速度，角速度是 `U_θ/r`。

令合成输运速度为

\[
\mathbf V=(V_X,V_Y),\qquad
V_X=U_1+c_lX+c_r,\quad V_Y=U_2+c_lY.
\]

转换到极坐标后，

\[
V_r=U_r+c_lr+c_r\cos\theta,\qquad
V_\theta=U_\theta-c_r\sin\theta.
\]

最终密度方程为

\[
\boxed{
R_\tau=-\left(-\frac{\Psi_\theta}{r}+c_lr+c_r\cos\theta\right)R_r
-\frac{\Psi_r-c_r\sin\theta}{r}R_\theta+c_\omega R.
}
\]

两条墙对应 `θ=0,π`，上面都有 `Ψ=0,V_Y=0`。外半径 `R_b` 是重标度坐标中的半径；对应物理圆心为 `(-S/C_l,0)`，物理半径为 `R_b/C_l`，随尺度变化。

连续方程也可等价写成

\[
R_\tau+\frac1r\partial_r(rV_rR)+\frac1r\partial_\theta(V_\theta R)
=(2c_l+c_\omega)R.
\]

当前代码采用前面的非保守输运形式，实际反应项只有 `c_ωR`；不能再额外加入 `2c_lR`，否则会重复计算膨胀贡献。

<a id="state-grid"></a>

## 2. 实际存储状态与网格

### 连续场与节点状态

连续密度仍满足 \(R_\tau=F(R)\)。内存采用原点标量与中心数组的分解：

\[
\boxed{R_\circ(\tau)=R(0,\tau),\quad \mathbf W=\mathbf R-R_\circ\mathbf1,\quad \mathbf W_{0j}=0.}
\]

给定 \((\mathbf W,R_\circ)\)，先根据空间变化求流场与规范，再计算

\[
\mathbf G=\mathcal T_h(\mathbf W;\mathbf V)+c_\omega\mathbf W,
\qquad g_0=\mathbf G_{00},
\]
\[
\boxed{\mathbf W_\tau=\mathbf G-g_0\mathbf1,\qquad R_{\circ,\tau}=g_0+c_\omega R_\circ.}
\]

完整密度右端和空间诊断为

\[
\mathbf F_h=\mathbf W_\tau+R_{\circ,\tau}\mathbf1,
\qquad \mathbf q=\mathsf D_X\mathbf W,
\qquad \mathbf Q=\mathsf D_X\mathbf W_\tau.
\]

空间常数的导数为零，所以在精确算术中 \(\mathbf q=\mathsf D_X\mathbf R\)、\(\mathbf Q=\mathsf D_X\mathbf F_h\)，完整密度仍满足原方程。分解避免大密度常数淹没较小的空间变化；没有额外空间背景或新的物理系统。\(\mathbf q,\mathbf Q\) 不独立推进。\(\mathsf D_X\) 是源与诊断的水平微分算子，不能与 WENO 上风导数混用。

初始连续密度记为 \(R_0\)，其节点采样为 \(\mathbf R^0\)。令

\[
\zeta=(X-X_{\rm G})/s_x,
\]
\[
R_0=\frac{A_{\rm G}s_x\sqrt\pi}{2}
\left[(1-w_{\rm G})\operatorname{erfc}(-\zeta)
+w_{\rm G}\operatorname{erfc}\left(-\frac{\zeta-.8}{1.25}\right)\right]
e^{-(Y/s_y)^2},
\]
\[
\partial_XR_0=A_{\rm G}\left[(1-w_{\rm G})e^{-\zeta^2}
+\frac{w_{\rm G}}{1.25}e^{-((\zeta-.8)/1.25)^2}\right]e^{-(Y/s_y)^2}.
\]

参数为 \(A_{\rm G}=.7,X_{\rm G}=0,s_x=.75,s_y=1,w_{\rm G}=.2\)。它是 Gaussian 混合源的水平原函数，固定 \(R_0(-\infty,Y)=0\)，在 \(X\to+\infty\) 有非零平台。

来源：[初态](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm/+field/initialDensity.m:70)、[初始化](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/initialize.m:1)。

### 节点索引与唯一原点

节点统一从零编号：

\[
0=r_0<r_1<\cdots<r_{N_r-1}=R_b,\qquad
\theta_j=\frac{j\pi}{N_\theta-1},\quad j=0,\ldots,N_\theta-1.
\]

\(\mathbf R_{ij}\) 表示 \((r_i,\theta_j)\) 上的密度值。逻辑轴与径向映射为

\[
\xi_i=\frac{i}{N_r-1},\qquad
r(\xi)=a_g\xi+(R_b-a_g)
\left[\frac{\sinh(\alpha_g\xi)}{\sinh\alpha_g}\right]^3.
\]

\(a_g\) 使第一正径节点满足 \(r_1=h_0\)。配置中的首间距请求值初始为 `.03`；当前构造器按核心半径自动选择实际间距：

\[
h_0=\max\left(h_{\rm requested},\frac{r_c}{.6(N_r-1)}\right),
\qquad \xi_c=\frac{r_c}{(N_r-1)h_0}\le .6.
\]

因此核心最多占用60%的逻辑径轴，余下至少40%的逻辑区间用于外尾；核心增长时允许同步增大首间距。\(\alpha_g\) 仍按核心位置处的尾部贡献确定，初始 \(r_c=2\)、尾部比例 `.01`，外端仍准确为 \(R_b=10^{10}\)。`build` 使用构造器返回的实际 \(h_0\)，不把它改回请求值。以下初始配置仍有 \(h_0=.03\)。区间宽只记作 \(\Delta r_i=r_{i+1}-r_i\)，不与时间步共用符号。

|配置|\(N_r\)|\(N_\theta\)|\(r\le2\) 正径环数|核心区间宽|
|---|---:|---:|---:|---:|
|默认|321|257|66|`.03–.0319100`|
|旧探索初态|385|128|66|`.03–.0316047`|

385×128是节点数，角向有127个间隔。密度数组有49280个条目，但第一行全代表同一个原点，唯一物理自由度为 \(1+(385-1)128=49153\)。该行值必须完全相同。

密度导数在正径节点使用七点局部多项式模板和极坐标变换。每一行实际按局部差值

\[
(Df)_i=\sum_jD_{ij}(f_j-f_i)
\]

求和，保留原模板系数，使常平台在乘加前精确消去。仅减一个全域常数不足以保护其他区域的平台。原点用前四个正径环做六次 Cartesian 多项式拟合。常数项固定，其他27个系数通过预先构造的 QR 映射求得，不在原点计算 `1/r`。

网格中的 `Dr/Drr` 服务于密度导数及相关几何，Poisson 径向离散另用后文的 B 样条空间。角向物理弧长约为 \(r\Delta\theta\)，远区固定角间隔很宽；覆盖到 \(R_b=10^{10}\) 不代表远处薄层也充分分辨。

来源：[网格](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+mesh/build.m:1)、[径向映射](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+mesh/uniformCoreAxis.m:1)、[密度导数](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+field/derivatives.m:1)。

<a id="poisson"></a>

## 3. Poisson：一步一步从密度导数求出流场

本节按求解的依赖关系讲流程。每一步都沿用同一个问题：**已知当前密度的水平导数，求当前流函数，再由流函数求速度。** 展开、表示、组装和求解各承担一步；具体基函数、积分公式和实现参数放在[Poisson 公式补充](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_POISSON_DERIVATION.md)。

先记住未知量的两次转换：

\[
\boxed{
\text{二维未知函数 }\Psi(r,\theta)
\ \longrightarrow\
\text{一组径向未知函数 }\phi_m(r)
\ \longrightarrow\
\text{一组未知系数 }z_{am}.
}
\]

第一箭头通过角向正弦展开完成，中间包含已知墙源修正；第二箭头通过径向 B 样条表示完成。接下来建立系数方程并求解，最后沿相反方向恢复流函数。

### 第一步：明确已知源、未知解和边界

在一个固定时间阶段，密度已经给定。先对密度作水平微分，得到节点源

\[
q=R_X=D_XW.
\]

Poisson 模块拿到的是这些源节点值及实际网格，要解的是

\[
-\left(\Psi_{rr}+\frac1r\Psi_r+
\frac1{r^2}\Psi_{\theta\theta}\right)=q.
\]

还要指定边界：两条直墙上 \(\Psi=0\)；外半圆弧上 \(\Psi=g_b(\theta)\)。程序由同一份源，通过半平面镜像 Green 核估算外弧值 \(g_b\)。这一步补齐了内部求解所需的边界数据。

**这一步结束后：源和边界都已知，流函数仍未知。** Poisson 本身不推进时间。

### 第二步：先让已知修正承担墙上的非零源

墙上流函数为零，但源 \(q\) 可以非零。这里有一个与后续表示有关的问题：有限正弦和在墙上的值和角向二阶导数都为零，单靠它无法保留非零墙源所要求的墙端曲率。

当前实现据此构造一个**已知**修正函数 \(L(r,\theta)\)：它在墙上的值为零，同时让 \(-\Delta L\) 在墙上等于给定墙源。令

\[
\Psi=L+\Phi,
\]

原问题便改写为

\[
-\Delta\Phi=q+\Delta L.
\]

这样，墙源要求的曲率由 \(L\) 承担，剩余量 \(\Phi\) 再交给正弦展开。外弧条件也相应变成 \(\Phi=g_b-L\)。

构造 \(L\) 及后面的积分都需要径向区间内部的源值，程序用局部 Hermite 插值从源节点值重构这些值和导数。**插值准备的是已知源；这里尚未求出未知流函数。**

**这一步结束后：已知修正 \(L\) 已备好，接下来只需求余项 \(\Phi\)。** 墙源修正是当前离散表示的处理；一般 Dirichlet Poisson 问题的可解性并不依赖必须采用这种修正。

这只修正了第一层墙端不相容：\(f=q+\Delta L\) 的墙值为零，使余项的二阶角导数为零；四阶角导数通常仍为 \(-r^2 f_{\theta\theta}|_{\rm wall}\)，不自动消失。当前有限正弦截断仍有高阶端点误差，不能据此保证谱收敛。推导、反例及精度含义见[高阶墙面相容性](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_POISSON_DERIVATION.md#wall-compatibility)。

在足够光滑、精确角投影等条件下，这项修正对应流函数四阶、一阶角导三阶、二阶角导二阶的角向截断界；它不是当前完整程序的全域阶数保证。假设与尾部估计见[条件性精度分析](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_POISSON_DERIVATION.md#angular-accuracy)。

### 第三步：展开角方向，把二维问题拆成一组一维问题

从 Cartesian 二维方程、极坐标链式法则、墙源修正到正弦投影的逐项推导，见[二维到径向一维的完整推导](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_POISSON_DERIVATION.md#angular-derivation)。

余项在两条墙上都是零，因此沿角方向写成

\[
\Phi(r,\theta)\approx\sum_m\phi_m(r)\sin(m\theta).
\]

这里 \(\sin(m\theta)\) 是选定的角向形状；\(\phi_m(r)\) 表示这种形状在不同半径上的幅度，仍然是待求函数。

将展开代回方程，利用
\(\partial_{\theta\theta}\sin(m\theta)=-m^2\sin(m\theta)\)，角向微分就变成乘以已知数 \(-m^2\)。对已知源和修正项使用相同的离散正弦变换后，每个角模分别满足

\[
-\phi_m''-\frac1r\phi_m'+\frac{m^2}{r^2}\phi_m=f_m(r).
\]

\(f_m\) 是该角模已知的源，包含墙源修正。外端值 \(\gamma_m\) 来自 \(g_b-L\) 的正弦变换；原点值固定为零。不同角模在这里互不耦合。

**这一步结束后：一个二维 Poisson 问题变成了一组一维径向边值问题。** 程序用 FFT 实现正弦变换；径向问题还没有解。

### 第四步：用 B 样条表示每个径向未知函数

对某个固定角模，未知仍是整条曲线 \(\phi_m(r)\)。程序预先选好一组由实际径向网格确定的 B 样条基 \(N_a(r)\)，写成

\[
\phi_{m,h}(r)=\sum_a z_{am}N_a(r).
\]

基函数是已知形状，\(z_{am}\) 是未知系数。它们通常不等于径向节点上的函数值。

边界条件先固定两端系数：原点对应零值，外端对应已知的 \(\gamma_m\)。剩余内部系数才需要求解。

**这一步结束后：求一条未知函数的问题，变成了求有限个系数的问题。** 这一步只确定解的表示方式，还没有确定系数应取多少。

### 第五步：用 Galerkin 条件建立系数方程

要确定系数，需要把微分方程转成代数条件。程序将径向方程乘以测试函数和极坐标权 \(r\)，沿半径积分，再分部积分。对两端为零的测试函数 \(v\)，得到

\[
\int r\,v'\phi_{m,h}'\,dr
+m^2\int\frac{v\phi_{m,h}}r\,dr
=\int r\,v f_m\,dr.
\]

取各个内部 B 样条作为测试函数，再代入上一步的样条表示，每个测试函数给出一条关于系数的线性方程。合在一起便是

\[
\boxed{\mathbf A_m\mathbf z_m=\mathbf b_m,
\qquad \mathbf A_m=\mathbf K+m^2\mathbf M.}
\]

\(\mathbf K\) 来自径向导数的积分，\(\mathbf M\) 来自角向模态项的积分；已知源、修正项和已移到右端的外边界贡献组成 \(\mathbf b_m\)。实际组装时，墙源修正的径向二阶项也经过分部积分。

**这一步结束后：每个径向边值问题都变成一个线性系统。** B 样条决定用什么函数表示解，Galerkin 决定这些系数必须满足什么方程。

### 第六步：解线性系统，得到真正的未知系数

现在才进入线性代数求解：对每个角模，用直接 LU 分解和回代求出 \(\mathbf z_m\)。

固定网格时，基函数及 \(\mathbf A_m\) 不变，因此矩阵分解可以缓存。进入新的时间阶段，源和边界改变，程序更新右端 \(\mathbf b_m\)，复用分解求新的系数。换网格后则重建相应矩阵与分解。

**这一步结束后：各个角模的径向样条系数已求出。** 缓存复用的是矩阵分解，每个新源仍需求出自己的解。

### 第七步：按相反顺序恢复流函数，再求速度

先把系数代回样条展开，得到各个径向模态；再合成角向正弦展开，最后加回第二步的已知修正：

\[
\Psi_h(r,\theta)
=L(r,\theta)+
\sum_m\left[\sum_a z_{am}N_a(r)\right]\sin(m\theta).
\]

这就是前面各步共同求出的流函数。在正半径处，速度由**同一个重构函数**求导：

\[
U_r=-\frac1r\partial_\theta\Psi_h,\qquad
U_\theta=\partial_r\Psi_h.
\]

径向求导作用于 B 样条，角向求导作用于正弦函数，并加上 \(L\) 的相应导数。再将径向、切向分量转回 Cartesian 速度。原点由第一角模的径向斜率指定一个 Cartesian 速度值；当前空间尚未消除高模态原点斜率缺陷，具体约定见[公式补充](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_POISSON_DERIVATION.md)。

**这一步结束后：当前阶段的流函数和 Darcy 速度已得到，可交给规范与输运模块。**

### 第八步：分别检查线性求解和空间方程

程序做两种检查：

- 将系数代回 \(\mathbf A_m\mathbf z_m-\mathbf b_m\)，检查组装出的线性系统是否解好。
- 在额外径向审计点上，用重构函数的导数计算 \(-\Delta\Psi_h-\widetilde q_h\)，检查离散流函数对重构源的 PDE 缺陷。

第一种误差很小，仍可能存在第二种空间缺陷；两者也都不能替代不同网格之间的收敛检验。具体定义见[空间审计](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ALGORITHM_DETAILS.md#audits)。

|流程中的动作|它解决的问题|完成后仍待求什么|
|---|---|---|
|准备源、外弧值与墙源修正|给出完整的已知数据|余项流函数|
|角向正弦展开|拆开二维问题的角方向|各个径向模态函数|
|径向 B 样条表示|把未知函数化为有限个系数|样条系数|
|Galerkin 组装|建立系数必须满足的方程|线性系统的解|
|LU 回代|求出系数|只需重构与求导|
|样条求值、正弦合成、加回修正|恢复当前流函数|由其求导得到速度|

这里按理解问题的依赖关系展开。实际代码将固定几何预处理集中缓存，并把多列源重构批量计算；完整变换、墙源修正、弱负载与边界消元公式见[公式补充](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_POISSON_DERIVATION.md)。

来源：[Poisson 完整实现](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+spatial/poisson.m:1)。


<a id="gauge"></a>

## 4. 规范：由墙观测求坐标速率

本模块输入中心节点密度 \(\mathbf W\) 和已经求出的 Darcy 速度，求三个标量 \(c_l,c_\omega,c_r\)。它通过固定墙观测确定坐标缩放和平移速率。

令 \(\mathcal I_w\) 表示：把两条极坐标墙射线拼成带符号的 \(X\) 轴，再用三次 `spline` 插值到原始固定 Cartesian 墙节点。令 \(\mathsf D_w\) 为这套固定墙节点上的原始微分矩阵，则

\[
\mathbf W^w=\mathcal I_w\mathbf W,
\qquad \mathbf q^w=\mathsf D_w\mathcal I_w\mathbf W.
\]

墙观测使用 \(\mathbf q^w\)，Poisson 使用 \(\mathbf q=\mathsf D_X\mathbf W\)；二者来自不同离散算子。

原墙 \(\mathsf D_w\) 也按局部差值 \(\sum_\mu(D_w)_{\nu\mu}(W_\mu^w-W_\nu^w)\) 计算。常数 \(R_\circ\) 不进入墙梯度；节点、模板、I/P/C 定义和原参考全部保留。

以下用 \(\nu\) 编号固定墙节点，\(q_\nu^w\) 表示 \(\mathbf q^w\) 的分量。三个离散观测函数定义为

\[
I(\mathbf R)=\sum_\nu\omega_\nu^I q_\nu^w,
\qquad
P(\mathbf R)=\max_{X\in\text{原始墙域}}
\mathcal S_{\rm natural}[\mathbf q^w](X).
\]

其中 \(\omega_\nu^I\) 是固定区间 \([-1,1]\) 上、对分段线性插值准确的裁剪梯形权重；\(\mathcal S_{\rm natural}\) 是自然三次样条。计算质心时，先在线性插值中加入 \(X=\pm1\) 端点并裁剪，得到 \((X_\nu^c,q_\nu^c)\) 及正梯形权重 \(\omega_\nu^c\)，再取

\[
E=\sum_\nu\omega_\nu^c(q_\nu^c)^2,
\qquad
C(\mathbf R)=\frac{\sum_\nu\omega_\nu^cX_\nu^c(q_\nu^c)^2}{E}.
\]

下文 \(I,P,C\) 是这些函数在当前 \(\mathbf R\) 上的值。对任意节点变化率 \(\mathbf G\)，记方向导数

\[
\dot I[\mathbf G]
=\left.\frac{d}{d\varepsilon}I(\mathbf R+\varepsilon\mathbf G)\right|_{\varepsilon=0},
\]

\(\dot P[\mathbf G]\)、\(\dot C[\mathbf G]\) 同理。例如，将 \(\mathsf D_w\mathcal I_w\mathbf G\) 按相同规则裁剪为 \((\mathsf D_w\mathcal I_w\mathbf G)_\nu^c\)，有

\[
\dot C[\mathbf G]
=\frac{2\sum_\nu\omega_\nu^c(X_\nu^c-C)q_\nu^c
(\mathsf D_w\mathcal I_w\mathbf G)_\nu^c}{E}.
\]

固定参考为

\[
I_0=.8340921004503201,\quad
P_0=.6369025096521013,\quad
C_0=.06283125189791006.
\]

需要满足的瞬时约束是

\[
\dot I[\mathbf F_h]=0,\qquad
\dot P[\mathbf F_h]=0,\qquad
\dot C[\mathbf F_h]=-\kappa(C-C_0),\qquad \kappa=.5.
\]

它约束观测的变化率，不会每步投影到 \(I=I_0,P=P_0\)，参考也不随换网格改变。

具体求解时，固定 \(U_{\rm ref}=U_1(C,0)\)，引入中心化平移参数 \(\delta=c_r+U_{\rm ref}\)。给定试探的 \((c_l,\delta)\)，令 \(c_r=\delta-U_{\rm ref}\)，据此构造输运速度 \(\mathbf V\)，调用下一节的无反应输运算子

\[
\mathbf T=\mathcal T_h(\mathbf W;\mathbf V).
\]

只需求两个方程：

\[
\boxed{
\frac{\dot I[\mathbf T]}{I}
-\frac{\dot P[\mathbf T]}{P}=0,
\qquad
\dot C[\mathbf T]+\kappa(C-C_0)=0.
}
\]

然后恢复

\[
c_\omega=-\frac{\dot P[\mathbf T]}P,
\qquad c_r=\delta-U_{\rm ref},
\qquad \mathbf G=\mathbf T+c_\omega\mathbf W.
\]

这一消元利用了墙采样对数据线性，以及幅度缩放使 \(I,P\) 同比例变化、使 \(C\) 不变。\(\mathbf G\) 与完整 \(\mathbf F_h\) 只差空间常数，所以对实际 \(\mathbf G\) 的梯度做独立审计，就在检查完整密度的三个约束。

两变量系统使用数值 Jacobian、Newton 和阻尼；上风方向随速度符号切换，因此它一般是分片非线性的。探测只重复 WENO 输运和墙观测，复用当前 Darcy 场，不重复 Poisson。根默认容差为 `10⁻⁶`，完整 RHS 的规范审计为 `10⁻⁵`。峰值有并列分支时，\(\dot P\) 采用代码选定最大分支的评价导数。

同一 RHS 的规范探测绑定当前 \(\mathbf W\) 的观测闭包，复用固定墙插值、I/P/C、峰位置和方向导数权重；后续仅采样试探 RHS 并求墙导数，不跨状态持久缓存。每个 RHS 仍有4次完整状态观测：空间入口、规范入口、重复回调审计、最终完整 RHS 审计。后两次从实际数组重新计算，保持双审计独立性；根的种子、容差和求解逻辑未删减。

来源：[规范测量](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+gauge/measure.m:1)、[规范求解](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+gauge/solve.m:1)。

<a id="transport"></a>

## 5. WENO：计算输运变化率

本模块输入 \(\mathbf W\)、Darcy 速度与规范速率，求无反应输运 \(\mathcal T_h(\mathbf W;\mathbf V)\)，再形成 \(\mathbf G\)、\(g_0\)、\(\mathbf W_\tau\) 与 \(R_{\circ,\tau}\)。中心数组已有唯一零原点。

径向 WENO 在均匀逻辑坐标 \(\xi\) 上计算，故使用逻辑速度 \(V_r/J_r\)，其中 \(J_r=dr/d\xi\) 在实现中由实际网格度量给出。角向 WENO 在均匀 \(\theta\) 轴上计算，使用角速度 \(V_\theta/r\)。记 \(D_\xi^+,D_\theta^+\) 为正速度所用的左侧上风导数，\(D_\xi^-,D_\theta^-\) 为负速度所用的右侧上风导数；上标标记速度方向，不表示导数本身的正负。

在正径向内部节点，无反应输运为

\[
\begin{aligned}
\mathcal T_h(\mathbf W;\mathbf V)={}&
-\max\!\left(\frac{V_r}{J_r},0\right)D_\xi^+\mathbf W
-\min\!\left(\frac{V_r}{J_r},0\right)D_\xi^-\mathbf W\\
&-\max\!\left(\frac{V_\theta}{r},0\right)D_\theta^+\mathbf W
-\min\!\left(\frac{V_\theta}{r},0\right)D_\theta^-\mathbf W.
\end{aligned}
\]

各乘除逐节点进行。生产配置为 `wenoEpsilon=10⁻¹²,alphaMode=local,wenoInflowDensityDegree=1`。这里的 `wenoEpsilon` 是归一化权重正则量的下限；当前 local 模式实际采用

\[
\epsilon_{\rm eff}=\max(\epsilon_{\rm cfg},1/512^2)
=3.814697265625\times10^{-6}.
\]

专用 MEX 与共享参考路径使用相同规则。生产配置将物理入流密度外推次数设为1；直接调用输运函数而不传该配置时，底层回退值为5。物理外边界的入流外推与下面的原点坐标延拓是不同处理。

两条墙拼成一条带符号的水平轴，另做 Cartesian WENO 墙输运，再拆回两端列；唯一原点的墙输运值赋给整行。因此上式的 \(1/r\) 不在原点求值。外边界采用外推及相应入流重构，密度不固定为零。

加上反应项得 \(\mathbf G=\mathcal T_h(\mathbf W;\mathbf V)+c_\omega\mathbf W\)，然后执行前述原点分离，按需恢复完整 \(\mathbf F_h\)。

固定本次 \(\mathbf W\) 后，上风导数只预备一次，规范探测仅按速度重新组装。专用 `weno5UnitPairMex` 同时求单位速度 \(+1,-1\)、单位度量的两套导数；实际速度和网格度量仍按原公式在外部施加。核保留原 WENO 权重归一化、六点幅度尺度、epsilon 规则和逐行入流外推次数。左右 face 使用共同参考值相减，ghost 的小偏移单独保留，使大常平台精确消去。它不改变 CFL 或入流规则；无 MEX 或显式关闭时仍可使用共享参考实现。

这里 WENO 的上风导数服务于输运；\(\mathsf D_X\) 是前文用于密度微分和 Poisson 源的算子，两者不混用。一次完整右端评估先做密度微分与 Poisson，再在固定 Darcy 场下耦合求规范与输运，最后作原点分离，得到 \(\mathbf W_\tau,R_{\circ,\tau}\) 及 \(\mathbf Q\)、CFL 和审计。

来源：[输运实现](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+field/transport.m:1)、[完整 RHS](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+spatial/evaluate.m:8)。

### 原点附近的负半径延拓

对每条径向射线，首圈 WENO 模板需要原点另一侧的值。负的有符号半径对应下半平面坐标，不能直接将密度对 \(Y\) 作偶延拓。当前代码在径向数组前补三个样本，分别对应 \(-r_3,-r_2,-r_1\)。令 \(\mathcal P_\circ(X,Y)\) 为已有原点 Cartesian 多项式拟合，定义其奇 \(Y\) 部分

\[
\mathcal P_{\circ,\mathrm{odd}\,Y}(X,Y)=\tfrac12[\mathcal P_\circ(X,Y)-\mathcal P_\circ(X,-Y)].
\]

在正半径 \(s\in\{r_1,r_2,r_3\}\) 处，构造

\[
\boxed{
W_{\rm ghost}(-s,\theta)
=W(s,\pi-\theta)
-2\mathcal P_{\circ,\mathrm{odd}\,Y}(-s\cos\theta,s\sin\theta).
}
\]

右端第一项读取上半圆的实际节点值；第二项修正穿过 \(Y=0\) 时的奇部分。若场恰为拟合空间中的 Cartesian 多项式，该式给出其正确的负半径值。两条墙的 \(Y=0\)，奇部分为零，因此恰好使用对侧实际墙数据。生产采用六次拟合，底层也支持四次拟合。

补入样本后再计算原 WENO，最后只保留原生节点的导数。这使首圈内角点的径向模板与穿过原点的墙模板衔接；角向和拼接墙算子沿用原实现。这里没有过滤 \(W\)、\(F\) 或 \(Q\)，也没有修改原点单值条件。该延拓修复不等于 Poisson 高模原点斜率约束，后者仍有独立的待验证问题。

来源：[延拓及输运](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+field/transport.m:173)、[MEX 权重规则](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+field/weno5_unit_pair_mex.cpp:95)、[共享参考权重规则](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm/+field/weno5FluxDerivative.m:151)。

<a id="time"></a>

## 6. 时间：直接 CFL 与 SSPRK3

本模块实际推进 \((\mathbf W,R_\circ)\) 与尺度 ODE。每个 SSPRK3 阶段使用对应的中心数组和标量；由于 \(\mathbf R=\mathbf W+R_\circ\mathbf1\) 是线性变换，在精确算术中与完整密度的 SSPRK3 等价，浮点差异由真实短轨道验收。

对正径节点，CFL 输运率为

\[
\lambda_{ij}=\frac{|V_r(r_i,\theta_j)|}{J_{r,i}/(N_r-1)}
+\frac{|V_\theta(r_i,\theta_j)|}{r_i\,\pi/(N_\theta-1)},
\qquad i\ge1.
\]

\(J_{r,i}\) 是实际映射 \(r(\xi)\) 的离散度量，不直接等于某个相邻节点差。两墙改用拼接水平轴上的 \(|V_X|/\Delta X_{\rm metric}\)，原点取该轴唯一原点的速率。把这些节点率共同取最大，再计入幅度反应，得到全局率 \(\Lambda\) 及时间步

\[
\Lambda=\max_{i,j}\lambda_{ij}+|c_\omega|,\qquad
\boxed{\Delta\tau_n=\min\left(.2,\frac{.4}{\Lambda},
\tau_{\rm end}-\tau_n,
\frac{t_{\rm end}-t_n}{C_\omega/C_l}\right).}
\]

式中各场和尺度均取步初接受态。当前实验 \(t_{\rm end}=\infty\)；`minDt=10⁻¹⁰`，低于它时停止。代码记录终点 CFL，但没有对所有内部阶段和终点重新选步；因此 `.4` 是步初规则，不是全部阶段已有的稳定性保证。

为了与连续方程对照，线性重构后的完整密度满足以下 SSPRK3 形式；实现对 `W` 和 \(R_\circ\) 分别使用相同阶段组合：

\[
\mathbf R^{[1]}=\mathbf R^n+\Delta\tau_n\mathbf F_h(\mathbf R^n),
\]
\[
\mathbf R^{[2]}=\tfrac34\mathbf R^n+
\tfrac14\left[\mathbf R^{[1]}+\Delta\tau_n\mathbf F_h(\mathbf R^{[1]})\right],
\]
\[
\boxed{\mathbf R^{n+1}=\tfrac13\mathbf R^n+
\tfrac23\left[\mathbf R^{[2]}+\Delta\tau_n\mathbf F_h(\mathbf R^{[2]})\right].}
\]

方括号上标表示本步内部阶段，\(n\) 表示接受步编号。\(R_{\circ,\tau}=g_0+c_\omega R_\circ\) 与中心数组同阶段积分，尺度和位移也使用相同组合；新端点再评估一次完整右端，用于检查、诊断和下一步缓存。

物理钟保存为 `physicalTime` 高位和 `physicalTimeLow` 低位。SSPRK3 的物理时间增量经补偿加法累计，先把低于当前高位分辨率的小增量保存在低位；没有重定义物理时间，也没有修改规范时间、CFL 或时间步。

固定网格且步初缓存命中时，每步需要3次新的完整 RHS/Poisson：\(\mathbf R^{[1]}\)、\(\mathbf R^{[2]}\) 和新端点上的评估。步初右端复用上一接受态缓存；无缓存时为4次。换网格候选另需完整右端。

当前推进没有隐式时间阶段、GMRES/ILU、整步对两半步误差控制或失败步长重试。数值阶段失败时保存最后接受态并停止。规范小系统阻尼、候选网格阻尼和分段恢复各有独立用途。

来源：[CFL 速率](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+spatial/transportRate.m:1)、[单步入口](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+evolve/advance.m:1)。

<a id="remesh"></a>

## 7. 换网格：同一时刻的表示转移

本模块已知同一接受时刻的 \(\mathbf W,R_\circ\)，把中心数组转移到新径向节点上；尺度、位移和两个时钟保持该时刻的值。它不推进时间。候选网格上重新计算 \(\mathbf q=\mathsf D_X\mathbf W\)、流场和中心 RHS，再恢复完整诊断用于审核转移。

对 \(\mathbf q\) 的正峰，在 `.1,.2,.5` 倍峰值处提取超水平集。面积通过原生节点构成的 \((X,Y)\) 三角形中的 P1 等值线裁剪求得，包含所有分离分量；原点使用唯一顶点扇形。

当前核心半径优先覆盖高峰区域。令 \(r_\alpha\) 为 \(\alpha\) 倍正峰超水平集的完整径向最外界，包含所有分离分量；配置 `uniformCorePadding=2` 给出

\[
r_c=\max\left(r_{c,\min},\min(r_{.1},2r_{.5})\right),
\qquad r_{c,\min}=\texttt{grid.coreRadius}.
\]

其中 \(r_{.1}\) 是 `.1` 水平集的最外径，不是将某个半径乘 `.1`。高峰区采用近均匀核心，较弱而更远的完整源尾部交给拉伸网格表示；没有从 Poisson 负载或输运中裁剪源。节点数、均匀角轴和 \(R_b=10^{10}\) 保留，首间距请求下限由配置给定（默认 `.03`，冻结06实验为 `.06`）。实际目标首间距仍满足

\[
h_{\rm target}=\max\left(h_{\rm requested},\frac{r_c}{.6(N_r-1)}\right),
\]

其中 \(h_{\rm requested}\) 取配置首间距与换网格首间距下限的较大值，保留至少40%的逻辑区间作为尾部。

换网格请求有两个独立来源：三个完整水平集面积中任一个相对上次成功换网格参考变化达到20%，记 `areaRequested`；若现网格首间距满足 \(r_1^{\rm old}>h_{\rm target}(1+10^{-8})\)，记 `resolutionRequested`。两者取逻辑“或”。因此即使面积尚未改变20%，也能主动细化由旧核心规则造成的过粗首间距；三个面积触发仍使用完整源，并未缩为最高水平集。

候选处理依次为：

1. 检查径轴单调、端点、首间距及相邻网格尺度比。
2. 在同一接受时刻转移 \(\mathbf W\)，标量 \(R_\circ\) 保持不变。共享五次插值的每个模板先减局部参考、插值后加回，以保常平台；不先形成大完整密度再插值。
3. 重建实际几何和空间算子，保持尺度、时钟和原始墙参考。
4. 复用转移前后的固定墙观测，检查离散积分与墙积分的相对跳变，默认上限均为 \(10^{-3}\)；墙梯度峰另要求 \( |P_{\rm new}-P_{\rm old}|/\max(|P_{\rm old}|,\mathrm{realmin})\le10^{-4}\)，配置为 `peakTransferTolerance`。三个门槛共同限制同一时刻的网格转移。
5. 重新评估 Poisson、规范、输运和 \(\mathbf Q=\mathsf D_X\mathbf F_h(\mathbf R)\)，要求全域缩放强残差不超过 `.02`。
6. 把旧 Darcy 速度和 \(\mathbf Q\) 插值到候选新坐标，与重新求出的场比较；两个缩放差分别不超过 `.05` 和 `.1` 才提交。速度差使用矢量模的最大值，\(\mathbf Q\) 差使用最大绝对值；分母分别为 \(\max(1,\|\mathbf U_{\rm old}\|_\infty)\) 和 \(\max(1,\|\mathbf Q_{\rm old}\|_\infty)\)。比较的是相同新坐标上的场，不是新旧数组同下标的数值。

候选几何可按 `old+2^(−j)(target−old)` 阻尼；`minimumImprovement=.02` 时通常最多尝试6个阻尼量。候选失败回滚几何和场，已经接受的时间端点保留；面积请求保持待处理，并按8、16、32、64、128次检查退避，明显特征变化可提前重新提出。成功提交且 `output.verbose=true` 时，在既有 MATLAB 日志打印一行时刻、步号、累计换网格数、前后 I/P、峰相对跳变和强残差，便于定位转移；没有增加时间步重试。

用于检查转移跳变的离散积分为

\[
\mathcal M_h=\sum_{i,j}(R_\circ+\mathbf W_{ij})\,r_iw_i^rw_j^\theta.
\]

它属于有限计算域和当前求积规则。重标度、开放外边界、非紧支撑初态及远区分辨率都影响其物理解读；它不是无限域物理质量守恒证书。该非负初态在无限半平面的总密度积分本就不有限，远区固定角间距也可能无法解析墙附近薄层。

旧版本固定 \(h_0=.03\)，因而对385节点存在 `11.52` 的核心半径上限。当前版本在核心增长时按前述规则增大 \(h_0\)，该固定上限已消除，同时保留至少40%的逻辑区间作为尾部。仍须满足 \((N_r-1)h_0<R_b\)、\(r_c<R_b\) 以及构造器的单调与可表示性条件；它并非任意大核心都可无条件接受。

全域强残差与新坐标场差是针对旧长实验暴露问题增加的空间缺陷门槛。它们保留原有离散积分、墙积分、弱线性求解和规范约束检查，但不构成长期精度或网格收敛证书。

来源：[面积诊断](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+diagnostics/sourceFeatures.m:1)、[换网格事务](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+remesh/apply.m:1)。

<a id="audits"></a>

## 8. 空间审计：求解误差与离散缺陷

线性系统检查的是实际组装矩阵：

\[
\mathbf e_m=\mathbf A_m\mathbf z_m-\mathbf b_m,\qquad
\eta=\max_m\frac{\|\mathbf e_m\|_\infty}
{\|\,|\mathbf A_m||\mathbf z_m|+|\mathbf b_m|\,\|_\infty}.
\]

硬门槛为 `η≤10⁻¹⁰`，还记录逐分量 backward error 和弱残差。

另在每个径向区间中点、完整原生角轴上，用 B 样条与提升的解析二阶导数检查

\[
\mathcal R_{\rm strong}=-\Delta\Psi_h-\widetilde q_h.
\]

其中 \(\widetilde q_h(r,\theta_j)=\sum_m\widetilde q_m(r)\sin(m\theta_j)\) 用于内部角节点，两墙分别取 \(\widetilde q_+(r),\widetilde q_-(r)\)。它与 Poisson 负载采用同一局部 Hermite 源重构。记录全域、\(r\le2\)、首4个径向中点等区域，按各区域 \(\max(1,\|\widetilde q_h\|_\infty)\) 缩放。当前默认硬门槛为

\[
\frac{\|\mathcal R_{\rm strong}\|_\infty}
{\max(1,\|\widetilde q_h\|_\infty)}\le .02
\quad\text{（全域审计点）}.
\]

超过门槛即停止该空间评估，不作失败时间步重试。这是采样的粗空间缺陷检测，与代数残差不同，也不保证未采样位置或完整轨道达到给定精度。

代数门检查“组装出来的系统是否解好”；强残差检查“重构场在所选审计点是否近似满足 PDE”。换网格的跳变门检查同一时刻表示转移是否改变观测。空间收敛则需要不同网格下的独立对照，这些门槛不能相互替代。

来源：[Poisson 审计实现](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+spatial/poisson.m:1)。历史误差数值见[验证记录](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_VALIDATION_NOTES.md)。

<a id="data"></a>

## 9. 数据：字段、物理量恢复与绘图

### 数据结构与物理量恢复

`data.mat` 保存逐步标量轨道和切片目录，大场保存在 `data_frames`；相同网格按实际 `r,theta` 的 SHA 身份只在 `grids` 存一次。切片使用相对路径，整体搬移后可继续读取。每个场必须与该帧的实际网格、尺度和时钟配套。

每个接受步保存20列标量：

|类别|代码字段与含义|
|---|---|
|步与时钟|`step,tau,t,H` 对应接受态编号、\(\tau,t\) 和刚完成的步长；`step=n+1` 行的 `H` 为 \(\Delta\tau_n\)，初态行为 NaN；`dt` 是 `H` 的读取别名|
|变化率范数|`Qinf,Finf` 对应 \(\|\mathbf Q\|_\infty,\|\mathbf F_h(\mathbf R)\|_\infty\)|
|规范速率|`cl,cw,cr` 对应 \(c_l,c_\omega,c_r\)；读取别名 `c_l,c_omega,c_r` 分别对应三列，`c=[cl,cw,cr]` 是三列组成的数组|
|场与审计|`rhoInf,poisson,gauge,remesh`：密度最大模、Poisson 审计、规范审计、累计换网格数|
|尺度与平移|`logC_l,logC_omega,X_shift` 对应 \(\log C_l,\log C_\omega,S\)|
|导数和速度|`rhoXinf,rhoYinf,uInf`：两个密度导数和 Darcy 速度的最大模|
|离散积分|`mass` 对应 \(\mathcal M_h\)|

读取时派生 `C_l=exp(logC_l),C_omega=exp(logC_omega),C_s=cl/(cl−cw)`；最后一个分母接近零时返回 NaN。

每帧保存以下节点场，以及标量时间、尺度、速率、实际网格身份和诊断。当前 `diagnostics` 另保存完整 `poissonSolveInfo`、CFL 极值位置 `cflPosition` 和规范根条件数 `gaugeRootCondition`；运行 `report.json` 同时记录全域及核心强残差、速度最大模和上述两项定位诊断。接口别名不代表新增未知量：

|保存字段|统一数学记号与定义|
|---|---|
|`rho`|兼容完整密度 \(\mathbf R=\mathbf W+R_\circ\mathbf1\)|
|`rho_centered,rho_origin`|权威状态 \(\mathbf W,R_\circ\)|
|`rho_tau_centered,rho_origin_tau`|\(\mathbf W_\tau,R_{\circ,\tau}\)|
|`transportOrigin`|\(g_0\)，用于稳定恢复物理时间导数|
|`rho_X`|Poisson 源 \(\mathbf q=\mathsf D_X\mathbf R\)|
|`rho_Y`|密度的垂向离散导数|
|`rho_tau`|完整离散右端 \(\mathbf F_h(\mathbf R)\)；读取别名 `F`|
|`Q`|\(\mathbf Q=\mathsf D_X\mathbf F_h(\mathbf R)\)；读取兼容名 `rho_xt`|
|`psi`|样条重构 \(\Psi_h\) 在原生节点上的值|
|`u1,u2`|由 \(\Psi_h\) 求导得到的重标度 Darcy 速度 \(U_1,U_2\)|
|`vx,vy`|合成输运速度 \(V_X,V_Y\)，含缩放和平移贡献|

当前 schema4 帧保存中心主态及两者变化率，完整 `rho/rho_tau` 是兼容展示；空间导数从中心数组稳定计算。检查点采用严格 v4 kind，不能静默读取旧检查点为新状态。

已停止的04试算目录保留父03全部历史标量、旧不可变帧和实际网格，通过 `meta.mixedFrameSchemas` 显式声明 schema3/schema4 混合。读取器按每帧实际 kind/schema 解释；选旧帧时清除中心状态字段，`selected.schemaVersion` 标明实际版本。迁移输出继续保留父检查点和迁移元数据；正式06从初态开始保存 schema4，不继承旧试算前缀。`t` 列是物理钟高位，同帧 `scale.physicalTimeLow` 是低位。

读取时另给出 `vmag=hypot(u1,u2)` 和 `transport_speed=hypot(vx,vy)`。兼容名 `rho_xt` 仍指 \(\mathbf Q\)，不表示物理 \(\rho_{xt}\)。

```matlab
addpath('/Users/cbjin/Documents/one_sided_IPM/time_evolution');
d = ipm_c.data('/Users/cbjin/Documents/one_sided_IPM/runs/ipm_c_tau200_fresh_v4_20261021_06/native');
d = ipm_c.data(d,4.5,'nearest');
% d.tau / d.t 是完整轨道向量。
% d.selected.tau / d.selected.t 是当前场切片的标量时间。
```

还支持 `exact` 和按帧 `index` 选择。选帧不做时间插值；跨过输出周期后保存的是实际接受端点，例如 `.50066456`，不是强行插值到 `.5`。

连续物理量使用同一切片的尺度恢复：

\[
x=(X-S)/C_l,\qquad y=Y/C_l,\qquad \rho=R/C_\omega,
\]
\[
\rho_x=(C_l/C_\omega)R_X,\qquad
\rho_y=(C_l/C_\omega)R_Y,\qquad
\mathbf u=\mathbf U/C_\omega,\qquad
\psi=\Psi/(C_lC_\omega).
\]

物理 Eulerian 时间导数满足

\[
\rho_t=\frac{C_l}{C_\omega^2}
\left[F(R)+(c_lX+c_r)R_X+c_lYR_Y-c_\omega R\right].
\]

这里 \(F(R)=R_\tau\) 是连续方程右端。存储及绘图接口以对应帧的 \(\mathbf F_h(\mathbf R)\) 和离散密度导数近似该公式，不能把 `rho_tau` 直接当作物理 \(\rho_t\)。

对 schema4 帧，实际使用已消去大常数项的等价恢复式

\[
\boxed{\rho_t\approx\frac{C_l}{C_\omega^2}
\left[\mathbf W_\tau+g_0\mathbf1-c_\omega\mathbf W
+(c_lX+c_r)\mathbf q+c_lY\mathsf D_Y\mathbf W\right].}
\]

它使用同一帧的中心 RHS 与导数，再重采样到方形区域，避免先形成大 `rho_tau` 再减 `c_omega*rho`。旧 schema3 帧仍按原完整场公式读取。

在速率趋于常数、重标度场趋于稳定且 \(c_l>c_\omega\) 等假设下，\(dt/d\tau\) 可能衰减，物理钟可能趋向有限值，`C_s` 才具有相应尺度指数的解释。单个时刻的 `C_s` 或有限一段拟合不能证明物理奇性。

来源：[帧与轨道输出](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+output/frame.m:1)、[数据接口](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/data.m:1)。

### 方形区域插值与绘图

本模块读取一帧节点场，在指定方形或矩形区域生成 Cartesian 采样结构体。schema4 的 `rho` 先插 `W` 再加 \(R_\circ\)，`rho_tau` 先插 `W_tau` 再加 \(R_{\circ,\tau}\)，保留原三角有限元方法。它不重新求 Poisson，也不推进时间：

```matlab
opts = struct('fields',{{'rho','rho_X','psi','u1','u2'}}, ...
              'coordinates','native','method','cubic');
cartesianSlice = ipm_c.cartesian(d,[-2 2 0 4],401,opts);
graphics = ipm_c.plot(cartesianSlice,'rho_X');
exportgraphics(graphics.figure,'rho_X.pdf','ContentType','vector');
```

区域参数是 `[xmin xmax ymin ymax]`，分辨率可为单个数或 `[nx ny]`。`cartesianSlice` 包含目标轴 `x,y`、二维坐标 `X,Y`、有效域 `mask`、所选场和方法元信息；坐标含义由 `coordinates` 指定。

绘图三角网由实际极坐标节点构成，原点使用单顶点扇形。`method='cubic'` 采用三角形上的 P3 Hermite/Bernstein 有限元：共享顶点值与梯度确定边系数，局部三次最小二乘重构的质心值确定内部系数。条件不足时只降低内部质心重构阶数，并记录次数。

该绘图重构在三角形间为 \(C^0\) 连续，不是 \(C^1\)，也不保证不超调；`method='linear'` 使用 P1，采样值不超出所在三角形的顶点值范围。它们与 Poisson 余项的径向 \(C^2\) 三次 B 样条分工不同：Poisson 从源求流函数，这里仅把已保存的节点场重采样到绘图网格；完整 Poisson 流函数含 \(C^1\) 提升，不能笼统称为全局 \(C^2\)。

墙以下、外半径之外，以及最外圆弧与直边三角网之间的窄区域均返回 NaN，不外推。

`coordinates='native'` 使用 \((X,Y)\) 和重标度场；`coordinates='physical'` 使用 \((x,y)\)，并按上一节公式转换密度、密度导数、Darcy 速度和流函数。两种选择下，`rho_tau/F`、`Q/rho_xt`、`vx/vy` 均保留重标度定义，元信息逐项注明。物理 \(\rho_t\) 必须显式请求字段 `rho_t`，并采用物理坐标选项。

来源：[Cartesian 重采样](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/cartesian.m:1)、[绘图](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/plot.m:1)。

<a id="implementation"></a>

## 10. 实现：检查点恢复与复杂度

### 检查点与恢复

检查点使用两槽轮换。每次写未被当前指针指向的槽，读回核验主载荷，再重命名、计算 SHA 并原子发布 `latest.json`。主状态、实际网格、规范参考和轨道保留，派生 RHS、空间算子和 LU 在恢复时重建。

源码清单绑定 `+ipm` 与 `+ipm_c`；实际换了数值代码后不能静默继续旧轨道。恢复只允许直接修改目标规范时间、单段墙钟预算和总步数上限。检查点指针和配置含绝对路径，因此数据目录可迁移读取，不代表检查点可直接搬移恢复。

物理钟低位 `physicalTimeLow` 用于接受端点的补偿累计，不是独立的 RK ODE 分量。恢复以检查点中保存的主态及历史为准；`flow`、`rhsCache` 和时间步诊断均重新构造。

### 工作量


设 `N=N_rN_θ`，三次次数固定。固定网格时，角向 FFT 约为 `O(N_rN_θ log N_θ)`；局部带宽的逐模径向回代、五点源斜率与 Hermite 稀疏映射约为 `O(N)`。首次几何组装和各模态分解也利用固定局部带宽；MATLAB 稀疏 LU 的具体常数和填充由实现决定。

Green 边界计算约为 `O(N_θN_s)`，其中 `N_s` 是压缩后源点数，另有 `O(N)` 的源扫描。工作数组与缓存按 `O(N)` 增长，8点积分引入常数因子。实现没有全二维稠密矩阵或径向稠密插值逆。

固定网格且首阶段缓存命中时，整个时间步包含3次新 RHS，另有规范探测、诊断和可能的换网格。规范探测会重复调用当前 WENO 组装，候选网格会重建算子和 Poisson 分解；失败候选回到旧几何也可能触发缓存重建。因此固定矩阵缓存消除了重复分解，却没有消除所有状态相关工作。

具体计时、历史性能对照和运行分段记录见[验证记录](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_VALIDATION_NOTES.md)。
