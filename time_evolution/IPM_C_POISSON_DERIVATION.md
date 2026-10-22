# 极坐标 IPM：二维到径向一维的完整推导

这里完整推导当前 Poisson 模块的降维过程：二维 Cartesian 方程 → 极坐标方程 → 已知墙源修正 → 角向正弦投影 → 径向边值方程。随后保留径向样条与 Galerkin 离散的实现公式。流程说明见[Poisson 分步讲解](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_ALGORITHM_DETAILS.md#poisson)。

<a id="angular-derivation"></a>

## 一、从二维方程推到径向边值方程

### 1. 从固定时间阶段的二维问题开始

固定一个时间阶段，省略时间变量。已知当前密度 \(R(X,Y)\)，定义源

\[
q(X,Y)=\partial_XR(X,Y).
\]

未知为流函数 \(\Psi(X,Y)\)，满足

\[
\boxed{-\partial_{XX}\Psi-\partial_{YY}\Psi=q.}
\]

计算域为上半圆，半径记为 \(R_b\)。直墙上流函数为零；外弧值记为已知函数 \(g_b\)。它在程序中由源的镜像 Green 积分近似计算，后面的降维只需要将其作为给定边界数据。

### 2. 用链式法则把二维算子写成极坐标形式

令

\[
X=r\cos\theta,\qquad Y=r\sin\theta,\qquad
0<r<R_b,\quad 0<\theta<\pi.
\]

逆变换给出

\[
r_X=\cos\theta,\quad r_Y=\sin\theta,\qquad
\theta_X=-\frac{\sin\theta}{r},\quad
\theta_Y=\frac{\cos\theta}{r}.
\]

所以，对任意光滑函数，

\[
\partial_X=\cos\theta\,\partial_r-\frac{\sin\theta}{r}\partial_\theta,
\qquad
\partial_Y=\sin\theta\,\partial_r+\frac{\cos\theta}{r}\partial_\theta.
\]

为展示二阶微分中的系数导数，暂记 \(c=\cos\theta,s=\sin\theta\)。再次使用链式法则，得到

\[
\begin{aligned}
\Psi_{XX}
={}&c^2\Psi_{rr}
-\frac{2sc}{r}\Psi_{r\theta}
+\frac{s^2}{r^2}\Psi_{\theta\theta}
+\frac{s^2}{r}\Psi_r
+\frac{2sc}{r^2}\Psi_\theta,\\
\Psi_{YY}
={}&s^2\Psi_{rr}
+\frac{2sc}{r}\Psi_{r\theta}
+\frac{c^2}{r^2}\Psi_{\theta\theta}
+\frac{c^2}{r}\Psi_r
-\frac{2sc}{r^2}\Psi_\theta.
\end{aligned}
\]

相加时，混合导数项和单独的 \(\Psi_\theta\) 项抵消；使用 \(c^2+s^2=1\)，得到

\[
\Delta_{X,Y}\Psi
=\Psi_{rr}+\frac1r\Psi_r+\frac1{r^2}\Psi_{\theta\theta}.
\]

同一个变换也作用于源：

\[
q(r,\theta)
=\cos\theta\,R_r-\frac{\sin\theta}{r}R_\theta.
\]

于是完整二维边值问题变成

\[
\boxed{
-\left(\Psi_{rr}+\frac1r\Psi_r+
\frac1{r^2}\Psi_{\theta\theta}\right)=q(r,\theta),
}
\]
\[
\Psi(r,0)=\Psi(r,\pi)=0,\qquad
\Psi(R_b,\theta)=g_b(\theta).
\]

以上含 \(1/r\) 的强方程在 \(r>0\) 推导。原点端条件稍后由单值性给出。

### 3. 构造已知墙源修正，并推导余项方程

记两条墙上的源为

\[
q_+(r)=q(r,0),\qquad q_-(r)=q(r,\pi).
\]

墙上 \(\Psi=0\) 不要求 \(q_\pm=0\)。但是，任何有限正弦和在墙上的值、径向导数和角向二阶导数均为零，其墙上 Laplacian 也为零。当前实现用已知函数 \(L\) 承担非零墙源对应的曲率，再用正弦表示余项。

选取

\[
\chi_+(\theta)=\frac14\sin^2\theta(1+\cos\theta),\qquad
\chi_-(\theta)=\frac14\sin^2\theta(1-\cos\theta).
\]

它们满足

\[
\chi_\pm(0)=\chi_\pm(\pi)=0,
\qquad \chi_\pm'(0)=\chi_\pm'(\pi)=0,
\]
\[
\chi_+''(0)=1,\quad\chi_+''(\pi)=0,\qquad
\chi_-''(0)=0,\quad\chi_-''(\pi)=1.
\]

定义

\[
\ell_\pm(r)=-r^2q_\pm(r),\qquad
L(r,\theta)=\ell_+(r)\chi_+(\theta)+\ell_-(r)\chi_-(\theta).
\]

在墙上，\(L=L_r=L_{rr}=0\)。因此

\[
\left.\Delta L\right|_{\theta=0}
=\frac{\ell_+(r)}{r^2}=-q_+(r),\qquad
\left.\Delta L\right|_{\theta=\pi}
=\frac{\ell_-(r)}{r^2}=-q_-(r).
\]

所以 \(L\) 保留零墙值，并满足墙上的 \(-\Delta L=q\)。写成

\[
\Psi=L+\Phi,
\]

代回原方程并移项：

\[
-\Delta(L+\Phi)=q
\quad\Longrightarrow\quad
\boxed{-\Delta\Phi=q+\Delta L=:f(r,\theta).}
\]

边界条件同步变为

\[
\Phi(r,0)=\Phi(r,\pi)=0,\qquad
\Phi(R_b,\theta)=g_b(\theta)-L(R_b,\theta).
\]

由刚才的墙上恒等式，修正源在墙上满足 \(f(r,0)=f(r,\pi)=0\)。这项修正服务于当前有限正弦表示及墙端强方程；一般 Dirichlet Poisson 问题也可以采用其他离散方法。

<a id="wall-compatibility"></a>

#### 二阶修正没有消除四阶及更高阶墙面相容性

上述修正只匹配第一层非零偶阶墙导数，并不意味着高阶墙源信息在二阶终止。下面在 \(r>0\) 的光滑墙段上讨论，假设涉及的导数存在；原点、外弧角点以及径向分片结点的正则性须另行处理。

记径向算子

\[
\mathcal A_r=\partial_{rr}+\frac1r\partial_r.
\]

余项方程为

\[
-\mathcal A_r\Phi-\frac1{r^2}\Phi_{\theta\theta}=f.
\]

因为墙值 \(\Phi|_{\rm wall}=0\) 对整段半径成立，其径向导数也为零，故

\[
\Phi_{\theta\theta}|_{\rm wall}=-r^2 f|_{\rm wall}=0.
\]

再对方程作两次角向微分：

\[
-\mathcal A_r\Phi_{\theta\theta}
-\frac1{r^2}\Phi_{\theta\theta\theta\theta}=f_{\theta\theta}.
\]

二阶修正已使 \(\Phi_{\theta\theta}|_{\rm wall}=0\) 沿径向恒成立，因此

\[
\boxed{\Phi_{\theta\theta\theta\theta}|_{\rm wall}
=-r^2 f_{\theta\theta}|_{\rm wall}.}
\]

右边一般非零。要让四阶墙导数也归零，还需要额外匹配 \(f_{\theta\theta}|_{\rm wall}=0\)。更一般地，令
\(T_{2j}(r)=\partial_\theta^{2j}\Phi|_{\rm wall}\)，有递推

\[
T_{2j+2}
=-r^2\left(
\mathcal A_r T_{2j}
+\partial_\theta^{2j}f|_{\rm wall}
\right),\qquad T_0=0.
\]

这些高阶关系由 PDE 和已知源决定，不是要额外指定的独立边界数据。

例如，取满足零直墙值的光滑函数

\[
\Phi(X,Y)=Y^4=r^4\sin^4\theta,\qquad
f=-\Delta\Phi=-12Y^2.
\]

它满足 \(f|_{\rm wall}=0\) 和 \(\Phi_{\theta\theta}|_{\rm wall}=0\)，但

\[
\Phi_{\theta\theta\theta\theta}(r,0)=24r^4\ne0.
\]

其墙源已经为零，当前二阶修正 \(L\) 也为零；高阶不相容仍然存在。

有限正弦和的所有偶阶角导数在墙上均为零，这是基空间的性质。真解可以由无限正弦级数表示，但不能因此将任意高阶导数都逐项求导后取墙值：相应导数级数未必在墙附近一致收敛。前面的径向方程推导只对角方向分部积分两次，并不要求四阶或所有高阶墙导数为零，所以降维依然成立。

在固定半径、角向足够光滑的条件下，墙上二阶导非零通常产生正弦系数的 \(m^{-3}\) 主尾；消去它后，非零四阶墙导一般留下 \(m^{-5}\) 主尾。二阶修正改善了角向衰减，却不保证指数谱收敛。是否足够解析当前流场和速度，仍需角向加密检验。

当前代码只实施上述二阶墙源修正，没有匹配四阶墙导的额外修正。若需要提高墙端光滑性，可继续构造高阶修正，使新的余项源满足 \(f_{\theta\theta}|_{\rm wall}=0\) 等条件；所需源导数的可用性和稳定性也需要验证。

<a id="angular-accuracy"></a>

#### 二阶墙源修正能给出什么条件性精度界

这里先只分析精确解的角向正弦投影截断。固定 \(r>0\)，记 \(v(\theta)=\Phi(r,\theta)\)，假设 \(v\in W^{5,1}(0,\pi)\)，且两端满足 \(v=v''=0\)。其正弦系数为

\[
b_m=\frac2\pi\int_0^\pi v(\theta)\sin(m\theta)\,d\theta.
\]

连续分部积分给出

\[
b_m=\frac{2}{\pi m^5}
\left[
v^{(4)}(0)-(-1)^m v^{(4)}(\pi)
+\int_0^\pi v^{(5)}(\theta)\cos(m\theta)\,d\theta
\right].
\]

因此存在与 \(m\) 无关的常数 \(C(r)\)，使 \(|b_m|\le C(r)m^{-5}\)。保留前 \(M\) 个角模，令
\(S_Mv=\sum_{m=1}^M b_m\sin(m\theta)\)。直接估计级数尾部：

\[
\begin{aligned}
\|v-S_Mv\|_\infty
&\le C(r)\sum_{m>M}m^{-5}=O(M^{-4}),\\
\|\partial_\theta(v-S_Mv)\|_\infty
&\le C(r)\sum_{m>M}m^{-4}=O(M^{-3}),\\
\|\partial_{\theta\theta}(v-S_Mv)\|_\infty
&\le C(r)\sum_{m>M}m^{-3}=O(M^{-2}).
\end{aligned}
\]

因为角间距与 \(M^{-1}\) 同阶，这分别是流函数四阶、一阶角导三阶、二阶角导二阶的角向最大范数误差界。它们是全角区间上的保守估计；特殊函数、特定采样点或远离墙的局部区域可能表现更高阶。

|修正后余项在两墙匹配的条件|一般正弦系数尾|流函数截断界|一阶角导截断界|二阶角导截断界|
|---|---|---|---|---|
|\(v=v''=0\)|\(O(m^{-5})\)|\(O(M^{-4})\)|\(O(M^{-3})\)|\(O(M^{-2})\)|
|\(v=v''=v^{(4)}=0\)，且更光滑|\(O(m^{-7})\)|\(O(M^{-6})\)|\(O(M^{-5})\)|\(O(M^{-4})\)|

速度分量 \(U_r=-\Psi_\theta/r\) 涉及一次角微分，因此在 \(r\ge r_*>0\)、上述常数对半径一致有界的区域，其角截断部分可得到三阶界。另一分量 \(U_\theta=\Psi_r\) 还需要单独控制径向导数的模态系数；仅由 \(|b_m|\le Cm^{-5}\) 不能自动推出这个控制。若相应径向导数也有统一角向正则性，才可进一步给出整个速度的相应区域估计。以上固定半径结论不能直接延伸为含原点的全域保证。

这些估计针对精确连续投影的截断，不是当前完整 DST／径向 Galerkin 实现的误差定理。程序还包含源微分与重构、有限角节点变换、径向离散与求积、Green 外弧近似、原点闭合等误差；各项的一致稳定性和随加密的误差都需要控制。

径向三次 B 样条 Galerkin 在光滑解、适当规则网格与一致稳定离散等条件下，标准目标是解的 \(L^2\) 四阶和能量／\(H^1\) 三阶；它们不是最大范数结论，也不能直接作为当前拉伸网格与原点处理的精度证明。现有特定制造解的节点测量见[验证记录](/Users/cbjin/Documents/one_sided_IPM/time_evolution/IPM_C_VALIDATION_NOTES.md)，当前代码尚无一般全域四阶流函数或四阶速度的保证。

### 4. 定义角向模态，并逐项投影二维方程

为了区分连续投影与程序的离散变换，连续投影记为

\[
\mathcal S_m[v](r)
=\frac2\pi\int_0^\pi v(r,\theta)\sin(m\theta)\,d\theta,
\qquad m=1,2,\ldots.
\]

正弦函数具有正交性

\[
\int_0^\pi\sin(k\theta)\sin(m\theta)\,d\theta
=\frac\pi2\delta_{km}.
\]

定义未知径向模态与已知模态源

\[
\phi_m(r)=\mathcal S_m[\Phi](r),\qquad
f_m(r)=\mathcal S_m[f](r).
\]

在适当收敛意义下，

\[
\Phi(r,\theta)=\sum_{m=1}^\infty\phi_m(r)\sin(m\theta).
\]

对余项方程

\[
-\Phi_{rr}-\frac1r\Phi_r-\frac1{r^2}\Phi_{\theta\theta}=f
\]

两边乘 \((2/\pi)\sin(m\theta)\)，沿角方向积分。因为积分区间不依赖 \(r\)，径向项给出

\[
\mathcal S_m[\Phi_{rr}]=\phi_m''(r),\qquad
\mathcal S_m[\Phi_r]=\phi_m'(r).
\]

角向二阶项要两次分部积分。完整计算为

\[
\begin{aligned}
\int_0^\pi\Phi_{\theta\theta}\sin(m\theta)\,d\theta
={}&[\Phi_\theta\sin(m\theta)]_0^\pi
-m\int_0^\pi\Phi_\theta\cos(m\theta)\,d\theta\\
={}&[\Phi_\theta\sin(m\theta)]_0^\pi
-m[\Phi\cos(m\theta)]_0^\pi\\
&-m^2\int_0^\pi\Phi\sin(m\theta)\,d\theta.
\end{aligned}
\]

第一项因为正弦在两端为零而消失，第二项因为 \(\Phi\) 在两条墙上为零而消失。因此

\[
\boxed{\mathcal S_m[\Phi_{\theta\theta}]=-m^2\phi_m.}
\]

将三个投影结果代回，注意原方程角向项前还有负号：

\[
-\phi_m''-\frac1r\phi_m'
-\frac1{r^2}(-m^2\phi_m)=f_m.
\]

整理得到每个角模的径向方程

\[
\boxed{
-\phi_m''(r)-\frac1r\phi_m'(r)
+\frac{m^2}{r^2}\phi_m(r)=f_m(r).
}
\]

也可乘以 \(r\)，写成后续径向弱式使用的形式

\[
\boxed{-(r\phi_m')'+\frac{m^2}{r}\phi_m=r f_m.}
\]

二维问题的角方向至此已经消去。每个 \(m\) 对应一个未知函数 \(\phi_m(r)\)，不同模态互不耦合；径向变量尚未离散。

### 5. 把模态右端及两端条件完整写出

定义

\[
q_m=\mathcal S_m[q],\qquad
\chi_{\sigma,m}=\mathcal S_m[\chi_\sigma],\qquad
\sigma\in\{+,-\}.
\]

已知修正的 Laplacian 为

\[
\Delta L=
\sum_{\sigma\in\{+,-\}}
\left[
\left(\ell_\sigma''+\frac{\ell_\sigma'}r\right)\chi_\sigma
+\frac{\ell_\sigma}{r^2}\chi_\sigma''
\right].
\]

因此，径向方程的已知右端为

\[
\boxed{
f_m=q_m+
\sum_{\sigma\in\{+,-\}}
\left[
\left(\ell_\sigma''+\frac{\ell_\sigma'}r\right)\chi_{\sigma,m}
+\frac{\ell_\sigma}{r^2}\mathcal S_m[\chi_\sigma'']
\right].
}
\]

连续投影下，\(\chi_\sigma\) 的墙值为零，同样可以两次分部积分，得到
\(\mathcal S_m[\chi_\sigma'']=-m^2\chi_{\sigma,m}\)。于是也可写成

\[
f_m=q_m+
\sum_{\sigma\in\{+,-\}}
\left(
\ell_\sigma''+\frac{\ell_\sigma'}r
-\frac{m^2}{r^2}\ell_\sigma
\right)\chi_{\sigma,m}.
\]

外端条件是对修正后的外弧值作同一个投影：

\[
\phi_m(R_b)=\gamma_m,\qquad
\gamma_m=\frac2\pi\int_0^\pi
[g_b(\theta)-L(R_b,\theta)]\sin(m\theta)\,d\theta.
\]

原点是一个物理点，其墙值为零。对正则解，\(\Psi(0,\theta)=0\)；源在原点有界时 \(L=O(r^2)\)，从而 \(\Phi(0,\theta)=0\)，给出

\[
\phi_m(0)=0.
\]

完整的一维边值问题最终为

\[
\boxed{
\begin{cases}
-\phi_m''-\dfrac1r\phi_m'
+\dfrac{m^2}{r^2}\phi_m=f_m,&0<r<R_b,\\[2mm]
\phi_m(0)=0,\qquad \phi_m(R_b)=\gamma_m,
&m=1,2,\ldots.
\end{cases}}
\]

原点端值是正则性的必要条件。更强的 Cartesian 光滑性还约束高模态的一阶斜率，不能单凭端值条件就断言已满足；也不能在这里普遍要求所有模态均为 \(O(r^m)\)。

### 6. 当前程序如何将连续投影换成离散正弦变换

上面的推导使用连续角积分。程序只保留有限个角模，并在角节点
\(\theta_j=j\pi/(N_\theta-1)\) 上使用 DST-I：

\[
\widehat v_m=
\frac2{N_\theta-1}
\sum_{j=1}^{N_\theta-2}v(\theta_j)\sin(m\theta_j),
\qquad m=1,\ldots,N_\theta-2.
\]

帽号在本文件其余实现公式中专指这个离散变换。对有限正弦展开 \(\Phi_N\)，其解析角向二阶导数依然逐模乘 \(-m^2\)；离散正交性让各个未知模态保持解耦。

已知源则先沿角向作 DST，再用局部径向 Hermite 插值，得到 \(\widetilde q_m(r)\)；两条墙源单独重构。以这些重构墙源构造 \(\ell_\pm=-r^2\widetilde q_\pm\) 后，代码对应的模态右端是

\[
\boxed{
f_m^{(d)}=\widetilde q_m+
\sum_{\sigma\in\{+,-\}}
\left[
\left(\ell_\sigma''+\frac{\ell_\sigma'}r\right)
\widehat{\chi_\sigma}_m
+\frac{\ell_\sigma}{r^2}\widehat{\chi_\sigma''}_m
\right],
\qquad
\gamma_m^{(d)}
=\widehat{\,g_b-L(R_b,\cdot)\,}_m.
}
\]

径向方程的左端形式保持不变，右端和外端值换成上述离散变换结果。这里 \(\chi_\sigma\) 与解析 \(\chi_\sigma''\) 分别采样、变换，有限角网格上一般有

\[
\widehat{\chi_\sigma''}_m\ne
-m^2\widehat{\chi_\sigma}_m.
\]

因此不能将连续投影中的简化恒等式直接替换进当前代码的已知修正项。当前径向源为 \(C^1\) 分片 Hermite 重构，修正的径向强式按区间内部理解；后续全局求解使用弱形式。

**到这里，二维到径向一维的推导结束。** 下一步才是用 B 样条表示每个 \(\phi_m(r)\)，并用 Galerkin 将一维方程变成系数线性系统。

来源：[Poisson 的源、修正与负载组装](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+spatial/poisson.m:35)。

## 二、源重构、径向样条与 Galerkin 离散的实现公式


已知节点源 \(\mathbf q=\mathsf D_X\mathbf W\)，即完整密度的水平导数，本模块求满足墙面零值和外弧值 \(g_b(\theta)\) 的流函数。下面依次把同一个未知量从二维函数变成径向模态，再变成样条系数，最后还原为流场。

### 1. 将节点源变成径向模态源

记保留的角模式数为 \(N_m=N_\theta-2\)。对角向节点数组 \((v_j)\)，统一定义 DST-I：

\[
\widehat v_m=\frac{2}{N_\theta-1}
\sum_{j=1}^{N_\theta-2}v_j\sin(m\theta_j),
\qquad m=1,\ldots,N_m.
\]

在内部角节点，其逆变换为 \(v_j=\sum_{m=1}^{N_m}\widehat v_m\sin(m\theta_j)\)。帽号始终表示这个**离散变换**，不表示连续 Fourier 积分。变换由长度 \(2(N_\theta-1)\) 的奇延拓 FFT 实现；128个角节点对应126个模式、长度254的 FFT。

对每个径向节点变换源，得到 \(\widehat q_{im}\)。再用物理径向坐标上的局部 \(C^1\) 三次 Hermite 插值算子 \(\mathcal I_r\) 定义

\[
\widetilde q_m(r)=\mathcal I_r[(\widehat q_{im})_i](r),\qquad
\widetilde q_+(r)=\mathcal I_r[(\mathbf q_{i0})_i](r),\qquad
\widetilde q_-(r)=\mathcal I_r[(\mathbf q_{i,N_\theta-1})_i](r).
\]

\(\widetilde q_m\) 是第 \(m\) 个内部角模的径向源函数；\(\widetilde q_+,\widetilde q_-\) 分别是 \(\theta=0,\pi\) 的墙源。它们是源的插值，尚不是流函数的有限元解。

\(\mathcal I_r\) 的定义对所有源列相同。给定任一节点列 \(\mathbf v=(v_i)\)，先用缓存的五点物理径向差分矩阵 \(\mathsf D_r^{(5)}\) 计算节点斜率。在区间 \([r_i,r_{i+1}]\) 内令 \(\eta=(r-r_i)/\Delta r_i\)，采用标准 Hermite 基

\[
\begin{aligned}
H_{00}&=2\eta^3-3\eta^2+1,&H_{01}&=-2\eta^3+3\eta^2,\\
H_{10}&=\eta^3-2\eta^2+\eta,&H_{11}&=\eta^3-\eta^2,
\end{aligned}
\]
\[
\mathcal I_r[\mathbf v](r)=H_{00}v_i+H_{01}v_{i+1}
+\Delta r_i\left[H_{10}(\mathsf D_r^{(5)}\mathbf v)_i
+H_{11}(\mathsf D_r^{(5)}\mathbf v)_{i+1}\right].
\]

相邻区间共享节点值和斜率，故重构为 \(C^1\)。单区间最多依赖6个原始节点，没有全局曲率逆；它在大尺度径向尾部避免了旧全局源插值的非局部传播。此局部性不等于无超调保证。

### 2. 提升墙源，定义要展开的流函数

墙上 \(\Psi=0\) 不意味着 \(q=0\)。为使墙源参与离散，定义提升函数

\[
\chi_+(\theta)=\tfrac14\sin^2\theta(1+\cos\theta),\qquad
\chi_-(\theta)=\tfrac14\sin^2\theta(1-\cos\theta),
\]
\[
\ell_\pm(r)=-r^2\widetilde q_\pm(r),\qquad
L(r,\theta)=\ell_+(r)\chi_+(\theta)+\ell_-(r)\chi_-(\theta).
\]

两端均有 \(\chi_\pm=\chi_\pm'=0\)，并且
\(\chi_+''(0)=1,\chi_+''(\pi)=0,\chi_-''(0)=0,\chi_-''(\pi)=1\)。因此 \(L\) 不改变零墙值，而 \(-\Delta L\) 在墙上承接相应源值。

将剩余流函数记为 \(\Phi=\Psi-L\)。连续分解为 \(-\Delta\Phi=q+\Delta L\)；离散时使用上面已定义的源重构。提升的径向导数直接由同一重构计算：

\[
\ell_\pm'=-2r\widetilde q_\pm-r^2\widetilde q_\pm'.
\]

没有对乘完 \(r^2\) 的节点值另做一次三次插值。

### 3. 角向分解得到径向微分方程

**展开的是剩余流函数 \(\Phi\)**，其第 \(m\) 个径向模态函数统一记为 \(\phi_m\)：

\[
\Phi(r,\theta)\approx\sum_{m=1}^{N_m}\phi_m(r)\sin(m\theta),\qquad
\Psi\approx L+\sum_{m=1}^{N_m}\phi_m\sin(m\theta).
\]

正弦基使 \(\partial_{\theta\theta}\sin(m\theta)=-m^2\sin(m\theta)\)。将同一个 DST 用于源及提升项后，每个模态满足

\[
\boxed{
-\phi_m''-\frac1r\phi_m'+\frac{m^2}{r^2}\phi_m=f_m(r),
\quad \phi_m(0)=0,\quad \phi_m(R_b)=\gamma_m.
}
\]

这里，径向强方程右端和外弧模态值分别是

\[
f_m=\widetilde q_m+
\sum_{\sigma\in\{+,-\}}
\left[
\left(\ell_\sigma''+\frac{\ell_\sigma'}r\right)\widehat{\chi_\sigma}_m
+\frac{\ell_\sigma}{r^2}\widehat{\chi_\sigma''}_m
\right],
\qquad
\gamma_m=\widehat{\,g_b-L(R_b,\cdot)\,}_m.
\]

\(\widehat{\chi_\sigma}_m\) 和 \(\widehat{\chi_\sigma''}_m\) 分别变换。有限角网格下，后者一般不等于 \(-m^2\widehat{\chi_\sigma}_m\)。到这里，二维 Poisson 已成为 \(N_m\) 个互不耦合的径向边值问题；径向未知函数尚未作有限元离散。这里的角向方法是配点变换，下一步只对径向作 Galerkin 离散。

### 4. 对同一个径向模态作 B 样条展开

在前面零起点编号的径向节点上，开放夹持结点向量为

\[
T=(0,0,0,0,r_1,\ldots,r_{N_r-2},R_b,R_b,R_b,R_b).
\]

由它构造四阶 \(k=4\)、三次 \(p=3\) 的 B 样条基 \(N_a(r)\)，共 \(N_b=N_r+2\) 个。余项 \(\Phi_h\) 的径向展开在内部单结点处为 \(C^2\)；提升 \(L\) 使用前述 \(C^1\) 源重构，故完整 \(\Psi_h=L+\Phi_h\) 在正半径内部结点一般仅为 \(C^1\)，不能称为全局 \(C^2\) 流函数。Cox–de Boor 递推及解析导数在本地实现，不依赖样条工具箱。

\(\phi_m\) 的有限元近似叫作 \(\phi_{m,h}\)，而不是另一个未知场：

\[
\phi_{m,h}(r)=\sum_{a=1}^{N_b}z_{am}N_a(r).
\]

\(z_{am}\) 是第 \(m\) 个模态、第 \(a\) 个基函数的**展开系数**，不是 \(\phi_m\) 的节点值。夹持基满足 \(N_1(0)=1\)、\(N_{N_b}(R_b)=1\)，故边界条件直接给出

\[
z_{1m}=0,\qquad z_{N_b,m}=\gamma_m.
\]

令内部未知的基函数索引为 \(\mathcal I=\{2,\ldots,N_b-1\}\)，则

\[
\boxed{\phi_{m,h}(r)=\sum_{a\in\mathcal I}z_{am}N_a(r)+\gamma_mN_{N_b}(r).}
\]

从角向分解到径向有限元，未知量始终是同一个 \(\phi_m\)；有限元仅把求函数变成求系数 \((z_{am})_{a\in\mathcal I}\)。

### 5. 从径向方程得到弱形式与矩阵

径向方程乘 \(rN_a\)、在 \([0,R_b]\) 积分并分部积分。对 \(a\in\mathcal I\)，测试函数两端为零，得到

\[
\int_0^{R_b}rN_a'\phi_{m,h}'\,dr
+m^2\int_0^{R_b}\frac{N_a\phi_{m,h}}r\,dr
=\int_0^{R_b}rN_af_m\,dr.
\]

将 \(f_m\) 中提升的径向二阶项也分部积分，右端记为弱负载 \(d_{am}\)：

\[
\begin{aligned}
d_{am}={}&\int_0^{R_b}rN_a\widetilde q_m\,dr
-\sum_{\sigma\in\{+,-\}}\widehat{\chi_\sigma}_m
\int_0^{R_b}rN_a'\ell_\sigma'\,dr\\
&+\sum_{\sigma\in\{+,-\}}\widehat{\chi_\sigma''}_m
\int_0^{R_b}N_a\frac{\ell_\sigma}{r}\,dr.
\end{aligned}
\]

**先去掉原点基函数 \(N_1\)**，再在 \(\mathcal J=\{2,\ldots,N_b\}\) 上组装矩阵 \(\mathbf K,\mathbf M\)，其元素为

\[
K_{ab}=\int_0^{R_b}rN_a'N_b'\,dr,\qquad
M_{ab}=\int_0^{R_b}\frac{N_aN_b}{r}\,dr,
\quad a,b\in\mathcal J.
\]

保留基函数在原点均为零，使上述积分可积。不能先组装包含 \(N_1\) 的完整 \(\mathbf M\) 再删行列，因为 \(M_{11}\) 发散。每个径向区间用8点 Gauss–Legendre 积分，积分点严格满足 \(r>0\)；含 \(1/r\) 的积分通常不具有多项式求积的精确性。

把已知的外弧系数移到右端，定义

\[
\begin{aligned}
\mathbf z_m&=(z_{am})_{a\in\mathcal I},\qquad
\mathbf d_m=(d_{am})_{a\in\mathcal I},\\
\mathbf A_m&=\mathbf K_{\mathcal I\mathcal I}+m^2\mathbf M_{\mathcal I\mathcal I},\\
\mathbf b_m&=\mathbf d_m-
\left(\mathbf K_{\mathcal I,N_b}+m^2\mathbf M_{\mathcal I,N_b}\right)\gamma_m.
\end{aligned}
\]

最终且唯一的模态线性系统是

\[
\boxed{\mathbf A_m\mathbf z_m=\mathbf b_m,\qquad m=1,\ldots,N_m.}
\]

因此 \(f_m\) 是强方程的函数右端，\(\mathbf d_m\) 是弱积分负载，\(\mathbf b_m\) 是扣除边界列后的线性系统右端，三者不再混用。385×128网格下，\(N_b=387\)，每模有385个内部未知量，共求126个385阶带状系统。

### 6. 缓存分解，直接求系数

令 \(\mathbf D_m=\operatorname{diag}((\mathbf A_m)_{aa}^{-1/2})\)。固定网格时缓存 \(\operatorname{LU}(\mathbf D_m\mathbf A_m\mathbf D_m)\)，每次求

\[
(\mathbf D_m\mathbf A_m\mathbf D_m)(\mathbf D_m^{-1}\mathbf z_m)
=\mathbf D_m\mathbf b_m,
\]

再恢复同一个 \(\mathbf z_m\)。实现执行缓存分解的回代，不显式求逆，也不作 GMRES 迭代。

第1步的源重构复用固定五点斜率映射和局部 Hermite 基映射，不再求全局曲率系统。所有 \(N_m\) 个模态源和两条墙源一起作为数据列，计算积分点、节点和审计点的值及解析导数；这一步不新增线性分解。

缓存按实际 \(r,\theta\) 精确匹配，只保留一个几何版本。B 样条基、积分映射、模态矩阵及 LU、源 Hermite 映射和 Green 积分几何随网格固定；源、边界值、负载、解和诊断随状态刷新。换网格或重启 MATLAB 后重建。

### 7. 重构流函数并求速度

把求出的 \(\mathbf z_m\) 放回第4步，再放回第3步，闭合整个求解链：

\[
\boxed{
\Psi_h(r,\theta)=L(r,\theta)+\sum_{m=1}^{N_m}
\left[\sum_{a\in\mathcal I}z_{am}N_a(r)+\gamma_mN_{N_b}(r)\right]\sin(m\theta).
}
\]

径向导数只把同一式中的 \(N_a\) 换成 \(N_a'\)，角向导数只把 \(\sin(m\theta)\) 换成 \(m\cos(m\theta)\)，并分别加上 \(L_r,L_\theta\)。随后在正径节点计算

\[
U_r=-\frac{\partial_\theta\Psi_h}{r},\quad U_\theta=\partial_r\Psi_h,
\quad U_1=U_r\cos\theta-U_\theta\sin\theta,
\quad U_2=U_r\sin\theta+U_\theta\cos\theta.
\]

原点不计算 \(1/r\)。记 \(\sigma_m=\phi_{m,h}'(0)\)，代码取唯一 Cartesian 速度 \(U_1(0)=-\sigma_1,U_2(0)=0\)。提升为 \(O(r^2)\)，不贡献一阶斜率。

当前空间只显式约束 \(\phi_{m,h}(0)=0\)，高模态的 \(\sigma_m\) 可能非零；实现记录其大小及相对 \(\sigma_1\) 的比值，尚未完全约束这些一阶缺陷。余项的径向 \(C^2\) 和完整流函数的径向 \(C^1\) 都不自动保证 Cartesian 原点具有单值一阶导数。一般半平面非零墙源下，也不能不加条件地要求所有模态满足 \(O(r^m)\)；高模态的一阶斜率缺陷仍需随加密检验。

### 外弧值 \(g_b\) 的来源

外弧值在上述边界消元之前计算。半平面镜像核记为

\[
\mathcal G((X,Y),(X_s,Y_s))=
\frac1{4\pi}\log\left(1+\frac{4YY_s}{(X-X_s)^2+(Y-Y_s)^2}\right).
\]

令 \(w_i^r,w_j^\theta\) 分别为实际径轴和均匀角轴的一维梯形权重。对原生内部节点，极坐标面积权重为 \(r_iw_i^rw_j^\theta\)，因此

\[
g_b(\theta)\approx\sum_{i,j}\mathcal G((R_b\cos\theta,R_b\sin\theta),(X_{ij},Y_{ij}))
\,\mathbf q_{ij}r_iw_i^rw_j^\theta.
\]

当前配置使用相对源阈值 \(10^{-12}\) 和最多48000个压缩源点。超预算时按径向/角向索引分箱，正负源分别聚合强度和质心；外弧目标按64个一批计算。底层函数在不传参数时有不同默认值。

阈值仅用于 Green 外边界压缩，内部负载使用完整源。有限外域、梯形面积求积和源压缩分别引入边界近似。

来源：[Poisson 完整实现](/Users/cbjin/Documents/one_sided_IPM/time_evolution/+ipm_c/+spatial/poisson.m:1)。

