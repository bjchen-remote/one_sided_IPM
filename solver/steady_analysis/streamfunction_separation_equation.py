"""Exact rank-one separation equations and a smooth-origin obstruction.

No fitted angular basis, radial power ansatz, numerical data, or PDE evolution.
"""
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_separation_20260929'


def main():
    t,th,r=sp.symbols('t theta r',real=True)
    lam,mu,v,sigma=sp.symbols('lambda mu v sigma',real=True)
    a=sp.Function('a')(t);F=sp.Function('F')(th);H=sp.Function('H')(t,th)
    Dx=lambda f:sp.exp(-t)*(sp.cos(th)*sp.diff(f,t)-sp.sin(th)*sp.diff(f,th))
    Dy=lambda f:sp.exp(-t)*(sp.sin(th)*sp.diff(f,t)+sp.cos(th)*sp.diff(f,th))
    psi=sp.exp(2*t)*a*F
    ux=-Dy(psi);uy=Dx(psi)
    q=-(sp.diff(a,t,2)+4*sp.diff(a,t)+4*a)*F-a*sp.diff(F,th,2)
    poisson=sp.trigsimp(sp.expand(-Dx(Dx(psi))-Dy(Dy(psi))-q))
    assert poisson==0
    rho=sp.exp(t)*H
    L=sp.cos(th)*(H+sp.diff(H,t))-sp.sin(th)*sp.diff(H,th)
    assert sp.trigsimp(sp.expand(Dx(rho)-L))==0
    direct=((ux+lam*sp.exp(t)*sp.cos(th)+v)*Dx(rho)
            +(uy+lam*sp.exp(t)*sp.sin(th))*Dy(rho)-mu*rho-v)*sp.exp(-t)
    reduced=(lam-a*sp.diff(F,th))*(H+sp.diff(H,t))+(2*a+sp.diff(a,t))*F*sp.diff(H,th)-mu*H+v*sp.exp(-t)*(L-1)
    transport=sp.trigsimp(sp.expand(direct-reduced))
    assert transport==0

    # Verify the characteristic integral inverse algebraically. I and It
    # denote its integral and t derivative before differentiating the lower
    # limit; this avoids asking CAS to evaluate an unspecified source.
    I,It,Q=sp.symbols('I I_t Q')
    K=sp.sin(th)*I;Kt=sp.sin(th)*It
    Kth=sp.cos(th)*I-Q/sp.sin(th)+sp.cos(th)*It
    inverse_certificate=sp.trigsimp(sp.cos(th)*(K+Kt)-sp.sin(th)*Kth-Q)
    assert inverse_certificate==0
    Cfun=sp.Function('C')
    Cpart=sp.exp(-t)*Cfun(sp.exp(t)*sp.sin(th))
    Ccertificate=sp.simplify(sp.cos(th)*(Cpart+sp.diff(Cpart,t))-sp.sin(th)*sp.diff(Cpart,th))
    assert Ccertificate==0

    # Smooth-origin angular rigidity, followed by exact wall parity identity.
    f2=sigma*sp.sin(th)*sp.cos(th)-sp.sin(th)**2/2
    assert sp.trigsimp(sp.diff(f2,th,2)+4*f2+1)==0
    h=sp.Function('h')(r)
    Ar=r*r*h
    qr=-((sp.diff(Ar,r,2)+sp.diff(Ar,r)/r)*f2+Ar/r**2*sp.diff(f2,th,2))
    qexpected=h-(r*r*sp.diff(h,r,2)+5*r*sp.diff(h,r))*f2
    assert sp.trigsimp(sp.expand(qr-qexpected))==0
    hv,S=sp.symbols('h S',real=True)
    Eplus=(r*(lam-sigma*hv)+v)*hv-mu*S-v
    Eminus=(-r*(lam-sigma*hv)+v)*hv+mu*S-v
    parity=sp.factor(Eplus+Eminus)
    assert sp.expand(parity-2*v*(hv-1))==0
    assert sp.simplify((Eplus-Eminus).subs({hv:1,S:r})-2*r*(lam-mu-sigma))==0
    # Rank-two angular-space audit: parity belongs to finite homogeneous
    # Taylor polynomials, not to an assumed analytic expansion of psi.
    k3=sp.symbols('k3',real=True)
    p3=-v*sp.cos(th)**2*sp.sin(th)+k3*sp.sin(th)**3
    assert sp.trigsimp(sp.expand_trig(f2.subs(th,th+sp.pi))-f2)==0
    assert sp.trigsimp(sp.expand_trig(p3.subs(th,th+sp.pi))+p3)==0

    # Necessary wall equations on an annulus retain the unknown endpoint jets.
    dp,dm,cp,cm,hp=sp.symbols('d_plus d_minus c_plus c_minus h_prime',real=True)
    Uplus=r*(lam-dp*hv)+v;Uminus=-r*(lam-dm*hv)+v
    qp=-cp*hv;qm=-cm*hv
    drr=lambda z:sp.diff(z,r)+sp.diff(z,hv)*hp
    right=sp.factor((drr(Uplus*qp)-mu*qp)/(-cp))
    left=sp.factor((drr(Uminus*qm)+mu*qm)/cm)
    expectedRight=(r*(lam-2*dp*hv)+v)*hp+(lam-mu)*hv-dp*hv**2
    expectedLeft=(r*(lam-2*dm*hv)-v)*hp+(lam-mu)*hv-dm*hv**2
    assert sp.expand(right-expectedRight)==0
    assert sp.expand(left-expectedLeft)==0

    x,y=sp.symbols('x y',real=True)
    Cy=sp.Function('C')(y)
    quadraticPsi=sigma*x*y-y*y/2
    quadraticRho=x+Cy
    quadraticTransport=sp.expand((-sp.diff(quadraticPsi,y)+lam*x+v)*sp.diff(quadraticRho,x)
        +(sp.diff(quadraticPsi,x)+lam*y)*sp.diff(quadraticRho,y)-mu*quadraticRho-v)
    finalPsi=quadraticPsi.subs(sigma,lam-mu)
    finalRho=x-y/(2*(lam-mu))
    finalResidual=sp.simplify((-sp.diff(finalPsi,y)+lam*x+v)*sp.diff(finalRho,x)
        +(sp.diff(finalPsi,x)+lam*y)*sp.diff(finalRho,y)-mu*finalRho-v)
    assert finalResidual==0
    result=dict(title='Unknown-angular-function rank-one separation: exact equations and origin obstruction',
        reparameterization='t=log r; a(t)=A(exp(t))*exp(-2t) is a definition, not a power ansatz; rho=r H(t,theta).',
        q=str(q),poissonRho='cos(theta)*(H+H_t)-sin(theta)*H_theta=q',
        exactTransport='(lambda-a F_prime)*(H+H_t)+(2a+a_t)*F*H_theta-mu H+v exp(-t)*(q-1)=0',
        certificates=dict(polarPoisson=str(poisson),transportReduction=str(transport),integralInverse=str(inverse_certificate),retainedC=str(Ccertificate)),
        originTheorem=dict(assumptions='Rank-one psi=A(r)F(theta) on complete half-circles 0<r<R; psi extends C2 to origin and wall, psi(x,0)=0, grad psi(0)=0, q(0)=1; rho is C1 with rho(0)=0; v!=0; steady Poisson and transport hold.',
            angularShapeForcedByTaylor=str(f2),radialNotation='A=r^2 h(r), h(r)->1 after the only scalar normalization of A,F.',
            generalQAfterTaylor=str(qexpected),wall='q(+r,0)=q(-r,0)=h(r); u1(+r,0)=-sigma*r*h(r); u1(-r,0)=sigma*r*h(r).',
            wallResidualSum=str(parity),conclusion='v!=0 forces h=1 identically, hence A=r^2, psi=sigma*x*y-y^2/2 and q=1. Wall antisymmetry then forces sigma=lambda-mu.',
            rhoRemainderEquation=str(quadraticTransport.subs(sigma,lam-mu)),
            targetRates='For lambda>mu>0, C1 regularity at the wall forces rho=x-y/[2(lambda-mu)].',
            exactAffineControlResidual=str(finalResidual),
            targetFailure='q is constant, so it cannot have the required nondegenerate wall maximum; velocity also grows linearly.'),
        annulusWallEquations=dict(assumptions='Only on R0<r<R1; F remains unknown. F(0)=F(pi)=0; d_plus=F_prime(0), d_minus=F_prime(pi); both c_plus=F_second(0),c_minus=F_second(pi) nonzero for division.',
            h='A/r^2, with no prescribed r->0 limit',right=str(expectedRight),left=str(expectedLeft),
            limitation='Necessary differentiated wall equations; lost integration constants must still satisfy undifferentiated wall transport. These do not by themselves certify or exclude the full annular solution.'),
        rankTwoAngularSpaceReview=dict(scope='Exact rank<=2 on full half-circles down to a smooth Cartesian origin; normalized nondegenerate wall peak and v!=0. Does not apply to an annulus or approximate low rank.',
            closedSpace='Use the actual angular row space S=span{psi(r,.)}; it lies in C([0,pi]), has dimension <=2 and is closed in the uniform norm. This avoids assuming that an initially supplied F,G basis is itself regular.',
            twoLimits='Uniform finite Taylor limits give p2=lim psi/r^2 in S and p3=lim(psi-r^2 p2)/r^3 in S.',
            nonzeroCubic='The wall second derivative identity with q_xx(0)=-2 gives P3=-v*x^2*y+k3*y^3; hence P3!=0.',
            independence='Restrictions p2,p3 of homogeneous polynomial degrees 2,3 are independent. If proportional on the upper semicircle, the finite trigonometric-polynomial identity extends to the full circle, contradicting their even/odd behavior under theta->theta+pi.',
            consequence='S=span{p2,p3}; the unknown angular space is thus derived, not preset.',
            finiteRegularity='C6 Taylor jets suffice for P4=u1*r^2*P2, P5=v1*r^2*P3, P6=u2*r^4*P2, by finite polynomial parity. No analyticity or convergence of an infinite Taylor series is needed.',
            radialSmoothnessQualification='Writing psi=U(r^2)P2+V(r^2)P3 is exact for r>0. C-infinity Cartesian regularity yields smooth U,V at s=0 via the smooth-even-function lemma; for finite C6 it is enough and safer to claim only the displayed finite jets.',
            companion='The further exact interior obstruction is recorded in RANK_TWO_JETS.md and rank_two_jet_certificate.json; it is distinct from the rank-one parity proof.'),
        noNumericalData=True)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'separation_equation_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    write_report()
    print('All exact separation and origin-rigidity certificates passed.')
    print('Wall parity identity:',parity)


def write_report():
    text=r'''# 未知角形的真正径角分离：精确方程与原点排除

研究对象仅为

\[
\psi(r,\theta)=A(r)F(\theta),\qquad 0\le\theta\le\pi.
\]

A 与 F 都是真正未知函数。以下不采用数值拟合，不换预设角基，不预设 A 是单一幂。结论分成两种范围：**包含平滑原点的完整 rank-one 分离，在 v≠0 和目标规范下只能给出常 q 的二次流函数；仅在离开原点的半环带分离，必须保留未知角形，不能套用该排除。**

## 1. 对任意未知 A 的直接 Poisson 消元

保持 `u=(-ψ_y,ψ_x)`，于是

\[
u_r=-\frac{A}{r}F',\qquad u_\theta=A'F,
\qquad q=\rho_x=-\left(A''+\frac{A'}r\right)F-\frac{A}{r^2}F''.
\]

为写得紧凑，定义

\[
t=\log r,\quad a(t)=e^{-2t}A(e^t),\quad
\rho(r,\theta)=rH(t,\theta).
\]

这里 `A=r²a(log r)` **只是任意径向函数的重参数化**，没有令 a 为常数，也没有丢弃任何非幂 A。记

\[
\mathcal Q[a,F]=-[a_{tt}+4a_t+4a]F-aF''.
\]

Poisson 的标量恢复方程精确为

\[
\boxed{\cos\theta\,(H+H_t)-\sin\theta\,H_\theta=\mathcal Q[a,F].}\tag{P}
\]

它沿固定 y 的水平线积分，必须保留积分函数 C(y)。若从当前点到 x=0 的水平段全在所研究区域内，令 `y=e^t sinθ`，便有

\[
\begin{aligned}
H(t,\theta)={}&e^{-t}C(e^t\sin\theta)\\
&+\sin\theta\int_\theta^{\pi/2}
\mathcal Q[a,F]\!\left(t+\log\frac{\sin\theta}{\sin\varphi},\varphi\right)
\frac{d\varphi}{\sin^2\varphi}.
\end{aligned}\tag{I}
\]

C(y)=ρ(0,y)。C 项严格在算子 `(P)` 的核中，不能设为零以制造分离。

**环带范围的域限制：**如果 `R0<r<R1` 且 y<R0，水平截面分成左右两块，x=0 不在区域内，不能使用跨过空洞的积分 `(I)`。应在每个连通水平段选一个域内基点 `x_*(y)`，写 `ρ=C_*(y)+∫_{x_*}^x q(s,y)ds`；其等价角积分只沿该段。左右段的积分函数不能未经延拓条件就视为同一个已知 C。后续所有方程都适用于这种局部恢复。

## 2. 未知 A、F、C 的精确分离函数方程

自相似速度的极分量为

\[
V_r=r(\lambda-aF')+v\cos\theta,\qquad
V_\theta=r(2a+a_t)F-v\sin\theta.
\]

用 `(P)` 合并平移项后，完整输运方程恰为

\[
\boxed{
(\lambda-aF')(H+H_t)
+(2a+a_t)F H_\theta-\mu H
+v e^{-t}(\mathcal Q[a,F]-1)=0.
}\tag{T}
\]

将 `(I)` 或其合法局部版本代入 `(T)`，便得到只含未知径向函数 a、未知角函数 F、以及必要积分函数 C 的**精确函数方程**。这里没有投影余量，没有把 A 改成幂律。反过来，只要 `(I)` 在合法区域恢复出的 ρ 满足 `(T)`，Poisson 与原输运同时成立。

角积分中的径向平移 `t↦t+log(sinθ/sinφ)` 说明一般 A 会混合多个径向尺度。可以研究 Mellin 变换，但不能在证明前把变换支持压成一个点。非线性输运还会耦合这些径向成分。

## 3. 包含平滑原点时，角形由 Taylor 定理强制产生

假设 rank-one 分离覆盖完整半圆族 `0<r<R`，ψ 在闭半圆盘上足够光滑，满足

\[
\psi(x,0)=0,\quad\nabla\psi(0,0)=0,\quad
q(0,0)=-\Delta\psi(0,0)=1.
\]

墙面恒等式给 `ψ_xx(0)=0`，故 `ψ_yy(0)=−1`。记 `σ=ψ_xy(0)`，则 Taylor 定理给

\[
\psi(r,\theta)=r^2\left[\sigma\sin\theta\cos\theta-\tfrac12\sin^2\theta\right]+o(r^2).
\]

在 θ=π/2，二次系数为 −1/2，故 F(π/2)≠0，且 A/r² 有非零极限。利用 A、F 唯一的常数倍自由度，规范为 `A/r²→1`。对每个固定 θ 取 r→0，得到

\[
\boxed{F(\theta)=\sigma\sin\theta\cos\theta-\tfrac12\sin^2\theta.}\tag{O}
\]

这不是先选有限 sin 基，而是完整 rank-one 假设和平滑原点共同强制出的未知角形。尤其 `F″(0)=F″(π)=−1`，绝没有把墙面二阶导数设成零。

## 4. 非零漂移再强制所有径向变化消失

写 `A=r²h(r)`，由上节 h→1，

\[
\psi=h(r)(\sigma xy-\tfrac12y^2),\quad
q=h-(r^2h''+5rh')F(\theta).
\]

在两侧墙面点 x=±r：

\[
q(r,0)=q(-r,0)=h(r),\quad
u_1(r,0)=-\sigma rh(r),\quad u_1(-r,0)=\sigma rh(r).
\]

取目标规范 ρ(0,0)=0，记 `S(r)=∫_0^r h(s)ds`，则 `ρ(±r,0)=±S(r)`。将墙面输运分别写在 ±r 并相加，奇函数项全部消去，精确剩下

\[
\boxed{2v[h(r)-1]=0.}
\]

因此 **v≠0 强制 h≡1**，于是

\[
\boxed{A=r^2,\qquad \psi=\sigma xy-\tfrac12y^2,\qquad q\equiv1.}
\]

墙面方程相减再给 `σ=λ−μ`。Poisson 此时仍要求保留 `ρ=x+C(y)`；代入输运得到

\[
(2\lambda-\mu)yC'-\mu C=-y.
\]

对目标 `λ>μ>0`，一般解为

\[
C=-\frac{y}{2(\lambda-\mu)}+K y^{\mu/(2\lambda-\mu)}.
\]

指数位于 (0,1)，C¹ 墙面正则性迫使 K=0。因此唯一这种正规 rank-one 分支为

\[
\boxed{\psi=(\lambda-\mu)xy-\tfrac12y^2,\qquad
\rho=x-\frac{y}{2(\lambda-\mu)}.}
\]

它确实精确满足原方程及原点中心化，但 q 恒定，没有目标 `q_xx(0,0)=−2` 的非退化峰，速度也线性增长。故**包含平滑原点的 rank-one 径角分离与目标分支不相容**。若 v=0，上述奇偶相加不再消去 h，不能声称同一结论。

## 5. 为什么不能把原点定理搬到局部环带

在 `R0<r<R1` 上没有 r→0 的 Taylor 极限，不能推出 `(O)`。只可由墙面无穿透取 F(0)=F(π)=0，而 F′、F″ 的端点值仍是未知的常数。令

\[
d_+=F'(0),\quad d_-=F'(\pi),\quad h=A/r^2.
\]

当两端 F″ 都非零时，对墙面输运沿 r 求一次导数，可得必要条件

\[
[r(\lambda-2d_+h)+v]h'+(\lambda-\mu)h-d_+h^2=0,
\]
\[
[r(\lambda-2d_-h)-v]h'+(\lambda-\mu)h-d_-h^2=0.
\]

这两式保留了未知端点斜率，不等同于原点情形 `d_+=d_-=σ`。还必须补回未微分墙面方程中的积分常数，并满足整个角区间的 `(I),(T)`。它们是精确必要条件，不是已构造的局部解，也不构成对任意环带分离的排除。

## 6. rank-two 的角空间：必须联合推导，不能逐项套用 rank-one

对 `ψ=A(r)F(θ)+B(r)G(θ)`，不能把上述 rank-one 推理逐项施加到两个未知角形。但在当前 v≠0、非退化壁面峰规范下，前三阶 Taylor 信息实际上足以**联合确定整个二维角空间**，过程不需要假设 ψ 解析。

令 `S=span{ψ(r,·):0<r<R}` 为实际角函数行空间。即使最初给出的 F、G 形式不规则，实际各行 ψ(r,·) 仍是连续函数；精确 rank≤2 使 S 是 `C([0,π])` 的至多二维闭子空间。有限阶 Taylor 定理的一致余项给

\[
p_2=\lim_{r\to0}\frac{\psi(r,\cdot)}{r^2}\in S,
\qquad
p_3=\lim_{r\to0}\frac{\psi(r,\cdot)-r^2p_2}{r^3}\in S.
\]

这里 `p_j=P_j(cosθ,sinθ)`。墙面规范和输运的一阶、二阶 x 导数依次强制

\[
P_2=(\lambda-\mu)xy-\tfrac12y^2,
\qquad P_3=-v x^2y+k y^3.
\]

v≠0 保证 P₃ 非零。p₂ 与 p₃ 不可能在上半圆成比例：若成比例，该有限三角多项式恒等式会延拓到整圆，而齐次二次、三次多项式在 `θ↦θ+π` 下分别为偶、奇，产生矛盾。**这里使用的解析性仅是有限 Taylor 多项式固有的性质，不是 ψ 的解析性假设。** 因此 `S=span{p2,p3}`；F、G 的角空间是由方程和精确分离假设推出的，不是预设角基。

对于接下来的排除证书，只需要到六阶的有限 Taylor jets。每个 `p_j∈S`，并具有齐次 j 次多项式的角奇偶性，所以

\[
P_4=u_1 r^2P_2,\qquad P_5=v_1 r^2P_3,
\qquad P_6=u_2 r^4P_2.
\]

这些结论在相应 C⁶ 正则性下成立，无需无穷 Taylor 级数收敛。对 r>0 可以精确写成 `ψ=U(r²)P2+V(r²)P3`。若进一步声称 U、V 在 s=0 是光滑函数，C∞ Cartesian 正则性配合光滑偶函数引理即可保证；只给有限 C⁶ 时，使用上述有限 jets 更准确，且足够完成矛盾。

后续墙面约束及内部输运系数的严格排除见 [RANK_TWO_JETS.md](RANK_TWO_JETS.md) 与 `rank_two_jet_certificate.json`。这仍只针对**覆盖平滑原点的精确 rank≤2**；不排除离开原点的环带分离、近似低秩、或随半径变化的角空间。

精确证书：`separation_equation_certificate.json`。复现脚本：`solver/steady_analysis/streamfunction_separation_equation.py`。所有断言仅使用符号微分、代数化简和边界恒等式；未读取数值源场。
'''
    (OUT/'SEPARATION_EQUATION.md').write_text(text)


if __name__=='__main__':main()
