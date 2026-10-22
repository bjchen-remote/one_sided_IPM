"""Exact certificates for psi=y F(x-Z(y)); no data fitting or numerics."""
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_symbolic_20260929'


def main():
    x,y,z,s=sp.symbols('x y z s',real=True)
    lam,mu,v=sp.symbols('lambda mu v',real=True)
    Z=sp.Function('Z')(y);C=sp.Function('C')(y);Ffun=sp.Function('F')
    fx=Ffun(x-Z)
    psi=y*fx
    u1=-sp.diff(psi,y);u2=sp.diff(psi,x)
    aa=2*sp.diff(Z,y)+y*sp.diff(Z,y,2)
    bb=y*(1+sp.diff(Z,y)**2)
    rho=aa*fx-bb*sp.diff(fx,x)+C
    q=-sp.diff(psi,x,2)-sp.diff(psi,y,2)
    poisson_certificate=sp.simplify(sp.diff(rho,x)-q)
    assert poisson_certificate==0
    full=(u1+lam*x+v)*sp.diff(rho,x)+(u2+lam*y)*sp.diff(rho,y)-mu*rho-v
    dd=lam*(Z-y*sp.diff(Z,y))+v
    compact=((lam*(x-Z)+dd-fx)*(aa*sp.diff(fx,x)-bb*sp.diff(fx,x,2))
        +(mu*bb-lam*y*sp.diff(bb,y))*sp.diff(fx,x)-y*sp.diff(bb,y)*sp.diff(fx,x)**2
        +y*sp.diff(aa,y)*fx*sp.diff(fx,x)+(lam*y*sp.diff(aa,y)-mu*aa)*fx
        +y*sp.diff(C,y)*sp.diff(fx,x)+lam*y*sp.diff(C,y)-mu*C-v)
    compact_certificate=sp.simplify(sp.expand(full-compact))
    assert compact_certificate==0

    # Formal autonomous-polynomial reduction, without treating F,F',F'' as
    # algebraically independent: insert the relation before any matching.
    a,b,ap,bp,c,cp,d=sp.symbols('a b a_prime b_prime C C_prime d',real=True)
    p0,p1,p2,p3=sp.symbols('p0 p1 p2 p3',real=True)
    P=p0+p1*s+p2*s**2+p3*s**3
    qP=P*(a-b*sp.diff(P,s))
    E=(lam*z+d-s)*qP+(mu*b-lam*y*bp)*P-y*bp*P**2+y*ap*s*P+(lam*y*ap-mu*a)*s+y*cp*P+lam*y*cp-mu*c-v
    assert sp.simplify(sp.diff(E,z)-lam*qP)==0
    w,kappa=sp.symbols('w kappa',real=True,nonzero=True)
    f0,z0=sp.symbols('f0 z0',real=True)
    powercert=[]
    for n in (2,3):
        PP=kappa*(s-f0)**n
        EN=(lam*z+d-s)*PP*(a-b*sp.diff(PP,s))+(mu*b-lam*y*bp)*PP-y*bp*PP**2+y*ap*s*PP+(lam*y*ap-mu*a)*s+y*cp*PP+lam*y*cp-mu*c-v
        inverse=z0+w**(1-n)/(kappa*(1-n))
        poly=sp.Poly(sp.expand(EN.subs({s:w+f0,z:inverse})),w)
        top=sp.factor(poly.coeff_monomial(w**(2*n)))
        assert sp.simplify(top-kappa**2*(n*b-y*bp))==0
        assert sp.simplify(sp.diff(inverse,w)-1/(kappa*w**n))==0
        k,kp,k0=sp.symbols('k k_prime k0',real=True)
        wall=sp.factor((n*b-y*bp).subs({b:y*(1+k*k),bp:1+k*k+2*y*k*kp})/y)
        walllimit=sp.factor(wall.subs({y:0,k:k0}))
        assert sp.simplify(walllimit-(n-1)*(1+k0*k0))==0
        powercert.append(dict(degree=n,inverse=str(inverse),residualPolynomial=str(poly.as_expr()),
            coefficientsDescending=[str(sp.factor(t)) for t in poly.all_coeffs()],
            highestCoefficient=str(top),wallCoefficientDividedByKappaSquaredY=str(wall),
            wallLimit=str(walllimit)))
    r1,r2,r3=sp.symbols('r1 r2 r3')
    simpleResidues={
        'quadratic_distinct':str(sp.simplify(sp.residue(1/(kappa*(s-r1)*(s-r2)),s,r1))),
        'cubic_distinct':str(sp.simplify(sp.residue(1/(kappa*(s-r1)*(s-r2)*(s-r3)),s,r1))),
        'cubic_double_plus_simple':str(sp.simplify(sp.residue(1/(kappa*(s-r1)*(s-r2)**2),s,r1)))}

    # Affine translation: an exact exclusion for arbitrary non-affine F.
    k,h,C0,C1,D0=sp.symbols('k h C0 C1 D0',real=True)
    F=sp.Function('F')(z);Fp=sp.diff(F,z);Fpp=sp.diff(F,z,2)
    constant=2*k*((lam*z+D0-F)*Fp-mu*F)-mu*C0-v
    linear=-h*(lam*z+D0-F)*Fpp+h*(mu-lam)*Fp-h*Fp**2+C1*Fp+(lam-mu)*C1
    eliminated=sp.simplify(sp.expand(linear+h*sp.diff(constant,z)/(2*k)))
    expected=-2*h*Fp**2+C1*Fp+(lam-mu)*C1
    assert sp.simplify(eliminated-expected)==0
    # Stronger target-specific obstruction, valid for arbitrary smooth Z,F.
    # Wall values depend only on k=Z'(0), regardless of higher Z derivatives.
    wall=2*k*(lam*z+v-F)*Fp-2*mu*k*F-v
    f1,f2,f3,qxx=sp.symbols('f1 f2 f3 qxx',real=True)
    jetmap={F:sp.Integer(0),Fp:f1,Fpp:f2,sp.diff(F,z,3):f3,z:sp.Integer(0)}
    walljets=[sp.expand(sp.diff(wall,z,j)).xreplace(jetmap) for j in (0,1,2)]
    normalizedjets=[sp.simplify(t.subs({f1:1/(2*k),f2:0,f3:qxx/(2*k)})) for t in walljets]
    assert normalizedjets[0]==0
    assert sp.simplify(normalizedjets[1]-(lam-mu-1/(2*k)))==0
    assert sp.simplify(normalizedjets[2]-v*qxx)==0
    affine_E=2*k*((lam*z+D0-F)*Fp-mu*F)+y*(-h*(lam*z+D0-F)*Fpp+h*(mu-lam)*Fp-h*Fp**2)+y*sp.diff(C,y)*Fp+lam*y*sp.diff(C,y)-mu*C-v
    yy=sp.expand(sp.diff(affine_E,y,2))
    yyexpected=sp.diff(y*sp.diff(C,y),y,2)*Fp+sp.diff(lam*y*sp.diff(C,y)-mu*C,y,2)
    assert sp.simplify(yy-yyexpected)==0

    # Positive-q affine control: it solves the equations, but has linear
    # velocity and constant q, hence is outside the desired decaying branch.
    K=lam-mu
    psia=K*x*y-y*y/2
    rhoa=x-y/(2*K)
    ua=-sp.diff(psia,y);va=sp.diff(psia,x)
    control=sp.simplify((ua+lam*x+v)*sp.diff(rhoa,x)+(va+lam*y)*sp.diff(rhoa,y)-mu*rhoa-v)
    assert control==0
    assert sp.simplify(-sp.diff(psia,x,2)-sp.diff(psia,y,2)-sp.diff(rhoa,x))==0
    certificate=dict(title='Exact translated-streamfunction exclusions',ansatz='psi=y F(z), z=x-Z(y)',
        parameters='a=2 Z_prime+y Z_second; b=y(1+Z_prime^2); d=lambda(Z-y Z_prime)+v',
        rho='a F-b F_prime+C(y)',poissonCertificate=str(poisson_certificate),
        compactTransportCertificate=str(compact_certificate),
        compactTransport='(lambda*z+d-F)(a*F_prime-b*F_second)+(mu*b-lambda*y*b_prime)*F_prime-y*b_prime*F_prime^2+y*a_prime*F*F_prime+(lambda*y*a_prime-mu*a)*F+y*C_prime*F_prime+lambda*y*C_prime-mu*C-v',
        autonomous=dict(E=str(sp.expand(E)),zCoefficient=str(sp.factor(sp.diff(E,z))),
            rationalPrimitiveObstruction='If lambda!=0 and q is not identically zero, E=0 solves z=R(F) rationally; differentiating gives R_prime=1/P. A simple complex zero of P yields a nonzero residue, impossible for a rational derivative.',
            simpleRootResidues=simpleResidues,purePowerCertificates=powercert,
            degree1='If P=p1 F+p0 with p1!=0, rational primitive obstruction forces a-b*p1=0 at every y>0, hence q=0.',
            degree0='P constant gives affine F; not a non-affine profile.'),
        wallPeakCertificate=dict(assumptions='Any smooth real F,Z; centre u1(0,0)=0; q(0,0)>0, q_x(0,0)=0, q_xx(0,0)!=0; v!=0.',
            identity='Wall rho=2k F+C0, q=2k F_prime, U=lambda*x+v-F; differentiate U*q-mu*rho-v twice: U_second*q+2 U_prime*q_prime+U*q_second-mu*q_prime=0.',
            peakReduction='q_prime(0)=0 implies F_second(0)=U_second(0)=0, hence v*q_second(0)=0: contradiction.',
            normalizedAssumptions='For the explicit jet certificate only: Z(0)=F(0)=rho(0,0)=0, q(0,0)=1, q_prime(0)=0; k=Z_prime(0)!=0.',
            residualJetsBeforeNormalization=[str(t) for t in walljets],
            residualJetsAfterNormalization=[str(t) for t in normalizedjets],
            normalizedCurvatureMinusTwoResidual=str(sp.simplify(normalizedjets[2].subs(qxx,-2))),
            strongerLocalUniqueness='With q(0)=1,rho(0)=F(0)=0, first differentiated wall equation gives F_prime(0)=lambda-mu and 2k(lambda-mu)=1. Since v!=0, F_prime=(mu F+v/(2k))/(lambda*x+v-F) is regular near 0 and uniquely selects F=(lambda-mu)x. Thus even a degenerate isolated peak cannot arise locally.'),
        affineTranslation=dict(assumptions='Z=ky+z0, k!=0; C smooth up to y=0; F_prime nonconstant on a connected interval',
            secondYDerivative=str(yy),consequence='(y C_prime)_second=0 implies C=C0+C1*y+C_log*log(y); smooth wall forces C_log=0.',
            constantYTerm=str(constant),linearYTerm=str(linear),
            elimination='linear_y_term + h/(2k) * derivative_z(constant_y_term)',
            eliminatedExpression=str(expected),h='1+k^2>0',
            contradiction='A continuous nonconstant F_prime cannot take values only in the at-most-two roots of a fixed nonzero quadratic.'),
        affineControl=dict(psi=str(psia),rho=str(rhoa),q='1',u1=str(ua),u2=str(va),
            transportCertificate=str(control),limitation='Linear velocity growth and nondecaying q; does not satisfy the intended sublinear far field.'),
        conclusion='For v!=0, the entire psi=y F(x-Z(y)) class fails the centred nondegenerate positive-wall-peak compatibility, without restrictions on F beyond smoothness. Independently, even without this target peak normalization, non-affine smooth-wall solutions are excluded for affine Z with k!=0 and arbitrary F, or arbitrary smooth Z with autonomous polynomial F_prime=P(F), degree<=3 and lambda!=0.',
        limits='The target-specific exclusion requires a centred wall peak and nonzero v. Without those assumptions, the two independent theorems have the stated narrower scopes. Additional y amplitude or multiple generators are not excluded here.')
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'translated_generator_certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    write_report()
    print('All exact symbolic assertions passed.')
    print(json.dumps({'purePowerTopCoefficients':[q['highestCoefficient'] for q in powercert],
                      'affineTranslationElimination':str(expected),'output':str(OUT)},indent=2))


def write_report():
    report=r'''# 平移流函数生成式的精确闭合检验

这次只做符号运算，不拟合 F，不把 F、F′、F″ 当成相互独立。最直接的结果是：**在所需中心化壁面正峰规范下，v≠0 已经排除整个 `ψ=yF(x−Z(y))` 类，不论 F 是何种母函数、Z 是否非线性。** 另保留两条不依赖目标峰规范的排除定理：线性横移 Z 排除所有非仿射 F；非线性 Z 配次数≤3自主多项式 F′ 也被排除。

## 0. 目标壁面正峰给出的最短排除证书

平移 F 的自变量后可令 Z(0)=0。中心化 `u1(0,0)=0` 给 F(0)=0。记 k=Z′(0)，则任意光滑非线性 Z 都给出同一墙面关系：

\[
\rho(x,0)=2kF(x)+C_0,\quad q(x,0)=2kF'(x),\quad
U(x)=u_1(x,0)+\lambda x+v=-F(x)+\lambda x+v.
\]

墙面正值 q(0,0)>0 强制 k≠0。墙面输运恰为 `Uq=μρ+v`。对 x 求两次：

\[
U''q+2U'q'+Uq''=\mu q'.
\]

若原点为 q 的峰，`q′(0)=0`，则 `F″(0)=0`、`U″(0)=0`。又 `U(0)=v`，故

\[
\boxed{v\,q''(0)=0.}
\]

所需非退化峰 `q″(0)<0` 与 v≠0 矛盾。项目的规范 `q(0)=1,q′(0)=0,q″(0)=−2` 使二阶残差直接等于 **−2v**。这个短证不限制 F 的微分代数关系，也不限制 Z 的高阶导数；基本矛盾甚至不需要规定密度积分常数 C₀。

在完整规范 `ρ(0)=F(0)=0,q(0)=1` 下还能排除退化的孤立峰：一次微分给 `F′(0)=λ−μ`、`2k(λ−μ)=1`，墙面 ODE 为

\[
F'=\frac{\mu F+v/(2k)}{\lambda x+v-F}.
\]

因 v≠0，它在原点正规，局部唯一解就是 `F=(λ−μ)x`。于是 q 在原点邻域恒为 1，不能形成孤立峰。后两节是独立的较宽条件证书，并非此短证所需附加假设。

## 1. 生成式与完整输运恒等式

采用项目符号

\[
\psi=yF(z),\quad z=x-Z(y),\qquad u=(-\psi_y,\psi_x),\quad q=-\Delta\psi=\rho_x.
\]

记

\[
a=2Z'+yZ'',\quad b=y(1+Z'^2),\quad d=\lambda(Z-yZ')+v.
\]

直接得到

\[
u_1=-F+yZ'F',\quad u_2=yF',\qquad
q=aF'-bF'',\quad \rho=aF-bF'+C(y).
\]

Poisson 关系是恒等式。输运方程
\((u+\lambda(x,y)+ve_1)\cdot\nabla\rho=\mu\rho+v\)
恰等价于

\[
\begin{aligned}
0=\mathcal E={}&(\lambda z+d-F)(aF'-bF'')
+(\mu b-\lambda yb')F'-yb'(F')^2\\
&+ya'FF'+(\lambda ya'-\mu a)F
+yC'F'+\lambda yC'-\mu C-v.
\end{aligned}
\]

脚本独立从 ψ、ρ 求偏导，并精确验证原残差减此式为 0。

## 2. 定理 A：线性横移不能支持非仿射 F

假设 `Z=ky+z0`、`k≠0`，C 光滑到墙面。F 可为任意充分光滑函数，**不附加自主 ODE 或代数关系**。若 F 非仿射，则 F′ 在某个连通区间非恒定。

此时 `a=2k`、`b=hy`、`h=1+k²>0`、`d=λz0+v`。残差对 y 求两次，只有 C 项留下：

\[
(yC')''F'+(\lambda yC'-\mu C)''=0.
\]

利用 F′ 非恒定，得到 `(yC′)″=0`，故 `C=C0+C1 y+C_log log y`。墙面光滑排除对数项，C 必须线性。

残差的 y⁰、y¹ 系数分别为

\[
A_0=2k[(\lambda z+d-F)F'-\mu F]-\mu C_0-v,
\]
\[
A_1=-h(\lambda z+d-F)F''+h(\mu-\lambda)F'
-h(F')^2+C_1F'+(\lambda-\mu)C_1.
\]

现在对 A₀ 求 z 导数并与 A₁ 组合，F″ 精确消去：

\[
\boxed{A_1+\frac h{2k}A_0'
=-2h(F')^2+C_1F'+(\lambda-\mu)C_1=0.}
\]

一个连续非恒定 F′ 不可能只取固定非零二次多项式的至多两个根。因此 F′ 必为常数，与非仿射假设矛盾。

若 k=0，则墙面 q=0，已经不属于所需墙上正峰分支。定理不宣称排除这一零墙面值分支的所有离墙解。

## 3. 定理 B：非线性 Z 也救不了低次自主多项式关系

假设 λ≠0，Z 在墙面附近实值、至少 C³，F 非恒定，且

\[
F'=P(F),\qquad P\in\mathbb R[F],\quad \deg P\le3.
\]

将关系及 `F″=P′P` **先代入** 完整恒等式。得到

\[
\mathcal E=\lambda z\,P(F)[a-bP'(F)]+S(F;y)=0,
\]

其中 S 在固定 y 时是 F 的多项式，`q=P(F)[a−bP′(F)]`。

在任意 q 非恒零的固定高度，可解出 `z=R(F)`，R 是有理函数。沿非恒定 F 的区间求导，必须有

\[
\boxed{R'(F)=1/P(F).}
\]

有理函数的导数不含简单极点留数。如果 P 有简单复根，则 1/P 在该根的留数非零，不可能等于 R′。这不是把 F、F′ 独立匹配，而是利用实际自主关系导出的微分代数必要条件。

对于二次或三次多项式，除纯重根情形外一定有简单复根。因此只剩

\[
P(F)=\kappa(F-f)^n,\qquad n=2\text{ 或 }3,\quad\kappa\ne0.
\]

令 `w=F−f`。非恒定支的精确逆关系为

\[
z=z_0+\frac{w^{1-n}}{\kappa(1-n)}.
\]

代回完整输运式后变成 w 的多项式，其最高次 w^(2n) 系数恰为

\[
\boxed{\kappa^2(nb-yb').}
\]

所以必要条件是 `yb′=nb`。代入 `b=y(1+Z′²)`，并在 y>0 除以 y：

\[
(n-1)(1+Z'^2)-2yZ'Z''=0.
\]

令 y→0，墙边光滑要求

\[
(n-1)(1+Z'(0)^2)=0,
\]

而 n>1、Z′(0) 实值，左侧严格正。矛盾。等价积分形式是 `1+Z′²=K y^(n−1)`，同样不能光滑抵达墙面。

次数 1 时，如果 P 非恒定，其简单根带来同一有理原函数障碍；唯一逃逸是 `a−bP′=0`，但这直接令 q=0。次数 0 只产生仿射 F。故此类不存在所需非仿射正 q 分支。

这个排除甚至是局部光滑墙边结论，不依赖远场。纯二次/三次剩余支分别是带实极点的倒一次式与倒平方根式，它们还另有全横轴光滑性障碍，但证明并未依赖该障碍。

## 4. 一个精确对照，说明并非所有解都被错误排除

令 κ=λ−μ≠0，取

\[
\psi=\kappa xy-\frac12y^2,\qquad
\rho=x-\frac{y}{2\kappa},\qquad q=1.
\]

它对应 `Z=y/(2κ), F(z)=κz`，且

\[
u_1=-\kappa x+y,\qquad u_2=\kappa y.
\]

它对任意 v（包括 v≠0）精确满足 Poisson 与输运，并满足墙面无穿透和原点中心化。脚本精确验证两残差均为 0。它的速度线性增长、q 不衰减，所以不是目标慢尾、次线性速度分支。

## 5. 结论的边界

对于当前 v≠0 的中心化壁面峰目标，已经排除整个 `ψ=yF(x−Z(y))` 类。若放弃该目标峰或令 v=0，则只能使用定理 A/B 各自明确的较窄条件，不能将结论扩大到任意非自主 F 配非线性 Z。**没有排除独立幅度 A(y) 或多个流函数生成项。** 本次不为绕过排除而继续扩展母式。

证书为 `translated_generator_certificate.json`；复现脚本 `solver/steady_analysis/translated_streamfunction_symbolic.py` 使用 SymPy 精确恒等式断言，未调用数值拟合或 PDE 时间推进。
'''
    (OUT/'TRANSLATED_GENERATOR.md').write_text(report)


if __name__=='__main__':main()
