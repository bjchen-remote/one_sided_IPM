"""Low-complexity structure probes for saved IPM steady profiles.

These are model-discrimination experiments, not corrections of the PDE solution.
The normal fits use the measured wall F, and therefore are NOT closed solutions.
Run: PYTHONPATH=/private/tmp/ipm_steady_deps OPENBLAS_NUM_THREADS=1 python ...
"""
from pathlib import Path
import json
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import least_squares, brentq
import sympy as sp

from fixed_point import fd_matrix

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'runs/steady_20260928_tau20'
OUT = ROOT / 'runs/structure_20260928'
NAMES = ('L4_129x65', 'L4_257x129', 'L6_193x97',
         'L8_257x129_continued')


def error_metrics(error, train):
    return dict(train_linf=float(np.max(abs(error[train]))),
                holdout_linf=float(np.max(abs(error[~train]))),
                rms=float(np.sqrt(np.mean(error**2))))


def svd_metrics(z):
    sigma = np.linalg.svd(z, compute_uv=False)
    total = np.linalg.norm(sigma)
    return dict(normalized_singular_values=(sigma[:6]/total).tolist(),
                rank1_relative_l2=float(np.linalg.norm(sigma[1:])/total),
                rank2_relative_l2=float(np.linalg.norm(sigma[2:])/total))


def polynomial_probe(values, a, b, train, degrees=(2, 3, 4)):
    result = []
    for degree in degrees:
        powers = [(i, j) for j in range(degree+1)
                  for i in range(degree+1-j)]
        matrix = np.column_stack([(a**i*b**j).ravel() for i, j in powers])
        coeff = np.linalg.lstsq(matrix[train.ravel()],
                               values[train], rcond=None)[0]
        error = (matrix@coeff).reshape(values.shape)-values
        result.append(dict(degree=degree, terms=powers, coefficients=coeff.tolist(),
                           **error_metrics(error, train)))
    return result


def one_profile(name):
    data = np.load(SOURCE/name/'state.npz')
    x, y, Q = data['x'], data['y'], data['Q']
    lam, mu, v, b = data['rates']
    dx, dy = fd_matrix(x, 1, 7), fd_matrix(y, 1, 7)
    qx, qy = Q@dx.T, dy@Q
    qxx, qxy = Q@fd_matrix(x, 2, 7).T, dy@qx
    qyy = fd_matrix(y, 2, 7)@Q
    i0 = int(np.argmin(abs(x)))

    # A genuinely common density rectangle, with NO x extrapolation.
    levels = np.linspace(-.7, .5, 81)
    heights = y[y <= 2]
    inverse = np.array([PchipInterpolator(Q[i], x, extrapolate=False)(levels)
                        for i in range(len(heights))])
    assert np.all(np.isfinite(inverse))
    S, Y = np.meshgrid(levels, heights)
    train = Y <= 1
    wall = inverse[0][None, :]
    zero_curve = np.array([PchipInterpolator(Q[i], x)(0.)
                           for i in range(len(heights))])
    displacement = inverse-wall
    centered = displacement-zero_curve[:, None]
    result = dict(name=name, rates=data['rates'].tolist(), alpha=float(mu/lam),
                  density_shift=float(v/mu),
                  inverse_domain=dict(Q=[-.7, .5], y=[0., 2.], train_y=[0., 1.]),
                  inverse_displacement_svd=svd_metrics(displacement[1:]),
                  inverse_centered_svd=svd_metrics(centered[1:]),
                  inverse_polynomial=polynomial_probe(inverse, S, Y, train))
    normal_models = []
    for power in (.5, 2/3, None):
        def model(parameters):
            k, a = parameters[:2]
            p = parameters[2] if power is None else power
            return wall+k*Y+S*((1+a*Y)**p-1)
        if power is None:
            initial, bounds = [1.4, 1.6, .5], ([0, 0, .05], [4, 20, 2])
        else:
            initial, bounds = [1.4, 1.6], ([0, 0], [4, 20])
        fit = least_squares(lambda p: (model(p)-inverse)[train], initial,
                            bounds=bounds)
        normal_models.append(dict(power=power if power is not None else 'free',
                                  parameters=fit.x.tolist(),
                                  **error_metrics(model(fit.x)-inverse, train)))
    result['normal_models_with_measured_wall'] = normal_models
    k = float(heights@zero_curve/(heights@heights))
    result['zero_density_curve'] = dict(linear_slope=k,
        line_linf=float(np.max(abs(zero_curve-k*heights))))

    # Test an affine change of the physical plane followed by analytic f(z).
    # Norm is componentwise, so it does not hide the size through cancellations.
    core = (abs(x)[None, :] <= 2) & (y[:, None] <= 2)
    htrain = (abs(x)[None, :] <= 1) & (y[:, None] <= 1) & (y[:, None] >= .125)
    htest = core & (y[:, None] >= .125) & ~htrain
    matrix = np.column_stack([qxy[htrain], qyy[htrain]])
    c, bb = np.linalg.lstsq(matrix, -qxx[htrain], rcond=None)[0]
    residual = qxx+c*qxy+bb*qyy
    def component_error(mask):
        return float(np.linalg.norm(residual[mask]) / np.sqrt(np.sum(
            qxx[mask]**2+(c*qxy[mask])**2+(bb*qyy[mask])**2)))
    result['constant_metric_harmonic'] = dict(c=float(c), b=float(bb),
        elliptic_discriminant=float(4*bb-c*c),
        train_component_relative_l2=component_error(htrain),
        holdout_component_relative_l2=component_error(htest),
        holdout_linf=float(np.max(abs(residual[htest]))),
        euclidean_laplacian_core_rms=float(np.sqrt(np.mean((qxx[core]+qyy[core])**2))))

    # Wall inclination K=X_Q-independent normal derivative X_y, with chain rule.
    K = -qy[0]/qx[0]
    K1 = (dx@K)/qx[0]
    K2 = (dx@K1)/qx[0]
    sample = np.linspace(-.7, .5, 121)
    kval = PchipInterpolator(Q[0], K)(sample)
    fits = []
    for degree in (1, 2):
        matrix = np.column_stack([sample**j for j in range(degree+1)])
        coeff = np.linalg.lstsq(matrix, kval, rcond=None)[0]
        fits.append(dict(degree=degree, coefficients=coeff.tolist(),
                         linf=float(np.max(abs(matrix@coeff-kval)))))
    xs = brentq(lambda xx: CubicSpline(x, Q[0])(xx)+v/mu, x[0], 0)
    jets = [float(CubicSpline(x, field)(xs)) for field in (K, K1, K2)]
    result['wall_inclination'] = dict(fit_Q_domain=[-.7, .5], polynomial_fits=fits,
        at_peak=[float(K[i0]), float(K1[i0])],
        peak_identity=float(v*K1[i0]+2*(lam-mu)*K[i0]),
        stagnation_x=xs, stagnation_jets_0_1_2=jets,
        resonance_test=float(2*jets[0]*jets[2]/(5*jets[1]**2)),
        r=float(lam/mu-1))

    # Canonical area coordinate t=integral X_s dy, s=Q+v/mu. Translation of
    # levels makes no difference to the polynomial spaces or rank tests.
    J = np.array([PchipInterpolator(Q[i], 1/qx[i], extrapolate=False)(levels)
                  for i in range(len(heights))])
    area = cumulative_trapezoid(J, heights, axis=0, initial=0)
    ts = np.linspace(0, np.min(area[-1]), 81)
    canonical_y = np.column_stack([PchipInterpolator(area[:, j], heights)(ts)
                                   for j in range(len(levels))])
    canonical_x = np.column_stack([PchipInterpolator(area[:, j], inverse[:, j])(ts)
                                   for j in range(len(levels))])
    SS, TT = np.meshgrid(levels, ts)
    ttrain = TT <= ts[-1]/2
    result['canonical_area'] = dict(tmax=float(ts[-1]),
        X_svd=svd_metrics(canonical_x[1:]-canonical_x[0]),
        Y_svd=svd_metrics(canonical_y[1:]),
        X_polynomial=polynomial_probe(canonical_x, SS, TT, ttrain),
        Y_polynomial=polynomial_probe(canonical_y, SS, TT, ttrain))
    return result


def symbolic_exclusions():
    q, k, a, mu, lam, d = sp.symbols('q k a mu lam d', real=True)
    J, K = 1+q*q, k+a*q/2
    R = sp.expand((mu*(q+d)*sp.diff(K, q)+2*(lam-mu)*K)*J
                  -2*mu*(q+d)*sp.diff(J, q)*K-1)
    assert sp.expand(R.subs({a: 0, d: 0, lam: 3*mu, k: 1/(4*mu)})) == 0
    # General K=A+Bz+Cz², r=2: extract the resonance condition.
    z, A, B, C, j0, j1, j2 = sp.symbols('z A B C j0 j1 j2')
    KK, JJ = A+B*z+C*z*z, j0+j1*z+j2*z*z
    res = sp.expand(2*z*KK*sp.diff(JJ, z)-(z*sp.diff(KK, z)+4*KK)*JJ)
    order1 = res.coeff(z, 1)
    order2 = sp.factor(res.coeff(z, 2).subs(j1, -5*B*j0/(2*A)))
    assert order1 == -2*A*j1-5*B*j0
    assert sp.simplify(order2-3*j0*(5*B**2-4*A*C)/(2*A)) == 0
    return dict(cubic_wall_sqrt_extension_wall_residual=str(R),
        residual_coefficients_descending=[str(c) for c in sp.Poly(R, q).all_coeffs()],
        only_nondegenerate_solution='a=0, d=0, lambda=3*mu, k=1/(4*mu)',
        affine_inclination_recurrence='2 A (n-r) j_n + B (2n-3-2r) j_(n-1)=0',
        r='lambda/mu-1',
        general_resonance_at_r2='2 K(0) K_second(0) = 5 K_prime(0)^2',
        affine_consequence='B=0 if J is C2 at zeta=0 and r=2',
        hypergeometric_particular='J=sqrt(1+B*zeta/A)/(2*mu*r*A) * hyp2f1(3/2,-r;1-r;-B*zeta/A)',
        homogeneous='C*zeta^r*sqrt(A+B*zeta); choose real branches locally',
        caveat='At integer r the hypergeometric representation needs a resonant limit; logarithms generally occur.')


def write_report(results):
    rows = []
    for r in results['profiles']:
        half, twothirds, free = r['normal_models_with_measured_wall']
        rows.append(f"| {r['name']} | {half['parameters'][0]:.5f} | {half['parameters'][1]:.5f} | "
                    f"{free['parameters'][2]:.5f} | {half['train_linf']:.5f} | {half['holdout_linf']:.5f} |")
    report = r'''# 从离散解寻找低复杂度二维结构

这次只做结构识别，不求更小的离散方程残差，也不进行新的 PDE 时间推进。数据为四个已存稳态。结果包含新的低参数几何线索、对它的严格排除，以及一条超几何候选方向；没有找到通过全部方程和全局条件的闭式解。

## 1. 等密度线的近分离结构

记 `x=X(Q,y)`，直接反演每个高度的单调 Q，不向原网格外插值。公共检验矩形为 `−0.7≤Q≤0.5, 0≤y≤2`。壁面记 `F(Q)=X(Q,0)`。以下模型仅在 `y≤1` 拟合，在 `1<y≤2` 外推：

\[
X(Q,y)\approx F(Q)+ky+Q\big[(1+ay)^p-1\big].
\]

F 使用测得的整条壁面函数，所以这不是有限参数精确解。新发现是自由拟合的 p 四组均接近 1/2；取 1/2 时有一个平方根展宽的简单几何模板。

| 数据 | k（固定 p=1/2） | a | 自由拟合 p | 固定 p=1/2 训练最大 x 误差 | 外推最大 x 误差 |
|---|---:|---:|---:|---:|---:|
''' + '\n'.join(rows) + r'''

误差单位是归一化 x，绝非 PDE 残差。这种近分离不是严格分离：减掉壁面和 Q=0 等值线后，位移矩阵仍有非零的第二、第三奇异值。拟合只用于定位应该寻找的坐标结构。

另一条视觉线索是 Q=0 的等密度线非常接近直线 `x=ky`，但壁面真实切线与拟合直线仍相差约几个百分点。低阶二维多项式虽能拟合近壁面，却在留出高度明显失效。精确的面积坐标 `t=∫X_Q dy` 也没有把 X、Y 化为低阶多项式：三、四阶多项式在留出高度的外推误差反而增大。

## 2. 一个漂亮新公式被精确排除

将上述平方根结构与原先的精确三次壁面结合：

\[
X(Q,y)=Q^3/3+Q\sqrt{1+ay}+ky.
\]

它的壁面 `J=F'=1+Q²`，离壁斜率 `K=X_y=k+aQ/2`。令 `ζ=Q+d`，`d=v/μ`。任何精确稳态都必须满足

\[
\mu\zeta K'+[2(\lambda-\mu)-2\mu\zeta J'/J]K=1/J.
\]

代入后乘 J，Q³、Q²、Q¹、Q⁰ 系数依次为

\[
a\lambda-5a\mu/2,\quad
-3ad\mu/2+2k\lambda-6k\mu,\quad
a\lambda-a\mu/2-4dk\mu,\quad
ad\mu/2+2k\lambda-2k\mu-1.
\]

对 μ>0，实参数联立只剩 `a=0,d=0,λ=3μ,k=1/(4μ)`，退化回已有报告排除的纯斜移延拓。因此平方根展宽 + 精确三次壁面的漂亮表达式不能是目标精确解。这个结论比拟合误差强，不依赖网格精度。

## 3. 仿射离壁斜率会导向超几何函数，但在 1/3 发生共振

若放开壁面 F，只假设初始等密度斜率 `K(ζ)=A+Bζ`，则 `J=F'` 精确满足一阶线性 ODE：

\[
2\mu\zeta(A+B\zeta)J'
-[\mu B\zeta+2(\lambda-\mu)(A+B\zeta)]J=-1.
\]

设 `r=λ/μ−1`，在适当的局部实分支、非共振参数下，一个特解为

\[
J(\zeta)=\frac{\sqrt{1+B\zeta/A}}{2\mu rA}
{}_2F_1\!\left(\frac32,-r;1-r;-\frac{B\zeta}{A}\right),
\]

并可加齐次项 `C ζ^r sqrt(A+Bζ)`。它说明“非多项式壁面 + 简单几何闭合”值得检查，但还没有满足二维 Poisson 及全局尾场；非零 B 的仿射 K 也不能在整个双向无穷密度轴上保持正。

其驻点 Taylor 系数 `J=Σj_n ζ^n` 满足

\[
2A(n-r)j_n+B(2n-3-2r)j_{n-1}=0,
\qquad j_0=\frac1{2\mu r A}.
\]

若 μ/λ 恰为 1/3，即 r=2，n=1 给 `j1=−5B j0/(2A)`，n=2 又强制 `−3B j1=0`。所以在驻点至少要求 J 为 C² 时，非恒定仿射 K 被排除；否则一般出现 `ζ² log|ζ|` 共振。

更一般、不假设 K 仿射，r=2 的光滑相容性给出新必要条件：

\[
\boxed{2K(0)K''(0)=5[K'(0)]^2.}
\]

局部函数 `K=A(1−cζ)^(−2/3)` 的导数比例恰好满足此条件，因此可以作为搜索正值非多项式闭合的线索；它本身并不是已构造出的全局解。

## 4. 复解析 / 常系数椭圆结构检验

Q 本身并不谐：公共物理核心的 ΔQ 均方根约 0.48–0.58。拟合 `Qxx+c Qxy+b Qyy=0` 得到正定椭圆型系数，提示可以尝试经过剪切伸缩的复变量结构。然而以各项平方和作尺度，训练误差仍约 3.5%–4.4%，留出核心为 7.5%–9.4%；所以常系数各向异性谐函数也没有达到精确结构证据。单报 Hessian 最小奇异值约 2% 会使印象过于乐观，JSON 同时保留完整的归一化残差。

## 5. 稳健性与复现

峰点的精确必要条件 `μd K'(d)+2(λ−μ)K(d)=1` 用七点独立导数得到 1.000075–1.000411。d 在四组间为 0.834–0.934，尚不足认定常数。高阶驻点 jets 对网格和 r 接近整数的程度敏感，不能拿它们直接猜精确常数。

详细量化、SVD、留出区域误差、壁面斜率和驻点 jets 均存于 `data_structure.json`，复现脚本为 `solver/steady_analysis/structure_data.py`。最有价值的下一步是寻找正值非仿射 K 或有限个非多项式等密度模态，并立即使用壁面驻点共振及全二维 Poisson 作代数筛选。
'''
    (OUT/'DATA.md').write_text(report)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    results = dict(symbolic=symbolic_exclusions(),
                   profiles=[one_profile(name) for name in NAMES])
    (OUT/'data_structure.json').write_text(json.dumps(results, indent=2)+'\n')
    write_report(results)
    for result in results['profiles']:
        print(result['name'], result['normal_models_with_measured_wall'])
    print(OUT/'DATA.md')
