"""Held-out wall geometry and 2D structure tests without fixing alpha=1/3.

No new PDE solution is computed. All r values are measured, except the one
explicitly labelled free-r validation. All model parameters and failed fits
are retained. Positivity certificates use mpmath interval arithmetic.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import mpmath as mp
from numpy.polynomial import Polynomial as Poly
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.integrate import cumulative_simpson
from fixed_point import fd_matrix

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'runs/steady_20260928_tau20'
OUT = ROOT/'runs/structure_free_alpha_20260928'
NAMES = ('L4_129x65', 'L4_257x129', 'L6_193x97', 'L8_257x129_continued')

z, r, c, eps, delta = sp.symbols('z r c eps delta', real=True)
w = 1+z*z
R = (1+z**4)**sp.Rational(1, 4)
chi = z/sp.sqrt(w)
SYMBOLIC = {
    'quartic': (1+c*R**r*(1+eps*z/R),
                 r+c*R**(r-4)*(r+eps*(r-1)*z/R)),
    'counterterm': (
        1+delta*chi+c*(w**(r/2)-r/2*w**(r/2-1)
            +eps*z*w**((r-1)/2)-eps*(r-1)/2*z*w**((r-3)/2)),
        r+delta*chi*(r-1/w)-c/2*w**(r/2-2)
            *(r*(r-2)+eps*(r-1)*(r-3)*chi))}
FUNCTIONS = {}
for name, (bb, hh) in SYMBOLIC.items():
    f = bb/hh**sp.Rational(1, 3)
    FUNCTIONS[name] = [sp.lambdify((z, r, c, eps, delta), expr, 'numpy')
                      for expr in (f, sp.diff(f, z), sp.diff(f, z, 2), hh)]


def normalize(family, rr, cc, ee, dd=0.):
    f, f1, f2, hh = FUNCTIONS[family]
    args = (rr, cc, ee, dd)
    side = 1 if f1(0., *args) <= 0 else -1
    hi = float(side)
    while side*f1(hi, *args) < 0 and abs(hi) < 1e4:
        hi *= 2
    zp = brentq(lambda zz: f1(zz, *args), min(0., hi), max(0., hi))
    fp = f(zp, *args)
    length = np.sqrt(f2(zp, *args)/(2*fp))
    hpeak = hh(zp, *args)
    return dict(family=family, r=float(rr), alpha=float(1/(rr+1)),
        c=float(cc), eps=float(ee), delta=float(dd), zpeak=float(zp),
        fpeak=float(fp), L=float(length), d=float(length*zp),
        muKpeak=float(fp/(2*hpeak**(2/3))),
        Jss_at_stagnation=float(f2(0., *args)/(length**2*fp)))


def positive_c_ceiling(rr, ee, dd):
    """Reduce H>0 on the entire axis to stationary roots on t in [-1,1].

    Roots are computed numerically for parameterization. Final parameter
    values receive a separate outward-rounded interval certificate.
    """
    aa, bb = rr*(rr-2), ee*(rr-1)*(rr-3)
    h = Poly([rr, dd*(rr-1), 0, dd])
    lin, one_minus_t2 = Poly([aa, bb]), Poly([1, 0, -1])
    stationary = (h.deriv()*one_minus_t2*lin
                  -h*(Poly([bb])*one_minus_t2-(4-rr)*Poly([0, 1])*lin))
    candidates = [float(t.real) for t in stationary.roots()
                  if abs(t.imag) < 1e-7 and abs(t.real) < 1
                  and aa+bb*t.real > 0]
    if not candidates:
        return None
    return min(2*h(t)/((1-t*t)**((4-rr)/2)*(aa+bb*t)) for t in candidates)


def unpack(family, measured_r, u):
    rr = u[3] if len(u) == 4 else measured_r
    ee = -1/(1+np.exp(-u[1]))
    dd = .98*np.tanh(u[2]) if len(u) > 2 else 0.
    ceiling = positive_c_ceiling(rr, ee, dd) if family == 'counterterm' else None
    cc = (ceiling/(1+np.exp(-u[0]))) if ceiling is not None else np.exp(u[0])
    result = normalize(family, rr, cc, ee, dd)
    result['c_positive_ceiling'] = None if ceiling is None else float(ceiling)
    return result


def certify_H(model, partitions=512):
    if model['family'] == 'quartic':
        return dict(method='analytic positive summands',
                    H_lower_bound=model['r'], certified=True)
    iv = mp.iv
    iv.dps = 40
    rr, cc, ee, dd = [iv.mpf(str(model[key])) for key in ('r', 'c', 'eps', 'delta')]
    lower = float('inf')
    for j in range(partitions):
        # Power-of-two mesh endpoints are exactly representable dyadics.
        t = iv.mpf([str(-1+2*j/partitions), str(-1+2*(j+1)/partitions)])
        hh = (rr+dd*t*(rr-1+t**2)-cc/2*(1-t**2)**((4-rr)/2)
              *(rr*(rr-2)+ee*(rr-1)*(rr-3)*t))
        lower = min(lower, float(np.nextafter(float(hh.a), -np.inf)))
    return dict(method='mpmath.iv on full compactified axis t in [-1,1]',
        partitions=partitions, decimal_precision=40,
        H_lower_bound=lower, certified=bool(lower > 0),
        interpretation='Certificate applies to the decimal parameters stored in this JSON; B positivity follows from H>0 and positive two-sided leading tails.')


def model_at_Q(model, qq, mu):
    args = tuple(model[key] for key in ('r', 'c', 'eps', 'delta'))
    zz = (qq+model['d'])/model['L']
    f, _, _, hh = FUNCTIONS[model['family']]
    jj = f(zz, *args)/model['fpeak']
    kk = model['fpeak']/(2*mu)*hh(zz, *args)**(-2/3)
    return 1/jj, kk, jj


def comparison(model, source):
    qq = np.linspace(-1.1, 1.1, 2001)
    predq, predk, jj = model_at_Q(model, qq, source['mu'])
    trueq, truek = source['q_of_Q'](qq), source['K_of_Q'](qq)
    train = abs(qq) <= .6
    xp = cumulative_simpson(jj, x=qq, initial=0)
    xp -= xp[len(qq)//2]
    tx = source['x_of_Q'](qq)
    common = (tx >= xp[0]) & (tx <= xp[-1]) & (abs(tx) <= 2)
    same_x_pred = PchipInterpolator(xp, predq)(tx[common])
    same_x_error = abs(same_x_pred-trueq[common])
    same_x_train = train[common]
    return dict(train_Q=[-.6, .6], comparison_Q=[-1.1, 1.1],
        train_q_linf=float(np.max(abs(predq[train]-trueq[train]))),
        holdout_q_linf=float(np.max(abs(predq[~train]-trueq[~train]))),
        train_K_linf_over_peak=float(np.max(abs(predk[train]-truek[train]))/source['Kpeak']),
        holdout_K_linf_over_peak=float(np.max(abs(predk[~train]-truek[~train]))/source['Kpeak']),
        same_x_q_linf=float(np.max(same_x_error)),
        same_x_train_q_linf=float(np.max(same_x_error[same_x_train])),
        same_x_holdout_q_linf=float(np.max(same_x_error[~same_x_train])),
        same_x_split='Training/holdout labels follow the true data density Q(x), using |Q|<=0.6.',
        same_x_interval=[float(tx[common][0]), float(tx[common][-1])],
        inverse_x_linf=float(np.max(abs(xp-tx))),
        d_relative_error=float(model['d']/source['d']-1),
        asymptotic_K_plus=float(model['fpeak']/(2*source['mu'])
            /(model['r']*(1+model['delta']))**(2/3)),
        asymptotic_K_minus=float(model['fpeak']/(2*source['mu'])
            /(model['r']*(1-model['delta']))**(2/3)))


def fit_model(source, family, delta_free=False, r_free=False, scalars_only=False):
    qq = np.linspace(-.6, .6, 25)
    qt, kt = source['q_of_Q'](qq), source['K_of_Q'](qq)
    def objective(u):
        model = unpack(family, source['r'], u)
        derror = model['d']/source['d']-1
        if scalars_only:
            return [derror, model['muKpeak']/(source['mu']*source['Kpeak'])-1]
        pq, pk, _ = model_at_Q(model, qq, source['mu'])
        return np.r_[derror, (pq-qt)/np.sqrt(len(qq)),
                     (pk-kt)/(source['Kpeak']*np.sqrt(len(qq)))]
    lower, upper = [-8, -8], [8, 10]
    if delta_free:
        lower += [-2.5]; upper += [2.5]
    if r_free:
        lower += [1.1]; upper += [2.9]
    solutions = []
    for initial in ([0, 3, -.5], [1, 4, -.2], [0, 2, 0]):
        u = initial if delta_free else initial[:2]
        if r_free:
            u = u+[source['r']]
        opt = least_squares(objective, u, bounds=(lower, upper), max_nfev=150,
                            xtol=1e-10, ftol=1e-10, gtol=1e-10)
        solutions.append(opt)
    opt = min(solutions, key=lambda item: np.linalg.norm(item.fun))
    model = unpack(family, source['r'], opt.x)
    return dict(parameters=model, comparison=comparison(model, source),
        positivity=certify_H(model),
        fit=dict(objective='d and peak inclination only' if scalars_only else
                 'd plus 25 q and K samples on -0.6<=Q<=0.6',
                 r_free=r_free, delta_free=delta_free, optimizer_status=int(opt.status),
                 message=opt.message, active_mask=opt.active_mask.tolist(),
                 objective_rms=float(np.sqrt(np.mean(opt.fun**2))),
                 parameters_transformed=opt.x.tolist(),
                 bounds=[lower, upper], scalar_defects=opt.fun.tolist() if scalars_only else None))


def source_and_features(name):
    data = np.load(SOURCE/name/'state.npz')
    x, y, Q = data['x'], data['y'], data['Q']
    lam, mu, v, _ = data['rates']
    dx, dy = fd_matrix(x, 1, 7), fd_matrix(y, 1, 7)
    q, qy = Q@dx.T, dy@Q
    K, J = -qy/q, 1/q
    ss = Q+v/mu
    source = dict(name=name, r=float(lam/mu-1), alpha=float(mu/lam),
        mu=float(mu), d=float(v/mu), Kpeak=float(K[0, len(x)//2]),
        q_of_Q=PchipInterpolator(Q[0], q[0]),
        K_of_Q=PchipInterpolator(Q[0], K[0]), x_of_Q=PchipInterpolator(Q[0], x))
    Js = (dx@J[0])/q[0]
    Jss = (dx@Js)/q[0]
    jets = dict(J=float(CubicSpline(ss[0], J[0])(0)),
                J_s=float(CubicSpline(ss[0], Js)(0)),
                J_ss=float(CubicSpline(ss[0], Jss)(0)))
    jets['J_ss_polynomial_windows'] = []
    for width in (.08, .15, .25):
        mask = abs(ss[0]) < width
        coef = np.polynomial.polynomial.polyfit(ss[0, mask], J[0, mask], min(4, sum(mask)-1))
        jets['J_ss_polynomial_windows'].append(dict(half_width=width,
            nodes=int(sum(mask)), value=float(2*coef[2])))

    # The B,H identity is a wall identity. Test rather than assume it off-wall.
    BB = 2*mu*J/np.sqrt(K)
    HH = K**(-1.5)
    Bs = (BB@dx.T)/q
    generator = []
    for height in (0., .5, 1., 2.):
        iy = int(np.argmin(abs(y-height)))
        mask = (Q[iy] >= -.6) & (Q[iy] <= .5)
        resid = source['r']*BB[iy]-ss[iy]*Bs[iy]-HH[iy]
        effective = np.dot(BB[iy, mask], HH[iy, mask]+ss[iy, mask]*Bs[iy, mask])/np.dot(BB[iy, mask], BB[iy, mask])
        generator.append(dict(y=float(y[iy]), effective_r=float(effective),
            identity_l2_over_H=float(np.linalg.norm(resid[mask])/np.linalg.norm(HH[iy, mask]))))

    # A derivative-level algebraic normal-width law, with its exact wall jet.
    levels = np.linspace(-.6, .5, 81)
    heights = y[y <= 2]
    _, yy = np.meshgrid(levels, heights)
    jtable = np.array([PchipInterpolator(Q[i], J[i])(levels) for i in range(len(heights))])
    Ks = (K@dx.T)/q
    kprime = PchipInterpolator(Q[0], Ks[0])(levels)
    jzero = jtable[0]
    def width_model(m):
        return jzero*(1+m*kprime*yy/jzero)**(1/m)
    training = yy <= 1
    best = minimize_scalar(lambda m: np.mean((width_model(m)[training]/jtable[training]-1)**2),
                           bounds=(.1, 5), method='bounded')
    width = []
    for label, m in [('square_root', 2.), ('m_equals_r', source['r']), ('fitted', float(best.x))]:
        error = width_model(m)/jtable-1
        width.append(dict(label=label, m=m,
            train_relative_rms=float(np.sqrt(np.mean(error[training]**2))),
            holdout_relative_rms=float(np.sqrt(np.mean(error[~training]**2))),
            holdout_relative_linf=float(np.max(abs(error[~training])))))

    tail = []
    for side in (-1, 1):
        pos = side*np.array([.5, .875])*max(x)
        values = CubicSpline(x, K[0])(pos)
        tail.append(dict(side=side, x=pos.tolist(), K=values.tolist(),
                         outward_relative_change=float(values[1]/values[0]-1)))
    return source, dict(name=name, measured_r=source['r'], measured_alpha=source['alpha'],
                        density_shift=source['d'], stagnation_jets=jets,
                        generator_by_height=generator, derivative_width_laws=width,
                        finite_box_tail_inclinations=tail)


def write_report(rows):
    lines = []
    for row in rows:
        mod = row['counterterm_delta_fit']
        p, e = mod['parameters'], mod['comparison']
        lines.append(f"| {row['name']} | {p['r']:.6f} | {p['c']:.5f} | {p['eps']:.5f} | {p['delta']:.5f} | {e['train_q_linf']:.5f} | {e['holdout_q_linf']:.5f} | {e['holdout_K_linf_over_peak']:.4f} |")
    free = rows[1]['counterterm_free_r_fit']
    text = r'''# 自由指数的数据结构检验

本次没有固定 α=1/3，也没有运行新的 PDE 求解。固定指数检验分别使用四组实测 `r=λ/μ−1`；另外仅对唯一细网格数据做了一次自由 r 检验。下面的低参数族都只是壁面和第一法向 jet 的候选，完整二维系统仍未闭合。

## 1. 先用几何特征筛选，而非只匹配峰点

记 `s=Q+v/μ`、`J=1/Qx`、`K=−Qy/Qx`。四次软根族虽能把 d 与峰倾角两个标量匹配到机器精度，却在驻点强制 `Jss(0)<0`。L4 细网格给 `Jss(0)≈5.31`，密度窗口 0.08、0.15、0.25 内的四次拟合分别给 5.54、6.50、6.84，稳定为正。该族因此在一个未参与拟合的局部几何特征上已经失败。

两个起始族还都强制左右无穷处 K 相同。L8 数据中右侧 x=3 到 7 时 K 仅从约 2.60 变到 2.71；左侧 x=−3 到 −7 则从约 0.725 降到 0.585。右尾较像平台，左尾仍未稳定；有限框不能证明真正尾极限，但不支持简单地把两个尾倾角设为同一个常数。

## 2. 允许不同倾角尾的自由指数族

使用主分析给出的 `w=1+z², χ=z/sqrt(w), z=s/L`：

\[
B=1+\delta\chi+c\left[w^{r/2}-\frac r2w^{r/2-1}
+\varepsilon z w^{(r-1)/2}-\frac{\varepsilon(r-1)}2z w^{(r-3)/2}\right],
\]
\[
H=r+\delta\chi\left(r-\frac1w\right)
-\frac c2w^{r/2-2}\left[r(r-2)+\varepsilon(r-1)(r-3)\chi\right].
\]

它满足 `rB−sB'=H`，并给出 `J=B/(2μA H^(1/3))`、`K=A H^(−2/3)`。峰规范确定 L、A，剩余结构参数为 c、ε、δ。两端 H 分别趋于 `r(1±δ)`，所以允许不等倾角尾。训练仅使用 d 和 `−0.6≤Q≤0.6` 内 25 个 q、K 值；留出 `0.6<|Q|≤1.1`。没有把留出曲线用于选参数。

| 数据 | 固定的实测 r | c | ε | δ | 训练 q 最大误差 | 留出 q 最大误差 | 留出 K 误差/峰 K |
|---|---:|---:|---:|---:|---:|---:|---:|
'''+'\n'.join(lines)+r'''

本表全部误差在同一个 Q 值比较。q 误差以规范化峰高 1 为单位，K 误差以各自峰倾角为单位，均不是逐点相对误差。唯一细网格的留出 q 误差约 0.61%，K 误差约 1.26%；这是这次最值得继续的壁面结构线索。将 J 积分成模型 x(Q) 后再作同 x 比较，固定 r 模型在公共区间约 [−1.814,1.419] 上的最大 q 误差约 0.486%，自由 r 模型约 0.483%；JSON 另按真实密度 |Q(x)|≤0.6 分开记录训练/留出误差。其他三个粗网格的同 Q 留出 q 误差仍为约 8%–10%，不能宣称此族统一命中了四个候选。细网格拟合还预测驻点 Jss 约 8.1，比直接离散计算约 5.3 大，说明高阶几何仍未匹配。

![细网格同密度的 q 与 K 留出检验](fine_heldout_geometry.png)

作为对照，限制 δ=0 时，细网格的留出 q 误差约 9.16%、K 误差约 27.3%，并且 ε 接近参数边界。改善来自允许不同倾角尾这一明确结构改变，而非把 α 固定成特殊分数。

## 3. 自由 r 的独立性范围

仅对细网格，进一步放开 `1.1<r<2.9`，保持同一训练/留出划分。得到

'''+f"- 拟合 r = {free['parameters']['r']:.9f}，α = {free['parameters']['alpha']:.9f}。\n- 实测速率 r = {rows[1]['measured_r']:.9f}，α = {rows[1]['measured_alpha']:.9f}。\n- 留出 q 最大误差 {free['comparison']['holdout_q_linf']:.6f}；K 最大误差/峰 K 为 {free['comparison']['holdout_K_linf_over_peak']:.6f}。\n"+r'''

几何拟合选出的 r 与实测速率相容，并未选择恰好 2。这仍使用了实测 μ 和 d，并非完全不依赖速率的指数识别，也没有证明连续极限或唯一性。

## 4. 全轴正性不是靠画图判断

粗充分界对最有希望的 δ 模型并不成立。因此程序先将 H 写为紧坐标 `t=z/sqrt(1+z²)∈[−1,1]` 上的显式函数，再对最终模型使用 `mpmath.iv` 40 位区间算术、512 个区间覆盖整个 [−1,1]。JSON 记录每个模型的外向舍入下界与认证状态。细网格固定 r 的三参数模型得到 H 下界约 0.2646；这是对保存的十进制参数所作的全轴正性检查，不是普通离散采样。

H>0 加上 `rB−sB'=H` 与 B 两尾的正首项保证 B>0：任何正半轴零点只能向负值穿过，无法再回到已知的正尾；负半轴同理。因此该参数模型 J、K 在全轴均为正。区间认证并不保证完整二维延拓存在。

## 5. 真正的二维检验

从二维场直接构造 `B=2μJ/sqrt(K)`、`H=K^(−3/2)` 后，壁面上的 `rB−sBs−H` 相对 L² 误差约 0.03%–0.30%。相同恒等式离壁后失效：y=0.5 时约 35%–40%，y=1 时约 62%–69%，y=2 时约 103%–111%。所以该壁面闭合不能原样复制到每个高度。

还检验了一个只含一个指数、保留精确初始法向导数的二维展宽式：

\[
J(s,y)^m=J(s,0)^m+mJ(s,0)^{m-1}K_s(s,0)y.
\]

训练 `0≤y≤1` 后得到 m 约 2.04–2.18；但在 `1<y≤2` 的相对 RMS 误差仍约 3.9%–5.7%，最大误差约 9%–13%。强行设 m=r 也没有明显优势。这个接近平方根的规律只是可追踪的几何线索，不是已经通过 Poisson 方程的二维精确公式。

所有参数、优化边界状态、普通与区间正性结果、同密度与同 x 误差以及二维检验均在 `free_data_checks.json`，脚本为 `solver/steady_analysis/structure_free_data.py`。
'''
    (OUT/'FREE_DATA.md').write_text(text)


def plot_finest(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    source, _ = source_and_features('L4_257x129')
    row = next(item for item in rows if item['name'] == 'L4_257x129')
    qq = np.linspace(-1.1, 1.1, 1001)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.1), constrained_layout=True)
    truths = [source['q_of_Q'](qq), source['K_of_Q'](qq)]
    for ax, truth, label in zip(axes, truths, ['wall gradient q', 'isodensity inclination K']):
        ax.axvspan(-.6, .6, color='#dce6ee', alpha=.65, label='training density interval')
        ax.plot(qq, truth, color='#111111', lw=2.4, label='saved fine-grid candidate')
        ax.set(xlabel='density Q (comparison at equal Q)', ylabel=label, xlim=(-1.1, 1.1))
        ax.grid(alpha=.2)
    for key, color, style, label in [
        ('counterterm_equal_tail_fit', '#cf8440', '--', 'equal inclination tails'),
        ('counterterm_delta_fit', '#2a689c', '-', 'different tails, measured r'),
        ('counterterm_free_r_fit', '#56956d', ':', 'different tails, fitted r')]:
        pq, pk, _ = model_at_Q(row[key]['parameters'], qq, source['mu'])
        axes[0].plot(qq, pq, style, color=color, lw=1.8, label=label)
        axes[1].plot(qq, pk, style, color=color, lw=1.8)
    axes[0].set_title('q: held-out max error 0.00613 of peak height')
    axes[1].set_title('K: held-out max error 1.26% of peak K')
    axes[0].legend(fontsize=7.5, loc='lower center')
    fig.suptitle('Free-exponent wall structure: L4, 257 x 129 | white regions are held out')
    fig.savefig(OUT/'fine_heldout_geometry.png', dpi=175)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in NAMES:
        source, row = source_and_features(name)
        row['quartic_two_scalar_fit'] = fit_model(source, 'quartic', scalars_only=True)
        row['counterterm_two_scalar_fit'] = fit_model(source, 'counterterm', scalars_only=True)
        row['counterterm_equal_tail_fit'] = fit_model(source, 'counterterm')
        row['counterterm_delta_fit'] = fit_model(source, 'counterterm', delta_free=True)
        if name == 'L4_257x129':
            row['counterterm_free_r_fit'] = fit_model(source, 'counterterm', delta_free=True, r_free=True)
        rows.append(row)
        print(name, row['counterterm_delta_fit']['comparison'], flush=True)
    result = dict(scope='Structure identification on saved finite-box profiles, not a new PDE solve.', data=rows)
    (OUT/'free_data_checks.json').write_text(json.dumps(result, indent=2)+'\n')
    write_report(rows)
    plot_finest(rows)
    print(OUT/'FREE_DATA.md')


if __name__ == '__main__':
    main()
