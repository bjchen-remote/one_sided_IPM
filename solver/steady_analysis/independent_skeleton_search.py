"""Independent, reproducible positive-hodograph model discrimination.

Reads raw tau20 and conditional roots; does not import any previous fitted shape.
All model Jacobians are positive for every real density coordinate.  Coefficients
are free functions of height: 'parameters' means parameters per section, not a
claim of a finite-constant exact steady solution.  A small projected residual is
never labelled a symbolic identity.

PYTHONPATH=/private/tmp/ipm_steady_deps OPENBLAS_NUM_THREADS=1 \
  /Users/cbjin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  solver/steady_analysis/independent_skeleton_search.py
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import RectBivariateSpline, PchipInterpolator, CubicSpline
from scipy.optimize import least_squares
from numpy.polynomial.legendre import leggauss

from fixed_point import fd_matrix

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "runs/steady_20260928_tau20"
OUT = ROOT / "runs/independent_skeleton_20260928"
LEVELS = np.linspace(-.7, .5, 97)
HEIGHTS = np.linspace(0., 1.5, 13)
GLX, GLW = leggauss(32)


def profiles():
    """Normalizations copied from export_seed.m, independently of old fits."""
    raw = loadmat(SOURCE / "run_analysis.mat", simplify_cells=True)
    meta = raw["metadata"]
    a, peak, ell = [float(meta[k]) for k in ("a", "P", "ell")]
    center = float(RectBivariateSpline(raw["yFull"], raw["xFull"], raw["rho"])(0, a)[0, 0])
    yield dict(name="native_tau20", kind="finite_time_raw_export", path=str(SOURCE / "run_analysis.mat"),
        x=(raw["xFull"]-a)/ell, y=raw["yFull"]/ell,
        Q=(raw["rho"]-center)/(peak*ell), q=raw["rhoX"]/peak,
        u1=(raw["u1"]-meta["uAtPeak"])/(peak*ell), u2=raw["u2"]/(peak*ell),
        rates=np.asarray(meta["innerRates"], float))
    for name in ("L4_129x65", "L4_257x129", "L6_193x97", "L8_257x129_continued"):
        path = SOURCE / name / "state.npz"
        z = np.load(path)
        x, y, psi = z["x"], z["y"], z["psi"]
        # Match the saved root's five-point velocity convention exactly.
        u1 = -fd_matrix(y, 1, 5) @ psi
        u2 = psi @ fd_matrix(x, 1, 5).T
        u1 -= float(RectBivariateSpline(y, x, u1)(0, 0)[0, 0])
        yield dict(name=name, kind="conditional_frozen_exterior_root", path=str(path),
            x=x, y=y, Q=z["Q"], q=z["qx"], u1=u1, u2=u2, rates=z["rates"])


def inverse_data(p):
    fields = {k: RectBivariateSpline(p["y"], p["x"], p[k]) for k in ("Q", "q", "u1", "u2")}
    xx, qq = [], []
    for h in HEIGHTS:
        # A dense physical mesh prevents native nonuniform-x interpolation
        # from being confused with a model's density-coordinate fit error.
        dense_x = np.linspace(max(-6., p["x"].min()), min(6., p["x"].max()), 12001)
        row = fields["Q"](h, dense_x)[0]
        assert np.all(np.diff(row) > 0), (p["name"], h, np.min(np.diff(row)))
        inv = PchipInterpolator(row, dense_x, extrapolate=False)(LEVELS)
        assert np.all(np.isfinite(inv)), (p["name"], h)
        xx.append(inv)
        qq.append(fields["q"](h, inv)[0])
    return fields, np.array(xx), np.array(qq)


def model(name, pars, s=LEVELS):
    """Return X(s) and globally positive J=X_s; all primitives vanish at s=0."""
    s = np.asarray(s)
    d, logc = pars[:2]
    c = np.exp(logc)
    if name == "positive_cubic":
        a, b = np.exp(pars[2]), np.tanh(pars[3])
        j = c*(1+2*b*np.sqrt(a)*s+a*s*s)
        x = d+c*(s+b*np.sqrt(a)*s*s+a*s**3/3)
    elif name == "sos_quintic":
        b, a, e = pars[2:]
        j = c*(1+(b+a*s+e*s*s)**2)
        x = d+c*((1+b*b)*s+b*a*s*s+(a*a+2*b*e)*s**3/3+a*e*s**4/2+e*e*s**5/5)
    elif name == 'saturated_sos':
        b, a = pars[2], np.exp(pars[3])
        e = a*np.tanh(pars[4])
        rt = np.sqrt(1+s*s)
        j = c*(1+(b+a*s+e*s*s/rt)**2)
        x = d+c*((1+b*b)*s+a*b*s*s+a*a*s**3/3+b*e*(s*rt-np.arcsinh(s))+
                 2*a*e*((s*s-2)*rt+2)/3+e*e*(s**3/3-s+np.arctan(s)))
    elif name == "algebraic_root":
        a, center = np.exp(pars[2]), pars[3]
        z, z0 = s-center, -center
        j = c*(1+2*a*z*z)/np.sqrt(1+a*z*z)
        x = d+c*(z*np.sqrt(1+a*z*z)-z0*np.sqrt(1+a*z0*z0))
    elif name == "positive_rational":
        a, b, center = np.exp(pars[2]), np.exp(pars[3]), pars[4]
        z, z0, rb = s-center, -center, np.sqrt(b)
        j = c+a*z*z/(1+b*z*z)
        x = d+c*s+(a/b)*(s-(np.arctan(rb*z)-np.arctan(rb*z0))/rb)
    elif name == "positive_power":
        a, center, power = np.exp(pars[2]), pars[3], pars[4]
        j = c*(1+a*(s-center)**2)**power
        abscissa = (GLX[:, None]+1)*s[None, :]/2
        primitive = s/2*np.sum(GLW[:, None]*(1+a*(abscissa-center)**2)**power, axis=0)
        x = d+c*primitive
    elif name == "positive_exponential":
        a, center = np.exp(pars[2]), pars[3]
        j = c*np.exp(a*(s-center)**2)
        abscissa = (GLX[:, None]+1)*s[None, :]/2
        primitive = s/2*np.sum(GLW[:, None]*np.exp(a*(abscissa-center)**2), axis=0)
        x = d+c*primitive
    else:
        raise ValueError(name)
    assert np.all(j > 0)
    return x, j


SPECS = {
    "positive_cubic": dict(initial=[0, 0, .1, 0], bounds=([-10,-5,-8,-4],[10,5,5,4])),
    "sos_quintic": dict(initial=[0, 0, 0, 1, 0], bounds=([-10,-5,-4,-8,-8],[10,5,4,8,8])),
    'saturated_sos': dict(initial=[0, 0, 0, .1, .1], bounds=([-10,-5,-4,-5,-3],[10,5,4,4,3])),
    "algebraic_root": dict(initial=[0, 0, -.4, 0], bounds=([-10,-5,-8,-2],[10,5,5,2])),
    "positive_rational": dict(initial=[0, 0, 0, -1, 0], bounds=([-10,-5,-8,-8,-2],[10,5,5,5,2])),
    "positive_power": dict(initial=[0, 0, 0, 0, 1], bounds=([-10,-5,-12,-2,.05],[10,5,5,2,100])),
    "positive_exponential": dict(initial=[0, 0, -.2, 0], bounds=([-10,-5,-8,-2],[10,5,2,2])),
}


def error(value, truth, mask):
    err = np.asarray(value)[mask]-np.asarray(truth)[mask]
    return dict(rms=float(np.sqrt(np.mean(err*err))), linf=float(np.max(abs(err))),
                relative_l2=float(np.linalg.norm(err)/max(np.linalg.norm(np.asarray(truth)[mask]),1e-30)))


def projection(name, pars, h, fields, rates, train, measured_parameter_prime=None):
    """Actual saved Darcy velocity stays fixed, evaluated at candidate X(s).

    The inverse transport equation is V1 - V2 X_y - (mu*s+b) X_s=0.
    X_y = sum_i X_{p_i} p_i'.  Fit only p_i', never velocity or omitted terms.
    This is a numerical screen, not self-consistent Biot-Savart validation.
    """
    if h == 0 and measured_parameter_prime is None:
        return None
    x, j = model(name, pars)
    lam, mu, drift, offset = rates
    v1 = fields["u1"](h, x)[0]+lam*x+drift
    v2 = fields["u2"](h, x)[0]+lam*h
    columns = []
    for k in range(len(pars)):
        step = 2e-5*max(1., abs(pars[k]))
        plus, minus = pars.copy(), pars.copy()
        plus[k] += step; minus[k] -= step
        columns.append((model(name, plus)[0]-model(name, minus)[0])/(2*step))
    matrix = v2[:,None]*np.array(columns).T
    forcing = (mu*LEVELS+offset)*j-v1
    scale = np.maximum(np.linalg.norm(matrix[train], axis=0), 1e-30)
    coef, _, rank, _ = np.linalg.lstsq(matrix[train]/scale, -forcing[train], rcond=1e-11)
    prime = coef/scale
    residual = matrix @ prime+forcing
    scale_terms = np.sqrt(np.linalg.norm(v1[~train])**2+np.linalg.norm((v2*(np.array(columns).T@prime))[~train])**2+
                          np.linalg.norm(((mu*LEVELS+offset)*j)[~train])**2)
    # A direct measured baseline, evaluated where the inverse model places X.
    measured_q = fields["Q"](h, x, dy=1)[0]
    measured_y = fields["Q"](h, x, dx=1)[0]
    measured_Q = fields["Q"](h, x)[0]
    baseline = v1*measured_q+v2*measured_y-mu*measured_Q-offset
    answer = dict(rank=int(rank), prime=prime.tolist(), residual_holdout_linf=float(np.max(abs(residual[~train]))),
        residual_holdout_rms=float(np.sqrt(np.mean(residual[~train]**2))),
        residual_holdout_term_relative_l2=float(np.linalg.norm(residual[~train])/max(scale_terms,1e-30)),
        measured_transport_baseline_rms=float(np.sqrt(np.mean(baseline[~train]**2))),
        measured_inverse_transport_baseline_rms=float(np.sqrt(np.mean((baseline[~train]/measured_q[~train])**2))),
        singular_values=np.linalg.svd(matrix[train]/scale,compute_uv=False).tolist())
    if measured_parameter_prime is not None:
        true_k = np.array(columns).T @ measured_parameter_prime
        true_residual = v2*true_k+forcing
        true_scale = np.sqrt(np.linalg.norm(v1[~train])**2+np.linalg.norm((v2*true_k)[~train])**2+
                             np.linalg.norm(((mu*LEVELS+offset)*j)[~train])**2)
        answer['geometry_derivative_check'] = dict(
            method='CubicSpline across fitted per-height parameters, no residual minimization',
            parameter_prime=measured_parameter_prime.tolist(),
            inverse_residual_holdout_rms=float(np.sqrt(np.mean(true_residual[~train]**2))),
            forward_residual_holdout_rms=float(np.sqrt(np.mean((true_residual[~train]/j[~train])**2))),
            residual_holdout_term_relative_l2=float(np.linalg.norm(true_residual[~train])/max(true_scale,1e-30)),
            implied_K_difference_rms=float(np.sqrt(np.mean((true_k-(np.array(columns).T@prime))**2))))
    return answer


def write_report(result):
    lines = ['# 独立正反函数骨架：数据竞争与未闭合证据', '',
        '本次不采用旧报告选定的横向形状，仅复用原始数据导出约定和实际速度差分约定。'
        '原生 tau20 与每个条件有限箱稳态分别计分。尚未找到同时满足数据、正性和精确余空间消失的骨架。', '',
        '## 数据与留出', '',
        '- 原生：`runs/steady_20260928_tau20/run_analysis.mat`，按 export_seed.m 的 a、P、ell 归一化；它是有限时间数据。',
        '- 条件根：L4_129x65、L4_257x129、L6_193x97、L8_257x129_continued 的 state.npz；冻结外部数据，不能等同无限半平面解。',
        '- 反函数坐标 s=Q（未做密度平移），x=X(s,y)，q=1/J，J=X_s。97个密度点 s∈[-0.7,0.5]，13个高度 y∈[0,1.5]；无x外插。',
        '- interlaced：密度点交错训练/留出。central_to_wings：只用 s∈[-0.4,0.2]训练，检验两翼 [-0.7,-0.4)∪(0.2,0.5]。',
        '- 高度各自拟合。下列4或5个参数均是随高度变化的函数，包含X平移D；并非全局4或5个常数。没有把实测整条壁面当一个参数。', '',
        '## 完整候选库', '',
        '|族|截面参数数|正Jacobian J|', '|---|---:|---|',
        '|positive_cubic|4|C[1+2 tanh(b)√A s+A s²]|',
        '|sos_quintic|5|C[1+(b+a s+e s²)²]|',
        '|saturated_sos|5|C[1+(b+a s+e s²/√(1+s²))²]，a>0且abs(e)<a|',
        '|algebraic_root|4|C[1+2A(s-c)²]/√[1+A(s-c)²]|',
        '|positive_rational|5|C+A(s-c)²/[1+B(s-c)²]|',
        '|positive_power|5|C[1+A(s-c)²]^p，p∈[0.05,100]|',
        '|positive_exponential|4|C exp[A(s-c)²]|', '',
        'C,A,B均正。各J在整个实密度轴严格为正，X=D+∫₀ˢJ(t)dt单调且满射，所以q=1/J>0。'
        '这一几何正性不提供Darcy–输运闭合。', '',
        '## 分开计分的翼留出', '',
        '数值为各高度X留出RMSE的中位数；所有模型列在JSON中。', '',
        '|数据|五次SOS|饱和SOS|4参数指数|翼留出几何优者|', '|---|---:|---:|---:|---|']
    names = list(dict.fromkeys(r['profile'] for r in result['aggregate']))
    for profile in names:
        rows = {r['family']:r for r in result['aggregate'] if r['profile']==profile and r['split']=='central_to_wings'}
        vals = [rows[f]['median_X_holdout_rms'] for f in ('sos_quintic','saturated_sos','positive_exponential')]
        winner = min(rows,key=lambda f:rows[f]['median_X_holdout_rms'])
        lines.append(f'|{profile}|'+'|'.join(f'{v:.6g}' for v in vals)+f'|{winner}|')
    lines += ['', '原生和L4细网格偏好五次SOS；L4粗网格、L6、L8偏好饱和SOS。'
        '这只是本候选库内的局部几何竞争，不能视作数据确定了唯一精确族。'
        '五次SOS的e≠0时X尾为s⁵；饱和SOS的两翼均为s³。后者的尾指数仍是候选先验，不能由核心拟合认证。', '',
        '自由幂族在原生tau20和L4细网格的所有高度均趋p=100上界，与正指数Jacobian极限接近。'
        '因此本密度区间没有选出稳定的有限p，也没有据此识别远场尾指数。', '',
        '## 余空间诊断与几何导数自洽', '',
        '实际测得的Darcy速度保留为独立输入，在候选位置X(s,y)取值。令V₁=u₁+λX+v、V₂=u₂+λy，'
        '反函数输运式为 V₁−V₂X_y−(μs+b₀)J=0。仅对参数导数pᵢ′投影，矩阵列为V₂∂X/∂pᵢ；'
        '未修改速度，未删除输运项。这个固定测量速度的投影仍不等于用候选q重算Biot–Savart后的自洽残差。', '',
        '同时对各高度拟合参数做CubicSpline，以真实参数轨道导数替代自由投影导数。'
        '下面只用内部高度，比较同单位前向输运RMS；中位数按高度取。', '',
        '|数据/模型|沿几何轨道前向RMS|原始测量输运baseline RMS|', '|---|---:|---:|']
    for profile in ('native_tau20','L4_257x129'):
        for family in ('sos_quintic','saturated_sos'):
            rows = [r for r in result['rows'] if r['profile']==profile and r['family']==family and
                    r['split']=='central_to_wings' and r['geometry_derivative_interior_height']]
            geom = np.median([r['actual_velocity_projection']['geometry_derivative_check']['forward_residual_holdout_rms'] for r in rows])
            base = np.median([r['actual_velocity_projection']['measured_transport_baseline_rms'] for r in rows])
            lines.append(f'|{profile}/{family}|{geom:.7g}|{base:.7g}|')
    lines += ['', '原生tau20的非零稳态残差大部分已经存在于原数据中，不能全归因于骨架。'
        'L4细网格条件根baseline约10⁻⁷，而骨架轨道约10⁻³，清楚显示“每层自由投影很小”并不代表沿实际几何轨道闭合。'
        '投影小残差不是R_perp≡0的符号证明。'
        '高度步长从0.125减到0.0625后，L4细网格SOS翼留出前向RMS由0.000779784变为0.000779721，'
        '饱和SOS由0.001038609变为0.001038754；非零量不由高度采样粗度解释。', '',
        '## 更宽物理域的独立数据门', '',
        '所有模型保持原核心拟合参数，在共同物理域x∈[-4,4]重新与实测q比较，不重新拟合。'
        '壁面q相对L2误差：', '',
        '|数据|五次SOS|饱和SOS|正三次|', '|---|---:|---:|---:|',
        '|原生tau20|49.73%|7.23%|15.87%|',
        '|L4细网格|23.28%|3.74%|13.71%|', '',
        '五次SOS的平方二次式有两个实根时，会产生两个等高q峰。原生拟合的第二峰位于x=3.529944，'
        '预测q=1.001483，而实测q=0.317160；原始数据完全覆盖这一位置。'
        '因此五次SOS在核心的领先不能推广为全局骨架证据。饱和SOS在扩展域中更稳健，但仍未满足方程闭合。'
        '完整140条扩展域结果见broad_domain_screen.json。', '',
        '## 精确核验的范围', '',
        '脚本内symbolic_model_checks逐式验证五次SOS和饱和SOS的X_s=J，并给出全实轴正性表示。'
        '它只认证坐标几何。另一个并行agent负责Darcy–输运壁面及二维符号兼容；本报告不以清分母后“有限个系数”冒充余空间消失。'
        '五次SOS已有并行符号排除证据；饱和SOS的规范壁面C=1、b=0、a=1、非零e且abs(e)<1，'
        '连同该族自然参数切向，也被壁面必要条件排除。该结论不排除一般隐式代数骨架，'
        '更不能因为有限候选失败而宣称目标不存在。独立符号证明见'
        '`solver/steady_analysis/fresh_symbolic_closure.py`及`runs/fresh_symbolic_closure_20260928/symbolic_certificates.json`。', '',
        '## 复跑', '', '```sh',
        'PYTHONPATH=/private/tmp/ipm_steady_deps OPENBLAS_NUM_THREADS=1 /Users/cbjin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 solver/steady_analysis/independent_skeleton_search.py',
        '```', '',
        'data_screen.json保留910个逐截面记录、全部参数、优化器状态、训练/留出误差、投影导数、真实轨道导数及奇异值。'
        'data_screen.npz保存反函数数据与每族拟合曲线；symbolic_model_checks.json保存精确几何恒等式。']
    (OUT/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n')


def symbolic_model_checks():
    import sympy as sp
    s,b,a,e = sp.symbols('s b a e',real=True)
    rt = sp.sqrt(1+s*s)
    quintic = (1+b*b)*s+b*a*s*s+(a*a+2*b*e)*s**3/3+a*e*s**4/2+e*e*s**5/5
    saturated = (1+b*b)*s+a*b*s*s+a*a*s**3/3+b*e*(s*rt-sp.asinh(s)) + \
                2*a*e*((s*s-2)*rt+2)/3+e*e*(s**3/3-s+sp.atan(s))
    checks = {
        'quintic_derivative':str(sp.factor(sp.diff(quintic,s)-(1+(b+a*s+e*s*s)**2))),
        'saturated_derivative':str(sp.simplify(sp.diff(saturated,s)-(1+(b+a*s+e*s*s/rt)**2))),
        'quintic_X':str(quintic), 'saturated_X':str(saturated),
        'positive_q_scope':'For C>0 both Jacobians are C*(1+real_square)>0, hence q=1/J>0; this is not a PDE closure certificate.'}
    assert checks['quintic_derivative'] == checks['saturated_derivative'] == '0'
    (OUT/'symbolic_model_checks.json').write_text(json.dumps(checks,indent=2)+'\n')


def height_refinement(result):
    """Check derivative diagnostics using twice as many sampled heights."""
    global HEIGHTS
    old_heights = HEIGHTS.copy()
    HEIGHTS = np.linspace(old_heights[0],old_heights[-1],2*len(old_heights)-1)
    refined = []
    for p in profiles():
        if p['name'] not in ('native_tau20','L4_257x129'):
            continue
        fields, truth, _ = inverse_data(p)
        for family in ('sos_quintic','saturated_sos'):
            for split, mask in result['train_masks'].items():
                train = np.array(mask,dtype=bool)
                old = [r for r in result['rows'] if r['profile']==p['name'] and r['family']==family and r['split']==split]
                params = []
                for k,h in enumerate(HEIGHTS):
                    initial = np.array(min(old,key=lambda r:abs(r['y']-h))['parameters'])
                    fit = least_squares(lambda z:(model(family,z)[0]-truth[k])[train],initial,
                        bounds=SPECS[family]['bounds'],max_nfev=600,ftol=1e-12,xtol=1e-12,gtol=1e-12)
                    if family == 'sos_quintic' and fit.x[3] < 0:
                        fit.x[2:] *= -1
                    params.append(fit.x)
                curve = CubicSpline(HEIGHTS,np.array(params),axis=0)
                before, after = [], []
                for r in old[1:-1]:
                    before.append(r['actual_velocity_projection']['geometry_derivative_check']['forward_residual_holdout_rms'])
                    proj = projection(family,curve(r['y']),r['y'],fields,p['rates'],train,curve(r['y'],1))
                    after.append(proj['geometry_derivative_check']['forward_residual_holdout_rms'])
                refined.append(dict(profile=p['name'],family=family,split=split,
                    coarse_height_spacing=float(old_heights[1]-old_heights[0]),
                    fine_height_spacing=float(HEIGHTS[1]-HEIGHTS[0]),
                    coarse_median_forward_residual_rms=float(np.median(before)),
                    fine_median_forward_residual_rms=float(np.median(after)),
                    max_change_in_forward_residual_rms=float(np.max(abs(np.array(before)-np.array(after))))))
    HEIGHTS = old_heights
    (OUT/'height_refinement.json').write_text(json.dumps(refined,indent=2)+'\n')


def broad_domain_screen(result):
    """Audit actual physical data outside the fitted density interval."""
    rows = []
    xs = np.linspace(-4,4,801)
    ss = np.linspace(-8,8,24001)
    for p in profiles():
        qfield = RectBivariateSpline(p['y'],p['x'],p['q'])
        for family in SPECS:
            for h in (0.,.5,1.,1.5):
                row = next(r for r in result['rows'] if r['profile']==p['name'] and
                           r['family']==family and r['y']==h and r['split']=='interlaced')
                pars = np.array(row['parameters'])
                xx, _ = model(family,pars,ss)
                assert np.all(np.diff(xx)>0)
                density = PchipInterpolator(xx,ss,extrapolate=False)(xs)
                assert np.all(np.isfinite(density))
                pred = 1/model(family,pars,density)[1]
                actual = qfield(h,xs)[0]
                info = dict(profile=p['name'],family=family,y=h,physical_domain=[-4,4],
                            q_error=error(pred,actual,np.ones(len(xs),dtype=bool)),
                            q_outer_error=error(pred,actual,abs(xs)>=2))
                if family == 'sos_quintic':
                    roots = np.roots(pars[[4,3,2]])
                    if np.all(np.isreal(roots)):
                        roots=np.sort(roots.real)
                        maxima,j=model(family,pars,roots)
                        measured = qfield.ev(np.full(len(maxima),h),maxima)
                        inside=(maxima>=p['x'].min())&(maxima<=p['x'].max())
                        info['equal_height_q_maxima'] = [dict(Q=float(q),x=float(x),model_q=float(1/jj),
                            measured_q=float(qm) if valid else None,inside_source_domain=bool(valid))
                            for q,x,jj,qm,valid in zip(roots,maxima,j,measured,inside)]
                rows.append(info)
    value=dict(scope='Models trained only on density [-.7,.5], now tested on actual physical x in [-4,4]. No refit.',rows=rows)
    (OUT/'broad_domain_screen.json').write_text(json.dumps(value,indent=2)+'\n')


def main():
    OUT.mkdir(exist_ok=True)
    allrows, curve_arrays = [], {}
    # Alternate-point validation tests interpolation. Contiguous density
    # holdout tests extrapolation, with both left and right wings excluded.
    splits = {"interlaced": np.arange(len(LEVELS))%2==0,
              "central_to_wings": (LEVELS>=-.4)&(LEVELS<=.2)}
    start = time.time()
    for p in profiles():
        fields, truth, qtruth = inverse_data(p)
        curve_arrays[p["name"]+"__truth_X"] = truth
        curve_arrays[p["name"]+"__truth_q"] = qtruth
        for split, train in splits.items():
            for name, spec in SPECS.items():
                previous = None
                fitted, params = [], []
                group_rows = []
                for k, h in enumerate(HEIGHTS):
                    initial = np.array(spec["initial"], float)
                    initial[0] = np.interp(0, LEVELS, truth[k])
                    initial[1] = np.log(1/np.interp(0, LEVELS, qtruth[k]))
                    starts = [initial]
                    if previous is not None:
                        starts.append(previous)
                    # Sum-of-squares has a nonconvex polynomial factorization.
                    if name == "sos_quintic":
                        starts += [initial+np.array([0, -.1, .2, 0, .2]), initial+np.array([0, -.1, -.2, 0, -.2])]
                    fits = [least_squares(lambda z: (model(name,z)[0]-truth[k])[train], z,
                                          bounds=spec["bounds"], max_nfev=1200,
                                          ftol=1e-11, xtol=1e-11, gtol=1e-11) for z in starts]
                    fit = min(fits,key=lambda f: np.linalg.norm(f.fun))
                    if not fit.success:
                        fit = least_squares(lambda z: (model(name,z)[0]-truth[k])[train], fit.x,
                            bounds=spec['bounds'], max_nfev=6000,ftol=1e-11,xtol=1e-11,gtol=1e-11)
                    if name == 'sos_quintic' and fit.x[3] < 0:
                        fit.x[2:] *= -1
                    previous = fit.x
                    xhat, jhat = model(name,fit.x)
                    fitted.append(xhat); params.append(fit.x)
                    group_rows.append(dict(profile=p["name"], kind=p["kind"], source=p["path"],
                        rates=p["rates"].tolist(), split=split, family=name, n_parameters=len(fit.x),
                        y=float(h), parameters=fit.x.tolist(), optimizer_success=bool(fit.success),
                        X_train=error(xhat,truth[k],train), X_holdout=error(xhat,truth[k],~train),
                        q_holdout=error(1/jhat,qtruth[k],~train),
                        actual_velocity_projection=projection(name,fit.x,float(h),fields,p["rates"],train)))
                param_curve = CubicSpline(HEIGHTS,np.array(params),axis=0)
                for row in group_rows:
                    row['actual_velocity_projection'] = projection(name,np.array(row['parameters']),row['y'],fields,p['rates'],train,
                                                                   param_curve(row['y'],1))
                    row['geometry_derivative_interior_height'] = bool(0 < row['y'] < HEIGHTS[-1])
                allrows.extend(group_rows)
                curve_arrays[p["name"]+"__"+split+"__"+name+"__X"] = np.array(fitted)
                curve_arrays[p["name"]+"__"+split+"__"+name+"__parameters"] = np.array(params)
                print(p["name"],split,name,"done",flush=True)
    aggregate = []
    profile_names = list(dict.fromkeys(r['profile'] for r in allrows))
    for split, profile_name in [(split,pn) for split in splits for pn in profile_names]:
        for name in SPECS:
            rows = [r for r in allrows if r["split"]==split and r["family"]==name and r['profile']==profile_name]
            projections = [r["actual_velocity_projection"] for r in rows if r["y"]>0]
            interior = [r["actual_velocity_projection"]['geometry_derivative_check'] for r in rows if r['geometry_derivative_interior_height']]
            aggregate.append(dict(split=split, profile=profile_name, kind=rows[0]['kind'], family=name,n_parameters=rows[0]["n_parameters"],
                median_X_holdout_rms=float(np.median([r["X_holdout"]["rms"] for r in rows])),
                max_X_holdout_linf=float(max(r["X_holdout"]["linf"] for r in rows)),
                median_q_holdout_relative_l2=float(np.median([r["q_holdout"]["relative_l2"] for r in rows])),
                median_projected_residual_relative_l2=float(np.median([r["residual_holdout_term_relative_l2"] for r in projections])),
                max_projected_residual_linf=float(max(r["residual_holdout_linf"] for r in projections)),
                median_geometry_derivative_residual_relative_l2=float(np.median([r['residual_holdout_term_relative_l2'] for r in interior])),
                median_geometry_derivative_forward_residual_rms=float(np.median([r['forward_residual_holdout_rms'] for r in interior]))))
    result = dict(scope="Independent data screening, positive coordinate maps, fixed measured Darcy-velocity diagnostic; no exact 2D closure is asserted.",
        positivity={"positive_cubic":"J=C[(sqrt(A)s+tanh(b))^2+sech(b)^2]>0",
            "sos_quintic":"J=C[1+(b+a*s+e*s^2)^2]>0",
            'saturated_sos':'J=C[1+(b+a*s+e*s^2/sqrt(1+s^2))^2]>0, a>0, |e|<a; cubic inverse tails on both sides',
            "algebraic_root":"J=C(1+2*A*(s-c)^2)/sqrt(1+A*(s-c)^2)>0",
            "positive_rational":"J=C+A*(s-c)^2/(1+B*(s-c)^2)>0",
            "positive_power":"J=C[1+A*(s-c)^2]^p>0",
            "positive_exponential":"J=C*exp(A*(s-c)^2)>0"},
        parameter_count="Section parameters are functions of height, not a finite set of global constants.",
        coordinate="s=Q (unshifted), x=X(s,y), q=Q_x=1/X_s",
        density_domain=LEVELS[[0,-1]].tolist(), heights=HEIGHTS.tolist(),
        train_masks={k:v.tolist() for k,v in splits.items()},
        aggregate=aggregate, rows=allrows, elapsed_seconds=time.time()-start)
    (OUT/"data_screen.json").write_text(json.dumps(result,indent=2)+"\n")
    np.savez_compressed(OUT/"data_screen.npz",levels=LEVELS,heights=HEIGHTS,**curve_arrays)
    symbolic_model_checks()
    height_refinement(result)
    broad_domain_screen(result)
    write_report(result)
    print(json.dumps(aggregate,indent=2),flush=True)


if __name__ == "__main__":
    main()
