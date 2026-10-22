"""Tail-compatible analytic masters: model discrimination, not PDE solutions.

Fits only the existing tau20 master on |xi| <= 1. Outer xi intervals and
individual rows are diagnostics. All geometry functions remain measured.
"""
from pathlib import Path
import json

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/geometry_loop_20260928"
BETA = 1.5


def master(xi, alpha, b, epsilon):
    p = BETA * (1 - alpha) / b
    z = np.sinh(b * xi) / np.sinh(b)
    odd = np.tanh(b * xi) * (z*z - 1) / (z*z + 1)
    M = np.expm1(2*np.log(2)/p)
    return np.exp(-p/2 * np.log1p(M*z*z*np.exp(epsilon*odd)))


def metrics(pred, obs):
    err = np.asarray(pred) - obs
    return dict(relative_l2=float(np.linalg.norm(err)/np.linalg.norm(obs)),
                max_absolute=float(np.max(abs(err))))


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    saved = loadmat(ROOT / "runs/construction_20260928/geometry_transform.mat",
                    simplify_cells=True)
    native = loadmat(ROOT / "runs/steady_20260928_tau20/run_analysis.mat",
                     simplify_cells=True)
    xx, yy, qq = (native[k] for k in ("xFull", "yFull", "rhoX"))
    field = RectBivariateSpline(yy, xx, qq)
    xi = saved['xi']; target = saved['learnedMaster']; train = abs(xi) <= 1
    outer_xi = np.linspace(-3, 3, 601)
    heights = np.r_[0., np.geomspace(.025, 1024., 60)]
    rows, valid, geometries = [], [], []
    for y in heights:
        row = field(y, xx)[0]
        s = CubicSpline(xx, row)
        roots = s.derivative().roots(extrapolate=False)
        roots = roots[(roots > xx[0]) & (roots < xx[-1])]
        m = float(roots[np.argmax(s(roots))]); A = float(s(m))
        half = CubicSpline(xx, row-A/2).roots(extrapolate=False)
        L = m - half[half < m][-1]; R = half[half > m][0] - m
        W = (L+R)*BETA/(2*np.sinh(BETA))
        k = (R-L)/(R+L)/np.tanh(BETA/2)
        if abs(k) >= 1:
            raise ValueError('Noninvertible sinh geometry')
        phys = m + W*(np.sinh(BETA*outer_xi)
                      + k*(np.cosh(BETA*outer_xi)-1))/BETA
        good = (phys >= xx[0]) & (phys <= xx[-1])
        rowout = np.full_like(phys, np.nan)
        rowout[good] = field.ev(np.full(np.sum(good), y), phys[good])/A
        rows.append(rowout); valid.append(good)
        geometries.append([m, A, W, k])
    rows, valid = np.asarray(rows), np.asarray(valid)
    fits = []
    for alpha in [.2, 1/3, .3654054862557143, .4, .5]:
        fit = least_squares(lambda p: (master(xi, alpha, *p)-target)[train],
                            [.4, -.3], bounds=([.02, -3], [3, 3]),
                            ftol=1e-12, xtol=1e-12, gtol=1e-12)
        b, eps = fit.x
        pred = master(xi, alpha, b, eps)
        ext = np.broadcast_to(master(outer_xi, alpha, b, eps), rows.shape)
        tests = []
        for lo, hi in [(0, 1), (1, 64), (64, 1024)]:
            ymask = (heights >= lo) & (heights <= hi)
            for xlo, xhi in [(0, 1), (1, 1.5), (1.5, 3)]:
                mask = (ymask[:, None] & (abs(outer_xi)[None, :] > xlo)
                        & (abs(outer_xi)[None, :] <= xhi) & valid)
                tests.append(dict(y_range=[lo, hi], abs_xi_range=[xlo, xhi],
                                  count=int(mask.sum()), **metrics(ext[mask], rows[mask])))
        fits.append(dict(alpha=alpha, b=float(b), epsilon=float(eps),
                         p=float(BETA*(1-alpha)/b),
                         master_train=metrics(pred[train], target[train]),
                         master_holdout=metrics(pred[~train], target[~train]),
                         native_row_audits=tests))
    # Exact certificate for Sol's prospective two-scalar tail coordinate.
    # h'=xi^2(xi^2-1)(11xi^2-3)/(1+xi^2)^(9/2), so this
    # positive polynomial certifies |h'|<1/2 on the entire real axis.
    import sympy as sp
    t = sp.symbols('t', nonnegative=True)
    certificate = sp.Poly((1+t)**9 - 4*t**2*(t-1)**2*(11*t-3)**2, t)
    root_count = certificate.count_roots(0, sp.oo)
    assert root_count == 0 and certificate.eval(0) == 1
    report = dict(
        source='tau20 geometry_transform.mat master + native run_analysis.mat',
        formula='f=[1+(2^(2/p)-1) z^2 exp(epsilon T)]^(-p/2)',
        definitions={'z': 'sinh(b xi)/sinh(b)',
                     'T': 'tanh(b xi)*(z^2-1)/(z^2+1)',
                     'p': '1.5*(1-alpha)/b'},
        exact_properties=['For b>0, 0<alpha<1 and finite epsilon: 0<f<=1.',
                          'f(0)=1, f_prime(0)=0, f(-1)=f(1)=1/2.',
                          'f ~ C_pm exp[-1.5*(1-alpha)*abs(xi)].',
                          'With |k|<1 sinh geometry, q ~ C_pm(y)*abs(x)^(alpha-1).'],
        nonidentifiability='At the current core measurement precision, widely different alpha fit almost equally well; no tail exponent is selected.',
        fits=fits, masked_out_of_domain_points=int(np.size(valid)-np.sum(valid)),
        tail_coordinate_certificate=dict(
            h='xi^3*(xi^2-1)^2/(1+xi^2)^(7/2)',
            phi_left='h*(1-tanh(xi))/2', phi_right='h*(1+tanh(xi))/2',
            chi='xi+ell_left*phi_left+ell_right*phi_right',
            squared_derivative_polynomial=str(certificate.as_expr()),
            nonnegative_real_root_count=int(root_count),
            sufficient_condition='abs(ell_left)+abs(ell_right)<=1 implies chi_prime>0 globally',
            status='Certificate for a candidate coordinate; data fitting and PDE closure are separate.'),
        caveats=['Every geometry row is measured, so this is not a finite-parameter 2D solution.',
                 'Master holdout is a transverse interval of the same saved field, not an independent run.',
                 'Native outer-row comparisons only check extension within the saved finite field.',
                 'No PDE identity or global Biot-Savart closure is asserted.',
                 'Algebraic-in-xi rational masters with sinh geometry have logarithmic physical tails; only use them in a finite core.'])
    (OUT/'tail_master_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(OUT/'tail_master_arrays.npz', xi=xi, learned=target,
                        outer_xi=outer_xi, heights=heights, rows=rows, valid=valid,
                        geometry=np.asarray(geometries))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.7), constrained_layout=True)
    ax[0].plot(xi, target, 'k', lw=2, label='measured master')
    for f in fits:
        curve = master(xi, f['alpha'], f['b'], f['epsilon'])
        label = f"alpha={f['alpha']:.3f}"
        ax[0].plot(xi, curve, label=label)
        ax[1].plot(xi, curve-target, label=label)
        ax[2].semilogy(outer_xi, master(outer_xi, f['alpha'], f['b'], f['epsilon']))
    for low, high in [(0, 1), (1, 64), (64, 1024)]:
        use = (heights >= low) & (heights <= high)
        ax[2].semilogy(outer_xi, np.nanmean(rows[use], axis=0), '--',
                       label=f'measured y={low}..{high}')
    ax[0].set(title='Different tails, nearly identical core', xlabel='xi', ylabel='q / peak')
    ax[0].legend(fontsize=7)
    ax[1].set(title='Error against measured master', xlabel='xi', ylabel='absolute normalized error')
    ax[2].set(title='Transverse extension (no refitting)', xlabel='xi', ylabel='q / peak')
    ax[2].legend(fontsize=7)
    for a in ax:
        a.grid(alpha=.2)
        a.axvspan(-1, 1, color='gray', alpha=.08)
    fig.savefig(OUT/'tail_master.png', dpi=170)
    plt.close(fig)
    print(json.dumps({'fits': [{k:v for k,v in f.items() if k != 'native_row_audits'} for f in fits],
                      'masked_points': report['masked_out_of_domain_points']}))


if __name__ == '__main__':
    main()
