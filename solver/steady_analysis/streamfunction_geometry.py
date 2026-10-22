"""Measure the geometry of centred psi, psi/y, and (psi/y)_x.

Read saved fields only (the source exporter may restore a static Poisson
field). No time evolution, domain extrapolation, or steady-state certification.
Run with the project's scipy/matplotlib dependencies on PYTHONPATH.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ipm_streamfunction_mpl")
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline

from streamfunction_sources import load_native, load_conditional, CONDITIONAL_NAMES

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/streamfunction_geometry_20260928"
NATIVE_HEIGHTS = np.array([0, .25, .5, 1, 2, 4, 8, 16, 32, 64, 128, 256.])
INNER_HEIGHTS = np.array([0, .025, .05, .1, .25, .5, .75, 1., 1.5, 2.])
BETA = 1.5


def scalar(x):
    return float(x) if x is not None and np.isfinite(x) else None


def relative(error, reference):
    den = np.linalg.norm(reference)
    return float(np.linalg.norm(error) / den) if den else None


class Fields:
    def __init__(self, p):
        self.p = p
        self.psi = RectBivariateSpline(p['y'], p['x'], p['psi'])
        self.q = RectBivariateSpline(p['y'], p['x'], p['q'])
        self.u1 = RectBivariateSpline(p['y'], p['x'], p['u1'])
        self.u2 = RectBivariateSpline(p['y'], p['x'], p['u2'])
        self.wall = CubicSpline(p['x'], -p['u1'][0], extrapolate=False)

    def H(self, x, y):
        return self.wall(x) if y == 0 else self.psi.ev(y + np.zeros_like(x), x) / y

    def D(self, x, y):
        return self.wall(x, 1) if y == 0 else self.u2.ev(y + np.zeros_like(x), x) / y

    def D_from_psi(self, x, y):
        return self.wall(x, 1) if y == 0 else self.psi.ev(y + np.zeros_like(x), x, dy=1) / y


def peak_and_widths(x, values, peak_window=None):
    s = CubicSpline(x, values, extrapolate=False)
    crit = s.derivative().roots(extrapolate=False)
    crit = crit[(crit > x[0]) & (crit < x[-1])]
    crit = crit[s(crit, 2) < 0]
    if peak_window is not None:
        crit = crit[(crit >= peak_window[0]) & (crit <= peak_window[1])]
    if not len(crit):
        return None
    m = float(crit[np.argmax(s(crit))]); A = float(s(m))
    if A <= 0:
        return None
    widths = {}
    for level in (.25, .5, .75):
        rr = CubicSpline(x, values - level * A).roots(extrapolate=False)
        left, right = rr[rr < m], rr[rr > m]
        widths[str(level)] = dict(left=scalar(m-left[-1]) if len(left) else None,
                                 right=scalar(right[0]-m) if len(right) else None)
    L, R = widths['0.5']['left'], widths['0.5']['right']
    k = (R-L)/(R+L)/np.tanh(BETA/2) if L is not None and R is not None else None
    W = (R+L)*BETA/(2*np.sinh(BETA)) if L is not None and R is not None else None
    return dict(m=m, A=A, curvature=float(s(m, 2)),
                curvatureWidth=float(np.sqrt(-2*A/s(m, 2))),
                widths=widths, W=W, k=k, peakSelectionWindow=peak_window,
                endpointExceedsCorePeak=bool(max(values[0],values[-1])>A),
                isGlobalOnSampledSupport=bool(A>=max(values)))


def measure(f, h):
    p = f.p
    # A measured window sufficiently wide for both tails of D; never extrapolate.
    radius = max(6., 8*h) if p['kind'] == 'finite_time_native' else np.inf
    lo, hi = max(p['x'][0], -radius), min(p['x'][-1], radius)
    # Do not union uniform samples with saved nodes: nominally identical
    # coordinates may differ by 1e-15, creating spurious cubic-spline extrema.
    x = np.linspace(lo, hi, 2401)
    H, D = f.H(x, h), f.D(x, h)
    q = f.q.ev(np.full_like(x, h), x)
    peak_window = [-3.,3.] if p['kind'] != 'finite_time_native' else None
    dg, qg = peak_and_widths(x, D, peak_window), peak_and_widths(x, q)
    dg_psi = peak_and_widths(x, f.D_from_psi(x,h), peak_window)
    roots = CubicSpline(x, H).roots(extrapolate=False)
    target = qg['m'] if qg is not None else 0.
    z = float(roots[np.argmin(abs(roots-target))]) if len(roots) else None
    u = f.u1.ev(np.full_like(x, h), x)
    ur = CubicSpline(x, u).roots(extrapolate=False)
    uz = float(ur[np.argmin(abs(ur-target))]) if len(ur) else None
    ell, P, a = p['ell'], p['P'], p['a']
    result = dict(y=float(h), Y=float(ell*h), xSupport=[float(lo), float(hi)],
                  psiZero=z, u1Zero=uz, allPsiZeros=roots.tolist(),
                  HMin=float(H.min()), HMax=float(H.max()),
                  DMin=float(D.min()), DMax=float(D.max()), D=dg, q=qg,
                  DpeakFromPsiDerivative=dg_psi['m'] if dg_psi is not None else None,
                  DpeakPositionDerivativeSensitivity=dg_psi['m']-dg['m'] if dg_psi is not None and dg is not None else None,
                  native=dict(psiZero=scalar(a+ell*z) if z is not None else None,
                              u1Zero=scalar(a+ell*uz) if uz is not None else None))
    if dg is not None and qg is not None:
        result['ratios'] = dict(DpeakOverQpeak=dg['A']/qg['A'],
            DatQpeakOverQpeak=float(f.D(np.array(qg['m']), h))/qg['A'],
            zeroMinusQpeak=scalar(z-qg['m']) if z is not None else None,
            DpeakMinusQpeak=dg['m']-qg['m'])
        for key, g in [('D', dg), ('q', qg)]:
            result['native'][key] = dict(m=a+ell*g['m'], A=P*g['A'],
                curvatureWidth=ell*g['curvatureWidth'],
                widths={lev:{side:scalar(ell*value) if value is not None else None
                        for side,value in w.items()} for lev,w in g['widths'].items()},
                W=scalar(ell*g['W']) if g['W'] is not None else None, k=g['k'])
    return result


def derivative_audit(f):
    x = np.linspace(max(-2., f.p['x'][0]), min(2., f.p['x'][-1]), 151)
    y = np.linspace(.05, min(1., f.p['y'][-1]), 101)
    X, Y = np.meshgrid(x, y)
    p1 = -f.psi.ev(Y, X, dx=1)
    p2 = f.psi.ev(Y, X, dy=1)
    qp = -f.psi.ev(Y, X, dx=2)-f.psi.ev(Y, X, dy=2)
    u1, u2, q = f.u1.ev(Y, X), f.u2.ev(Y, X), f.q.ev(Y, X)
    # Same scalar interpolant gives zero divergence analytically, so don't report
    # that identity as an independent accuracy test.
    return dict(region=dict(x=[float(x[0]),float(x[-1])], y=[float(y[0]),float(y[-1])]),
                psiSplineVsSavedU1=relative(p1-u1,u1),
                psiSplineVsSavedU2=relative(p2-u2,u2),
                psiSplinePoissonVsSavedQ=relative(qp-q,q),
                wallPsiMax=float(abs(f.p['psi'][0]).max()),
                HOrigin=float(f.wall(0.)),
                wallSplineMinusSolverH=float(f.psi.ev(0.,0.,dx=1)-f.wall(0.)))


def wall_geometry(f):
    s = float(f.wall(0., 1)); q0 = float(f.q.ev(0., 0.))
    heights = np.array([0., .005, .01, .02, .04])
    rows = [measure(f, h) for h in heights]
    z, u = ([r[k] for r in rows] for k in ('psiZero', 'u1Zero'))
    # The spline wall derivative and saved FD wall trace differ slightly.
    # Fit positive heights with a free intercept, avoiding an artificial jump
    # between these two numerical representations at y=0.
    zcoef = np.polynomial.polynomial.polyfit(heights[1:], z[1:], 2)
    ucoef = np.polynomial.polynomial.polyfit(heights[1:], u[1:], 2)
    return dict(s=s, qOrigin=q0,
        predictedPsiZeroSlope=q0/(2*s), predictedU1ZeroSlope=q0/s,
        fittedPsiZeroSlope=float(zcoef[1]), fittedU1ZeroSlope=float(ucoef[1]),
        fitHeights=heights[1:].tolist(), fitIntercept=float(zcoef[0]),
        HyWallFromPsi=float(f.psi.ev(0.,0.,dx=2)/2), HyWallFromQ=-q0/2,
        quadraticPsi='psi = s*x*y - q(0,0)*y^2/2 + O((abs(x)+y)^3)',
        scope='A local Poisson identity, not a steady transport assumption')


def power_fits(rows):
    r = [z for z in rows if 4-1e-9 <= z['Y'] <= 256+1e-9]
    Y = np.array([z['Y'] for z in r]); train = np.array([0,2,4,6]); test = np.array([1,3,5])
    cols = dict(psiZero=np.array([z['native']['psiZero'] for z in r]),
                Dpeak=np.array([z['native']['D']['m'] for z in r]),
                Damplitude=np.array([z['native']['D']['A'] for z in r]),
                DleftHalfWidth=np.array([z['native']['D']['widths']['0.5']['left'] for z in r]),
                DrightHalfWidth=np.array([z['native']['D']['widths']['0.5']['right'] for z in r]))
    out = {}
    for key, values in cols.items():
        slope, intercept = np.polyfit(np.log(Y[train]), np.log(values[train]), 1)
        pred = np.exp(intercept)*Y**slope
        out[key] = dict(coefficient=float(np.exp(intercept)), exponent=float(slope),
            trainRelativeL2=relative(pred[train]-values[train],values[train]),
            heldoutRelativeL2=relative(pred[test]-values[test],values[test]),
            heldoutMaxRelative=float(max(abs(pred[test]/values[test]-1))),
            observed=values.tolist(), predicted=pred.tolist())
    return dict(form='c*Y^p; empirical on native Y in [4,256], not an asymptotic exponent',
                Y=Y.tolist(), trainY=Y[train].tolist(), heldoutY=Y[test].tolist(), fits=out)


def topology(f):
    x = np.linspace(max(-3, f.p['x'][0]), min(3,f.p['x'][-1]),301)
    y = np.linspace(0,1,101)
    D = np.array([f.D(x,float(h)) for h in y])
    return dict(region=dict(x=[float(x[0]),float(x[-1])], y=[0.,1.]),
        minimumD=float(D.min()), sampledDStrictlyPositive=bool(np.all(D>0)),
        interpretation='Positive H_x means monotone H and u2>0 for y>0 in this measured core; no interior stagnation there. The wall origin is a boundary saddle.')


def plots(profiles_by_tau, records):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                         'savefig.facecolor':'white'})
    p = profiles_by_tau[20]; f = Fields(p)
    fig, ax = plt.subplots(1,3,figsize=(15.6,4.9),layout='constrained')
    x=np.linspace(-3,4,401); y=np.linspace(0,2,241)
    Z=f.psi(y,x); lim=np.quantile(abs(Z),.96)
    ax[0].pcolormesh(x,y,Z,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto',rasterized=True)
    ax[0].contour(x,y,Z,levels=np.linspace(-lim,lim,17),colors='#333333',linewidths=.45,alpha=.65)
    yy=np.linspace(0,2,81); rr=[measure(f,float(h)) for h in yy]
    for key,label,color,style in [('psiZero',r'$\psi=0$','#75214d','-'),
                                  ('u1Zero',r'$u_1=0$','#ffbd54','--')]:
        ax[0].plot([r[key] for r in rr],yy,style,color=color,lw=2.2,label=label)
    ax[0].plot([r['q']['m'] for r in rr],yy,':',color='#75ed9f',lw=2.4,label=r'$q$ peak')
    ax[0].plot([r['D']['m'] for r in rr],yy,'-',color='#f598db',lw=2.2,label=r'$u_2/y$ peak')
    ax[0].set(xlim=(-3,4),ylim=(0,2),xlabel='inner x',ylabel='inner y',title='Centred streamfunction: near-wall geometry')
    ax[0].legend(fontsize=8,loc='upper left',facecolor='#454545',labelcolor='white')

    yy=np.geomspace(4,256,151); eta=np.linspace(-3,3,241)
    rr=[measure(f,float(Y/p['ell'])) for Y in yy]
    H=np.array([f.H((eta*Y-p['a'])/p['ell'],Y/p['ell'])*p['P']*p['ell'] for Y in yy])
    A=np.array([r['native']['q']['A'] for r in rr])
    normalized=H/(yy[:,None]*A[:,None]); lim=np.quantile(abs(normalized),.97)
    im=ax[1].pcolormesh(eta,yy,normalized,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto',rasterized=True)
    for key,label,color,style in [('psiZero',r'$\psi=0$','#75214d','-'),('u1Zero',r'$u_1=0$','#ffbd54','--')]:
        ax[1].plot([r['native'][key]/Y if r['native'][key] is not None else np.nan for r,Y in zip(rr,yy)],yy,style,color=color,lw=2)
    ax[1].plot([r['native']['q']['m']/Y for r,Y in zip(rr,yy)],yy,':',color='#75ed9f',lw=2.5)
    ax[1].plot([r['native']['D']['m']/Y for r,Y in zip(rr,yy)],yy,color='#f598db',lw=2)
    ax[1].set(xlabel='native X/Y',ylabel='native Y (log scale)',yscale='log',title=r'$H/(Y A_q)$ across measured heights')
    fig.colorbar(im,ax=ax[1],shrink=.65,pad=.025)

    for Y in (4,16,64,256):
        r=measure(f,Y/p['ell']); g=r['D']; w=g['widths']['0.5']; width=(w['left']+w['right'])/2
        t=np.linspace(-2,2,401); xx=g['m']+width*t
        ax[2].plot(t,f.D(xx,Y/p['ell'])/g['A'],label=f'Y={Y}')
    ax[2].axhline(.5,color='#aaa',lw=.7);ax[2].axvline(0,color='#aaa',lw=.7)
    ax[2].set(xlabel=r'$(x-n_D)/[(L_D+R_D)/2]$',ylabel=r'$(u_2/y)/B$',title='Velocity slope: broad left wing',ylim=(0,1.05))
    ax[2].legend(fontsize=9)
    fig.suptitle('Observed streamfunction skeleton | native tau = 20',fontsize=15)
    fig.savefig(OUT/'streamfunction_geometry.png',dpi=180)
    fig.savefig(OUT/'streamfunction_geometry.pdf');plt.close(fig)

    fig, ax=plt.subplots(2,3,figsize=(14.5,8),layout='constrained')
    colors={10:'#9296bf',15:'#41968b',20:'#192e54'}
    for tau,p in profiles_by_tau.items():
        r=records[f'native_tau{tau}']['nativeRows']; r=[v for v in r if v['Y']>=4-1e-9]
        Y=np.array([v['Y'] for v in r]); label=f'tau={tau}'
        ax[0,0].semilogx(Y,[v['native']['psiZero']/v['Y'] for v in r],'-o',color=colors[tau],label=label)
        ax[0,1].semilogx(Y,[v['native']['D']['m']/v['Y'] for v in r],'-o',color=colors[tau],label=label)
        ax[0,2].semilogx(Y,[v['ratios']['DpeakOverQpeak'] for v in r],'-o',color=colors[tau],label=label)
        ax[1,0].semilogx(Y,[v['D']['widths']['0.5']['left']/v['D']['widths']['0.5']['right'] for v in r],'-o',color=colors[tau],label=label)
        ax[1,1].semilogx(Y,[v['D']['k'] for v in r],'-o',color=colors[tau],label=label)
        ax[1,2].semilogx(Y,[(v['D']['widths']['0.5']['left']+v['D']['widths']['0.5']['right'])/(v['q']['widths']['0.5']['left']+v['q']['widths']['0.5']['right']) for v in r],'-o',color=colors[tau],label=label)
    titles=[r'Zero-line slope $Z_H/Y$',r'Velocity-ridge slope $n_D/Y$',r'Amplitude ratio $B/A_q$',
            r'Left/right half-width ratio of $D$',r'Sinh skew of $D$ ($\beta=1.5$)',r'Width of $D$ / width of $q$']
    for a,title in zip(ax.ravel(),titles):
        a.set(xlabel='native Y',title=title);a.grid(alpha=.15)
    ax[0,0].legend()
    fig.suptitle('Geometry changes with height and time; no far-field limit assumed',fontsize=14)
    fig.savefig(OUT/'streamfunction_parameters.png',dpi=180);plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    native={tau:load_native(tau) for tau in (10,15,20)}
    records={}
    for p in [*native.values(),*[load_conditional(n) for n in CONDITIONAL_NAMES]]:
        f=Fields(p)
        r=dict(name=p['name'],kind=p['kind'],source=p['source'],gauge=p['gauge'],
            a=p['a'],P=p['P'],ell=p['ell'],U0=p['U0'],
            innerRows=[measure(f,float(h)) for h in INNER_HEIGHTS],
            derivativeAudit=derivative_audit(f),wallGeometry=wall_geometry(f),topology=topology(f))
        if p['kind']=='finite_time_native':
            r['nativeRows']=[measure(f,float(Y/p['ell'])) for Y in NATIVE_HEIGHTS]
            r['powerFits']=power_fits(r['nativeRows'])
        records[p['name']]=r
        print(p['name'],json.dumps(r['wallGeometry']),flush=True)
    report=dict(definitions=dict(H='centred psi/y; H(x,0)=-u1(x,0)',D='H_x=u2/y',
        identities=['u1=-H-y*H_y','u2=y*H_x','q=-2*H_y-y*(H_xx+H_yy)'],
        coordinates='Native X=a+ell*x,Y=ell*y; psi_native=P*ell^2*psi_inner, H_native=P*ell*H_inner, D_native=P*D_inner.'),
        methods=dict(interpolation='bicubic field splines; cubic horizontal root finding; H=psi/y with saved velocity wall trace; D uses saved u2/y and psi derivative is an independent sensitivity check',
            noExtrapolation=True,beta=BETA,
            fullSupport='native windows |x|<=max(6,8*y) intersect saved domain; conditional roots use saved box',
            missingWidths='null if crossing lies outside measured support'), records=records,
        limitations=['Finite-time geometry is not a converged steady state.',
            'Conditional roots freeze external response; box dependence is part of the findings.',
            'Power laws fit a finite height interval and are not asymptotic exponents.',
            'No independent full-domain Biot-Savart recomputation from a fitted skeleton has been done.',
            'Contours are Darcy streamlines in the centred frame; the effective similarity velocity has nonzero divergence and is not generated by this streamfunction.'])
    (OUT/'geometry.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    plots(native,records)
    print('Saved',OUT/'geometry.json',flush=True)


if __name__=='__main__':
    main()
