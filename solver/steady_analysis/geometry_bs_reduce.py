"""Geometry-adapted Biot--Savart reduction and angular equation solves.

The integer master is used exactly on |xi| <= 3/2. Positive outer
continuations are constructed by two scalar wall-regularity equations,
not fitted to the unconverged outer data. All output is a conditional
asymptotic construction, not a global smooth steady profile.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_simpson, quad, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import gamma, hyp2f1

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/geometry_bs_20260928"
HALF = 1.5
K = 0.75
OMEGA = 0.5
GAMMA = 1.0


def g(xi):
    return (np.sinh(HALF*xi) + K*(np.cosh(HALF*xi)-1))/HALF


def gp(xi):
    return np.cosh(HALF*xi) + K*np.sinh(HALF*xi)


def xi_of_t(t):
    theta = np.arctanh(K)
    return (np.arcsinh(np.sinh(theta) + HALF*(t-GAMMA)/OMEGA*np.cosh(theta))-theta)/HALF


def rational(xi):
    return 32/(xi**6-3*xi**5+3*xi**4+3*xi**3+28*xi**2+32)


def blend(xi, end):
    """C4 cutoff; enough smoothness for the two corrections solved here."""
    z = np.clip((np.abs(xi)-HALF)/(end-HALF), 0, 1)
    return 1-(126*z**5-420*z**6+540*z**7-315*z**8+70*z**9)


def primitive(t, alpha):
    return t*hyp2f1(0.5, (1-alpha)/2, 1.5, -t*t)


class Seed:
    """P=y^alpha F(x/y), F(1)=1/4, F'(1)=1 with smooth wall traces."""

    def __init__(self, alpha, end):
        self.alpha, self.end = alpha, end
        self.constant = np.sqrt(np.pi)*gamma(-alpha/2)/(2*gamma((1-alpha)/2))
        self.ends = {s: GAMMA+OMEGA*g(s*end) for s in (-1, 1)}
        self.amplitudes, self.solutions, self.scalar_defects = {}, {}, {}
        for side in (-1, 1):
            end_xi = side*end
            def terms(z):
                t = GAMMA+OMEGA*g(z)
                w = blend(z, end)
                jac = OMEGA*gp(z)
                return w*rational(z)*jac, (1-w)*(1+t*t)**((alpha-1)/2)*jac
            i0 = quad(lambda z: terms(z)[0], 0, end_xi, epsabs=2e-12,
                      points=[side*HALF])[0]
            i1 = quad(lambda z: terms(z)[1], 0, end_xi, epsabs=2e-12,
                      points=[side*HALF])[0]
            target = primitive(self.ends[side], alpha)-side*self.constant
            amp = (0.25+i0)/(target-i1)
            if amp <= 0:
                raise ValueError(f"Nonpositive continuation: {alpha=}, {end=}, {side=}, {amp=}")
            self.amplitudes[side] = amp
            self.scalar_defects[side] = float(0.25+i0+amp*i1-amp*target)
            sol = solve_ivp(lambda z, _: [sum((terms(z)[0], amp*terms(z)[1]))],
                            [0, end_xi], [0.25], rtol=3e-12, atol=3e-13,
                            dense_output=True, max_step=.025)
            if not sol.success:
                raise RuntimeError(sol.message)
            self.solutions[side] = sol.sol

    def values(self, theta):
        alpha = self.alpha
        theta = np.asarray(theta)
        aa, ap, dd = (np.empty_like(theta) for _ in range(3))
        inside = (theta > 1e-4) & (theta < np.pi-1e-4)
        tt = 1/np.tan(theta[inside]); xx = xi_of_t(tt)
        ff = np.empty_like(tt)
        amps = np.where(xx < 0, self.amplitudes[-1], self.amplitudes[1])
        ww = blend(xx, self.end)
        qq = ww*rational(xx)+(1-ww)*amps*(1+tt*tt)**((alpha-1)/2)
        for side in (-1, 1):
            choose = (xx < 0) if side < 0 else (xx >= 0)
            middle = choose & (abs(xx) <= self.end)
            if np.any(middle):
                ff[middle] = self.solutions[side](xx[middle])[0]
            outer = choose & ~middle
            ff[outer] = self.amplitudes[side]*(primitive(tt[outer], alpha)-side*self.constant)
        sine, cosine = np.sin(theta[inside]), np.cos(theta[inside])
        aa[inside] = sine**alpha*ff
        dd[inside] = sine**(alpha-1)*qq
        ap[inside] = (alpha*cosine*aa[inside]-dd[inside])/sine
        # Outside the patch D is exactly constant, and A has a convergent
        # even wall series. This avoids subtracting large terms near theta=0,pi.
        for side in (-1, 1):
            edge = (~inside) & ((theta < np.pi/2) if side > 0 else (theta >= np.pi/2))
            distance = theta[edge] if side > 0 else np.pi-theta[edge]
            c = self.amplitudes[side]
            aa[edge] = side*c*(1/alpha + distance**2/(2*(alpha-2)))
            ap[edge] = c*distance/(alpha-2)
            dd[edge] = c
        return aa, ap, dd


def angular_green(theta, source, p):
    """Full Dirichlet inverse of -d_theta^2-(p+1)^2, nonresonant p."""
    nu = p+1
    if abs(nu) < 1e-12:
        left = cumulative_simpson(theta*source, x=theta, initial=0)
        right0 = cumulative_simpson((np.pi-theta)*source, x=theta, initial=0)
        right = right0[-1]-right0
        return ((np.pi-theta)*left+theta*right)/np.pi, (-left+right)/np.pi
    if abs(np.sin(np.pi*nu)) < 1e-9:
        raise ValueError("Resonance requires a compatibility/logarithmic solve")
    left = cumulative_simpson(np.sin(nu*theta)*source, x=theta, initial=0)
    right0 = cumulative_simpson(np.sin(nu*(np.pi-theta))*source, x=theta, initial=0)
    right = right0[-1]-right0
    b = (np.sin(nu*(np.pi-theta))*left+np.sin(nu*theta)*right)/(nu*np.sin(np.pi*nu))
    bp = (-np.cos(nu*(np.pi-theta))*left+np.cos(nu*theta)*right)/np.sin(np.pi*nu)
    return b, bp


def nonlinear(ai, api, bi, bpi, pi, aj, apj, pj):
    return -pj*bpi*aj+(pi+1)*bi*apj


def solve_hierarchy(alpha, end, ntheta):
    theta = np.linspace(0, np.pi, ntheta)
    seed = Seed(alpha, end)
    a0, ap0, d0 = seed.values(theta)
    a, ap, d, b, bp = [a0], [ap0], [d0], [], []
    powers = [(n+1)*alpha-n for n in range(3)]
    for n, p in enumerate(powers):
        if n:
            forcing = np.zeros_like(theta)
            for i in range(n):
                j = n-1-i
                forcing += nonlinear(a[i], ap[i], b[i], bp[i], powers[i], a[j], ap[j], powers[j])
            an = forcing/(n*(1-alpha))  # lambda=1 fixes the rate unit.
            anp = CubicSpline(theta, an)(theta, 1)
            a.append(an); ap.append(anp)
            d.append(p*np.cos(theta)*an-np.sin(theta)*anp)
        bn, bnp = angular_green(theta, d[n], p)
        b.append(bn); bp.append(bnp)
    a, ap, d, b, bp = map(np.asarray, (a, ap, d, b, bp))
    idx = slice(4, -4)
    poisson = []
    for n, p in enumerate(powers):
        # Independent differentiated evaluation, not the defining algebraic identity.
        defect = -CubicSpline(theta, bp[n])(theta, 1)-(p+1)**2*b[n]-d[n]
        poisson.append(float(np.max(abs(defect[idx]))/max(np.max(abs(d[n][idx])), 1e-14)))
    # Root inside the integer core, used to test exact homogeneous closure.
    root = brentq(lambda z: seed.solutions[-1](z)[0], -HALF, 0)
    tr = GAMMA+OMEGA*g(root); ang = np.arctan2(1, tr)
    br = float(CubicSpline(theta, b[0])(ang))
    apr = float(seed.values(np.array([ang]))[1][0])
    contradiction = (alpha+1)*br*apr
    nonlinear00 = nonlinear(a[0], ap[0], b[0], bp[0], alpha, a[0], ap[0], alpha)
    core = (xi_of_t(1/np.tan(theta[1:-1])) >= -HALF) & (xi_of_t(1/np.tan(theta[1:-1])) <= HALF)
    core = np.r_[False, core, False]
    audits = []
    for radius in (16., 32., 64., 128., 256.):
        levels = []
        for upto in range(3):
            residual = sum((powers[n]-alpha)*radius**powers[n]*a[n] for n in range(upto+1))
            for i in range(upto+1):
                for j in range(upto+1):
                    residual += radius**(powers[i]+powers[j]-1)*nonlinear(
                        a[i], ap[i], b[i], bp[i], powers[i], a[j], ap[j], powers[j])
            scale = radius**alpha*np.sqrt(np.trapezoid(a[0]**2, theta))
            levels.append(dict(order=upto, relative_transport_l2=float(
                np.sqrt(np.trapezoid(residual**2, theta))/scale),
                core_relative_transport_l2=float(np.linalg.norm(residual[core])/
                    (radius**alpha*np.linalg.norm(a[0,core])))))
        audits.append(dict(radius=radius, levels=levels))
    report = dict(alpha=alpha, transition_end_xi=end, angular_nodes=ntheta,
                rates=dict(lambda_=1, mu=alpha, effective_constant_drift=0),
                continuation_amplitudes=seed.amplitudes,
                wall_regularity_scalar_defects=seed.scalar_defects,
                density_powers=powers, poisson_relative_linf=poisson,
                pure_homogeneous_obstruction=dict(xi=root, x_over_y=float(tr),
                    theta=float(ang), streamfunction_angular_value=br,
                    density_angular_derivative=apr, nonlinear_at_density_zero=contradiction,
                    nonlinear_relative_l2=float(np.sqrt(np.trapezoid(nonlinear00**2,theta)) /
                        np.sqrt(np.trapezoid(a[0]**2,theta)))),
                transport_audits=audits,
                norms=dict(a=[float(np.linalg.norm(v)) for v in a],
                           b=[float(np.linalg.norm(v)) for v in b]))
    report['cartesian_reconstruction'] = cartesian_audit(theta,a,b,powers)
    report['minimum_scaled_q_r16'] = float(np.min(sum(16.**(p-alpha)*d[n] for n,p in enumerate(powers))))
    return report, (theta, a, ap, d, b, bp)


def cartesian_audit(theta, a, b, powers):
    """Differentiate reconstructed Cartesian fields independently in x,y."""
    aa=[CubicSpline(theta,v) for v in a]
    bb=[CubicSpline(theta,v) for v in b]
    def field(x,y,stream=False):
        r=np.hypot(x,y);t=np.arctan2(y,x)
        return sum(r**(p+int(stream))*(bb if stream else aa)[n](t) for n,p in enumerate(powers))
    first=np.array([1.,-8.,0.,8.,-1.])/12
    second=np.array([-1.,16.,-30.,16.,-1.])/12
    rows=[]
    r=64.;step=r*1e-3;offset=step*np.arange(-2,3)
    for angle in (.35,.7,.95,1.2,1.8,2.5):
        x,y=r*np.cos(angle),r*np.sin(angle)
        qx=np.dot(first,field(x+offset,y))/step
        qy=np.dot(first,field(x,y+offset))/step
        px=np.dot(first,field(x+offset,y,True))/step
        py=np.dot(first,field(x,y+offset,True))/step
        lap=(np.dot(second,field(x+offset,y,True))+np.dot(second,field(x,y+offset,True)))/step**2
        transport=(-py+x)*qx+(px+y)*qy-powers[0]*field(x,y)
        rows.append(dict(theta=angle,relative_poisson=float(abs(lap+qx)/max(abs(qx),1e-12)),
                         transport_scaled_by_r_alpha=float(transport/r**powers[0])))
    return rows


def symbolic_checks():
    # Check the constant-shape Poisson ODE from the saved report against
    # direct Cartesian differentiation in the coordinate t=x/y.
    x, y, p, t = sp.symbols('x y p t', real=True)
    f, h = sp.Function('f')(t), sp.Function('h')(t)
    psi = y**(p+1)*h.subs(t, x/y)
    lap = sp.diff(psi,x,2)+sp.diff(psi,y,2)
    expected = (1+t*t)*sp.diff(h,t,2)-2*p*t*sp.diff(h,t)+p*(p+1)*h
    assert sp.simplify(lap.subs(x,t*y)/y**(p-1)-expected) == 0
    S = sp.Function('S')(t); bb = sp.Function('b')(t); ss=sp.diff(S,t)
    saved_ode = -(sp.diff((S*S+4)/(2*ss)*sp.diff(bb,t)-(p+1)*S*bb/2,t)
                  +(p+1)/2*(-S*sp.diff(bb,t)+(p+1)*ss*bb))/(2*ss)
    direct_ode = -((S*S+4)*sp.diff(bb,t,2)/ss**2
                   -((S*S+4)*sp.diff(S,t,2)/ss**3+2*p*S/ss)*sp.diff(bb,t)
                   +p*(p+1)*bb)/4
    assert sp.simplify(saved_ode-direct_ode) == 0
    P = y**p*f.subs(t,x/y)
    adv = -sp.diff(psi,y)*sp.diff(P,x)+sp.diff(psi,x)*sp.diff(P,y)
    assert sp.simplify(adv.subs(x,t*y)/y**(2*p-1)
                      -(p*sp.diff(h,t)*f-(p+1)*h*sp.diff(f,t))) == 0
    nu, q, r = sp.symbols('nu q r', positive=True)
    radial = r**(p+1)/(nu**2-(p+1)**2)
    assert sp.simplify(-(sp.diff(radial,r,2)+sp.diff(radial,r)/r-nu**2/r**2*radial)-r**(p-1)) == 0
    z, alpha = sp.symbols('z alpha', real=True)
    A, B = sp.Function('A')(z), sp.Function('B')(z)
    N = -alpha*sp.diff(B,z)*A+(alpha+1)*B*sp.diff(A,z)
    a1,a2,a3,b2,b3 = sp.symbols('a1 a2 a3 b2 b3')
    at_zero = {A:0,B:0,sp.diff(B,z):0,sp.diff(A,z):a1,
               sp.diff(A,z,2):a2,sp.diff(A,z,3):a3,
               sp.diff(B,z,2):b2,sp.diff(B,z,3):b3}
    second = sp.expand(sp.diff(N,z,2)).xreplace(at_zero)
    assert sp.simplify(second-(1-alpha)*b2*a1) == 0
    return dict(cartesian_poisson="exact zero", saved_C12_ode="exact zero", cartesian_advection="exact zero",
                radial_mellin_multiplier="exact zero",
                homogeneous_zero_obstruction="N''=(1-alpha) B'' A' at A=B=B'=0")


def direct_green_manufactured(alpha, order):
    """Independent 2D centered half-plane Green integral, q=r^(alpha-1)sin(theta)."""
    angle = .9
    x, y = np.cos(angle), np.sin(angle)
    nodes, weights = leggauss(order)
    total = 0.
    # Split at the observation point. The logarithmic singularity is integrable.
    for alo, ahi in ((0.,angle), (angle,np.pi)):
        th = (alo+ahi)/2+(ahi-alo)/2*nodes
        aw = (ahi-alo)/2*weights
        for lo, hi in ((-72.,-32.),(-32.,-8.),(-8.,0.),(0.,8.),(8.,32.)):
            lograd = (lo+hi)/2+(hi-lo)/2*nodes
            rw = (hi-lo)/2*weights
            rr = np.exp(lograd)[:,None]
            sx, sy = rr*np.cos(th)[None,:], rr*np.sin(th)[None,:]
            dist = (x-sx)**2+(y-sy)**2
            # Subtraction removes precisely the center's uniform horizontal velocity.
            green = np.log1p(4*y*sy/dist)/(4*np.pi)-y*np.sin(th)[None,:]/(np.pi*rr)
            source_area = rr**(alpha+1)*np.sin(th)[None,:]
            total += float(np.sum(rw[:,None]*aw[None,:]*green*source_area))
    exact = np.sin(angle)/(1-(alpha+1)**2)
    return dict(alpha=alpha, order=order, direct=total, exact=float(exact),
                absolute_error=float(abs(total-exact)))


def radial_mode_check():
    """General source Green formulas: regular origin and decaying source."""
    out = []
    radius = np.geomspace(.02,20,4001)
    for n in (1,2,3):
        def sol(r):
            low = quad(lambda s:s**(n+1)*np.exp(-s),0,r,epsabs=1e-12)[0]
            if n == 1:
                return (low/r-r*(1-np.exp(-r)))/2
            high = quad(lambda s:s**(1-n)*np.exp(-s),r,np.inf,epsabs=1e-12)[0]
            return (r**(-n)*low+r**n*high)/(2*n)
        values=np.array([sol(r) for r in radius])
        spline=CubicSpline(radius,values)
        residual=-(spline(radius,2)+spline(radius,1)/radius-n*n*values/radius**2)-np.exp(-radius)
        out.append(dict(mode=n, poisson_absolute_linf=float(np.max(abs(residual[4:-4])))))
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=OUT)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    report=dict(status="conditional angular solves; no global smooth steady state found",
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                symbolic=symbolic_checks(), solutions=[], refinement=[],
                direct_bs_checks=[], general_radial_green_checks=radial_mode_check())
    arrays={}
    for alpha in (2/5,3/8):
        for end in (3.,4.):
            coarse, carr=solve_hierarchy(alpha,end,2001)
            fine, farr=solve_hierarchy(alpha,end,4001)
            report['solutions'].append(fine)
            differences={}
            for k,name in ((1,'density'),(4,'streamfunction'),(5,'radial_velocity')):
                differences[name]=[float(np.linalg.norm(carr[k][n]-farr[k][n,::2])/
                    np.linalg.norm(farr[k][n,::2])) for n in range(3)]
            report['refinement'].append(dict(alpha=alpha,transition_end_xi=end,relative_l2=differences))
            prefix=f'a{alpha:g}_end{end:g}'
            for key,value in zip(('theta','A','A_prime','D','B','B_prime'),farr):
                arrays[prefix+'_'+key]=value
            print(json.dumps(dict(alpha=alpha,end=end,obstruction=fine['pure_homogeneous_obstruction'],
                                  refinement=differences)),flush=True)
        for order in (96,192,384):
            report['direct_bs_checks'].append(direct_green_manufactured(alpha,order))
    report['scope']=[
        'Exact integer master only on |xi|<=1.5; two positive wall-regular outer completions.',
        'All angular solves include both wall endpoints, not arbitrary finite-xi Dirichlet data.',
        'Higher negative radial degrees are exterior terms: they require core matching and possible harmonic multipoles.',
        'Two corrections solve their coefficient equations; finite series is not a global steady solution.',
        'Full radial mode Green formulas give a route to impose core and tail together; nonlinear radial system remains unsolved.',
        'No time evolution, original-data rewrite, or finite-box Newton solve was performed.']
    report['tail_sensitivity'] = []
    for alpha in (2/5,3/8):
        left,right=[s for s in report['solutions'] if s['alpha']==alpha]
        key0=f'a{alpha:g}_end3';key1=f'a{alpha:g}_end4'
        b0=arrays[key0+'_B'];b1=arrays[key1+'_B']
        report['tail_sensitivity'].append(dict(alpha=alpha,
            leading_streamfunction_relative_l2=float(np.linalg.norm(b1[0]-b0[0])/np.linalg.norm(b0[0])),
            zero_ray_nonlinearity_relative_change=float(abs(
                right['pure_homogeneous_obstruction']['nonlinear_at_density_zero']/
                left['pure_homogeneous_obstruction']['nonlinear_at_density_zero']-1))))
    (args.out/'checks.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.out/'angular_solutions.npz',**arrays)
    from scipy.io import savemat
    # MATLAB field names may not contain dots.
    savemat(args.out/'angular_solutions.mat',{k.replace('.','p'):v for k,v in arrays.items()},
            long_field_names=True,do_compression=True)
    from scipy.io import loadmat
    loaded=loadmat(args.out/'angular_solutions.mat')
    assert all(np.allclose(loaded[k.replace('.','p')].squeeze(),v) for k,v in arrays.items())
    make_plot(report,arrays,args.out)
    print(json.dumps(dict(output=str(args.out),direct_bs_checks=report['direct_bs_checks'])),flush=True)


def make_plot(report,arrays,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for row,alpha in enumerate((.4,.375)):
        for end,style in ((3,'-'),(4,'--')):
            key=f'a{alpha:g}_end{end}'
            theta=arrays[key+'_theta'];bb=arrays[key+'_B']
            axs[row,0].plot(theta,bb[0],style,label=f'outer join |xi|={end}')
        axs[row,0].axvline(.94845931985,color='gray',lw=.8,label='density zero ray')
        axs[row,0].axhline(0,color='gray',lw=.6)
        axs[row,0].set(xlabel='theta',ylabel='B0(theta)',title=f'Full angular BS, alpha={alpha:g}')
        axs[row,0].legend(fontsize=8)
        result=next(s for s in report['solutions'] if s['alpha']==alpha and s['transition_end_xi']==3)
        radii=[z['radius'] for z in result['transport_audits']]
        for order in range(3):
            error=[z['levels'][order]['relative_transport_l2'] for z in result['transport_audits']]
            axs[row,1].loglog(radii,error,'o-',label=f'{order} nonlinear corrections')
        axs[row,1].set(xlabel='r (lambda=1, peak coefficient=1)',ylabel='scaled transport L2',
                       title='Reconstructed residual, conditional outer series')
        axs[row,1].set_xticks(radii,labels=[str(int(r)) for r in radii])
        axs[row,1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        axs[row,1].legend(fontsize=8)
    for ax in axs.ravel():ax.grid(alpha=.2)
    fig.savefig(out/'angular_closure.png',dpi=180);plt.close(fig)


if __name__ == '__main__':
    main()
