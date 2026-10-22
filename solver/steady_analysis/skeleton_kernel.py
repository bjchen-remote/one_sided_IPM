"""One-dimensional geometric Cauchy kernels and a velocity verification.

A full transverse profile must be supplied. No continuation of the measured
core-only integer master is silently chosen. P(.,0) has the point value 0;
its one-sided limits are +/- pi F, and height integration is split at h=y.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad, quad_vec

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/skeleton_closure_20260928'
BETA = 1.5


def G(xi, k):
    return (np.sinh(BETA*xi)+k*(np.cosh(BETA*xi)-1))/BETA


def inverse_G(z, k):
    if abs(k) >= 1:
        raise ValueError('A global coordinate requires abs(k)<1')
    theta = np.arctanh(k)
    return (np.arcsinh(np.sinh(theta)+BETA*z*np.cosh(theta))-theta)/BETA


def cauchy_kernel(delta, t, transverse):
    """Return P,H and quadrature error estimates for F(z)=transverse(z).

    Paired radial integration implements the Hilbert principal value at t=0.
    This scalar reference implementation favors explicit error reporting over
    speed; a nonlinear boundary solver should tabulate these kernels first.
    """
    d, a = float(delta), abs(float(t))
    if a == 0:
        P, ep = 0., 0.
    else:
        # r=a*tan(theta) integrates the narrow Poisson kernel exactly.
        P, ep = quad(lambda theta: transverse(d+a*np.tan(theta))
                     + transverse(d-a*np.tan(theta)), 0, np.pi/2,
                     epsabs=2e-10, epsrel=2e-10, limit=200)
        P *= np.sign(t)
    def paired(r):
        if r == 0:
            return 0.
        return (transverse(d-r)-transverse(d+r))*r/(r*r+a*a)
    cut = max(1., abs(d), a)
    H1, eh1 = quad(paired, 0, cut, epsabs=2e-10, epsrel=2e-10, limit=200)
    H2, eh2 = quad(paired, cut, np.inf, epsabs=2e-10, epsrel=2e-10, limit=200)
    return np.array([P, H1+H2]), np.array([ep, eh1+eh2])


def geometry_kernel(delta, t, k, full_master):
    return cauchy_kernel(delta, t, lambda z: full_master(inverse_G(z,k), k))


def lorentz_kernel(delta, t):
    a = 1+abs(t)
    return np.pi*np.array([np.sign(t)*a, delta])/(delta*delta+a*a)


def relative_velocity(x, y, skeleton, kernel=lorentz_kernel):
    """Integrate the geometric kernels over source height, with subtraction.

    skeleton(h) returns (m,A,W). This validation uses a transverse Lorentzian.
    The center subtraction remains inside the same integrand, as required
    when individual absolute velocities would have nonintegrable tails.
    """
    if y < 0:
        raise ValueError('Half-plane y>=0 required')
    def integrand(h):
        m,A,W = skeleton(h)
        d, d0 = (x-m)/W, -m/W
        pm, hm = kernel(d, (y-h)/W)
        pp, hp = kernel(d, (y+h)/W)
        p0,_ = kernel(d0,h/W)
        return A*np.array([pm-pp+2*p0, hp-hm])/(2*np.pi)
    if y == 0:
        value, error = quad_vec(integrand,0,np.inf,epsabs=2e-10,epsrel=2e-10)
    else:
        a, ea = quad_vec(integrand,0,y,epsabs=1e-10,epsrel=1e-10)
        b, eb = quad_vec(integrand,y,np.inf,epsabs=1e-10,epsrel=1e-10)
        value, error = a+b,ea+eb
    return value,float(error)


def run_checks():
    cases=[]
    for k in [0., .75, -.4]:
        # Nonlinear coordinate whose physical transverse density is known.
        f=lambda xi,kk: 1/(1+G(xi,kk)**2)
        for delta,t in [(.3,.01),(1.2,.7),(-.8,-.4),(.7,0.)]:
            got,err=geometry_kernel(delta,t,k,f)
            exact=lorentz_kernel(delta,t)
            cases.append(dict(k=k,delta=delta,t=t,computed=got.tolist(),
                              exact=exact.tolist(),max_error=float(max(abs(got-exact))),
                              quadrature_error=err.tolist()))
    # Bounded, nonhomogeneous skeleton; no PDE time evolution.
    skeleton=lambda h:(.3*h,np.exp(-h),np.sqrt(1+h))
    checks=[]
    for x,y in [(.3,.2),(-.7,.6),(1.2,1.1)]:
        m,A,W=skeleton(y);source=A/(1+((x-m)/W)**2)
        derivatives=[]
        for e in [1e-3,5e-4]:
            ux=(relative_velocity(x+e,y,skeleton)[0]-relative_velocity(x-e,y,skeleton)[0])/(2*e)
            uy=(relative_velocity(x,y+e,skeleton)[0]-relative_velocity(x,y-e,skeleton)[0])/(2*e)
            derivatives.append(dict(step=e,divergence=float(ux[0]+uy[1]),
                                    curl_plus_source=float(ux[1]-uy[0]+source)))
        checks.append(dict(x=x,y=y,q=source,derivatives=derivatives))
    wall=[dict(x=x,velocity=relative_velocity(x,0,skeleton)[0].tolist())
          for x in [-1,0,1]]
    result=dict(kernel_tests=cases,velocity_derivative_tests=checks,wall_tests=wall,
                maximum_kernel_error=max(c['max_error'] for c in cases),
                scope='Synthetic exact Lorentzian source in a nonhomogeneous skeleton; verifies the 1D velocity operator, not the unknown IPM steady profile.',
                restrictions=['A full transverse completion is mandatory.',
                              'Pointwise zero at t=0 differs from one-sided Poisson limits.',
                              'Height integration is split at the target height.',
                              'The wall and center subtraction are imposed by the kernel.'])
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'kernel_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    run_checks()
