"""Independent transform identities; no discrete steady-state fitting.

Verify the centered multiplier, radial integrating factor, pole cancellation,
and the exact Airy/Bessel transform of a cubic inverse wall profile.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import airy, kv


def symbolic_checks():
    px, py, kx, ky = sp.symbols('px py kx ky', real=True)
    m = sp.Matrix([px*py, -px**2])/(px**2+py**2)
    assert sp.simplify(sp.Matrix([px, py]).dot(m)) == 0
    assert sp.simplify(sp.Matrix([kx-px, ky-py]).dot(m)
                       - sp.Matrix([kx, ky]).dot(m)) == 0

    r, lam, alpha = sp.symbols('r lambda alpha', positive=True)
    v, nx = sp.symbols('v nx', real=True)
    h = sp.Function('h')(r)
    factor = r**(2+alpha)*sp.exp(-sp.I*v*nx*r/lam)
    lhs = lam*r*sp.diff(h, r)+(2*lam+alpha*lam-sp.I*v*nx*r)*h
    assert sp.simplify(sp.diff(factor*h, r)-factor*lhs/(lam*r)) == 0

    a, s, y = sp.symbols('a s y', positive=True)
    # Unit horizontal source amplitude, source exp(-y), |k|=a.
    psi = (sp.exp(-y)-sp.exp(-a*y))/(a**2-1)
    assert sp.simplify(a*a*psi-sp.diff(psi, y, 2)-sp.exp(-y)) == 0
    w = sp.simplify(-sp.diff(psi, y).subs(y, 0))
    assert sp.simplify(w+1/(a+1)) == 0
    transformed = sp.cancel((1/(s+1)+w)/(a*a-s*s))
    assert sp.simplify(transformed-1/((a+1)*(s+1)*(s+a))) == 0
    assert sp.simplify((1/(s+1)+w).subs(s, a)) == 0

    z, mu = sp.symbols('z mu', real=True)
    F = z+z**3/3
    wall_velocity = sp.expand(mu*z*sp.diff(F, z)-3*mu*F)
    assert wall_velocity == -2*mu*z
    tilt, k = sp.symbols('tilt k', real=True, positive=True)
    mismatch = 1/(k+sp.I*k*tilt)-2*mu/(sp.I*k)
    assert sp.simplify(sp.re(mismatch)-1/(k*(1+tilt**2))) == 0
    # The same positivity holds for all real tilt, not only positive tilt.

    ell = sp.symbols('ell', positive=True)
    jump_transform = -2*sp.I*ell/(1+ell**2)
    assert sp.limit(ell*jump_transform, ell, sp.oo) == -2*sp.I
    return dict(centered_multiplier='exact symbolic zero',
                radial_integrating_factor='exact symbolic zero',
                laplace_pole_cancellation='exact symbolic zero',
                wall_jump_leading_coefficient='-2i',
                cubic_required_wall_velocity=str(wall_velocity),
                constant_tilt_mismatch_real_part='1/(abs(k)*(1+tilt**2)) > 0')


def cubic_transform_checks():
    # x=F(z), q=dz/dx. The x integral is checked independently of the
    # oscillatory cubic-phase representation used to derive the Airy formula.
    def q(x):
        z = 2*np.sinh(np.arcsinh(1.5*x)/3)
        return 1/(1+z*z)

    rows = []
    for k in (0.1, 0.25, 0.5, 1., 2., 4.):
        value, err = quad(q, 0, np.inf, weight='cos', wvar=k,
                          epsabs=2e-11, limlst=250, limit=250)
        direct = 2*value
        ai = 2*np.pi*k**(-1/3)*airy(k**(2/3))[0]
        bessel = 2/np.sqrt(3)*kv(1/3, 2*k/3)
        assert abs(direct-ai) < 2e-9
        assert abs(ai-bessel) < 2e-12
        rows.append(dict(k=k, direct_horizontal_fourier=direct,
                         airy=ai, bessel=bessel, absolute_error=abs(direct-ai),
                         quadrature_error_estimate=2*err))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path('runs/fourier_20260928'))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    result = dict(symbolic=symbolic_checks(), cubic_transform=cubic_transform_checks(),
                  scope='Transform identities and a wall trial only; not an IPM steady solution.',
                  references=['https://dlmf.nist.gov/9.5.E1',
                              'https://dlmf.nist.gov/9.6.E2'])
    (args.out/'symbolic_checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
