"""Small exact checks for the three-function closure research note.

No saved fields are changed and no PDE evolution is run.
"""
from pathlib import Path
import json
import sympy as sp
import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/three_function_special_20260928'


def main():
    s, p, mu, lam = sp.symbols('s p mu lambda', real=True)
    R = 1 + s**2
    k = 1 / (mu * (6*p-2))
    J, K = R**p, k*R**(1-p)
    first = (mu*s*sp.diff(K, s)
             + (2*(lam-mu)-2*mu*s*sp.diff(J, s)/J)*K-1/J)
    first = sp.simplify(sp.expand_power_base(first.subs(lam, 3*p*mu), force=True))
    assert first == 0

    Jp, Jpp = sp.diff(J, s), sp.diff(J, s, 2)
    Kp, Kpp = sp.diff(K, s), sp.diff(K, s, 2)
    E1 = (mu*(1+K**2)*(s*Jp**2-J*Jp-s*J*Jpp)
          + (4*mu-2*lam)*J**2*K*Kp + 2*mu*s*J*K*Kp*Jp
          + 2*mu*s*J**2*K*Kpp - 3*mu*s*J**2*Kp**2 + 2*J*Kp)
    E1_expected = -4*mu*p*s*R**(2*p-2) + 4*mu*k**2*(1+p-3*p**2)*s
    E1_error = sp.simplify(sp.expand_power_base(
        E1.subs(lam, 3*p*mu)-E1_expected, force=True))
    assert E1_error == 0

    kk = sp.symbols('k', real=True)
    a2 = -(1+kk**2)*s/(5*R**2)
    mu1 = 1/(4*kk)
    rhs = (-4*mu1*s-4*mu1*kk**2*s)/(2*R**3)
    a2_error = sp.factor(mu1*s*sp.diff(a2, s)
                        + (9*mu1-6*mu1*s**2/R)*a2-rhs)
    assert a2_error == 0

    x, a, b, c, v = sp.symbols('x a b c v', real=True)
    Q = x-x**3/3
    q = sp.diff(Q, x)
    r = c-b+a*x+b*x**2
    T = (mu*Q+v)/q
    Ewall = T*sp.diff(r, x)+(2*lam-mu-sp.diff(T, x))*r+q**2
    e0 = sp.factor(Ewall.subs(x, 0))
    e1 = sp.factor(sp.diff(Ewall, x).subs(x, 0))
    wall_solution = sp.solve([e0, e1], [b, c])
    b_expected = -(1+v*a)/(2*(lam-mu))-(2*lam-mu)*a/(2*v)
    assert sp.factor(wall_solution[b]-b_expected) == 0

    # Second-normal-layer obstruction for the centered cubic wall and a
    # height-independent physical affine mother; a,w denote second jets.
    aa, ww = sp.symbols('a2 w2')
    obstruction = sp.Poly(sp.expand(
        aa*R**3 + sp.Rational(2, 3)*ww*s**2*R**2
        -sp.Rational(2, 5)*(1+kk**2)), s)
    high_solution = sp.solve([obstruction.coeff_monomial(s**6),
                              obstruction.coeff_monomial(s**4)], [aa, ww])
    remaining = sp.factor(obstruction.coeff_monomial(1).subs(high_solution))
    assert sp.expand(remaining + sp.Rational(2, 5)*(1+kk**2)) == 0

    mp.mp.dps = 60
    rows = []
    for pp, ratio in [('0.6', '2.4'), ('1', '3'), ('1.3', '3.9')]:
        pp, ratio = mp.mpf(pp), mp.mpf(ratio)
        mm = mp.mpf('0.7')
        gamma = 2*ratio-2
        for ss in [mp.mpf('0'), mp.mpf('0.3'), mp.mpf('2'), mp.mpf('-3')]:
            integral = (1+ss**2)**(2*pp)/mm * mp.quad(
                lambda t: t**(gamma-1)/(1+ss**2*t*t)**(3*pp), [0, 1])
            hyper = (1+ss**2)**(2*pp)/(mm*gamma)*mp.hyp2f1(
                3*pp, gamma/2, gamma/2+1, -ss**2)
            row = dict(p=str(pp), ratio=str(ratio), s=str(ss),
                       relative_error=str(abs(integral-hyper)/abs(integral)))
            if abs(ratio-3*pp) < mp.mpf('1e-55'):
                elementary = (1+ss**2)**(1-pp)/(mm*(6*pp-2))
                row['resonance_relative_error'] = str(abs(integral-elementary)/abs(integral))
            rows.append(row)
    assert max(mp.mpf(r['relative_error']) for r in rows) < mp.mpf('1e-55')
    result = dict(
        first_resonance_residual=str(first), second_source_residual=str(E1_error),
        cubic_second_layer_residual=str(a2_error),
        wall_peak_equations=[str(e0), str(e1)],
        wall_peak_slope=str(sp.factor(wall_solution[b])),
        affine_cubic_obstruction=str(remaining),
        resonance_and_tail_intersection=sp.solve(sp.Eq(3*p, 2*p+1), p)[0].__str__(),
        hypergeometric_checks=rows,
        scope='Local symbolic identities and quadrature checks; no existence, all-order convergence, or global steady solution claimed.')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: val for k, val in result.items() if k != 'hypergeometric_checks'}, indent=2))


if __name__ == '__main__':
    main()
