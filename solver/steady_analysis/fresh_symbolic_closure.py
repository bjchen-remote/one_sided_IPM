"""Exact local closure certificates, independent of fitted profile labels.

Run with SymPy available.  No fitted decimal is promoted to an exact constant.
This script proves local transport/div/curl identities and wall obstructions;
it does NOT certify the prescribed centered half-plane Biot--Savart operator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


def full_inverse_residual(X, s, y, lam, mu):
    """Return the full two-variable rational certificate for x=X(s,y).

    J>0 is a separate hypothesis.  When X is rational and the y primitive
    remains rational, numerator(E)=0 is a finite coefficient certificate.
    Physical curl(V)+q equals E/J**3, not merely its wall restriction.
    """
    J = sp.diff(X, s)
    K = sp.diff(X, y)
    H = sp.integrate(J, (y, 0, y))
    F = (2 * lam - mu) * H - mu * s * sp.diff(H, s)
    E = (
        J * (1 + K**2) * sp.diff(F, s)
        + F * (2 * J * K * sp.diff(K, s) - (1 + K**2) * sp.diff(J, s)
               - J**2 * sp.diff(K, y))
        - 2 * (lam - mu) * J**3 * K
        + 2 * mu * s * J**2 * K * sp.diff(J, s)
        - mu * s * J**3 * sp.diff(K, s)
        + J**2
    )
    return {"J": J, "K": K, "H": H, "F": F, "E": E,
            "V1": mu * s * J + F * K / J, "V2": F / J}


def cleared_coefficients(expr, variables):
    """Finite sufficient/necessary gate where the denominator is nonzero."""
    numerator, denominator = sp.cancel(expr).as_numer_denom()
    poly = sp.Poly(numerator, *variables)
    return denominator, [(list(monom), str(sp.factor(coeff)))
                         for monom, coeff in poly.terms()]


def verify_full_formula():
    """Independent coordinate-differentiation check on a nontrivial map."""
    s, y, lam, mu = sp.symbols("s y lambda mu")
    X = s + s**3 / 3 + y / 2 + s*y**2/7 + s**2*y/11
    obj = full_inverse_residual(X, s, y, lam, mu)
    J, K, V1, V2 = [obj[k] for k in ("J", "K", "V1", "V2")]
    dx = lambda f: sp.diff(f, s)/J
    dy = lambda f: sp.diff(f, y)-K*sp.diff(f, s)/J
    identities = {
        "transport": V1/J-V2*K/J-mu*s,
        "divergence": dx(V1)+dy(V2)-2*lam,
        "curl_certificate": dx(V2)-dy(V1)+1/J-obj["E"]/J**3,
        "wall_no_penetration": V2.subs(y, 0),
    }
    return {name: bool(sp.cancel(value) == 0)
            for name, value in identities.items()}


def sos_wall_certificates():
    """Exhaust the real SOS-quartic J / polynomial-degree<=5 K wall gate.

    J=C[1+(p0+p1*s+p2*s**2)**2], C>0, p2!=0.
    Set kappa=mu*C*K. Highest degree fixes lambda/mu=(10-d)/2,
    where d=degree K in {0,...,5}. Low-to-high elimination only divides by
    positive j0 and positive n+2*lambda/mu-2, so no real branch is lost.
    Saturation p2*w-1 enforces a genuine quartic J.
    """
    s = sp.Symbol("s")
    j = sp.symbols("j0:5")
    p0, p1, p2, w = sp.symbols("p0 p1 p2 w")
    substitution = dict(zip(j, [1+p0**2, 2*p0*p1, p1**2+2*p0*p2,
                                2*p1*p2, p2**2]))
    J = sum(j[i]*s**i for i in range(5))
    reports = []
    for degree in range(6):
        ratio = sp.Rational(10-degree, 2)
        A = 2*ratio-2
        kappa = [1/(A*j[0])]
        for n in range(1, degree+1):
            kappa.append(sp.factor(-sum(
                j[i]*kappa[n-i]*(n+A-3*i)
                for i in range(1, min(n, 4)+1)
            )/((n+A)*j[0])))
        K = sum(kappa[i]*s**i for i in range(degree+1))
        residual = sp.Poly(sp.expand(s*J*sp.diff(K,s)+A*J*K
                                    -2*s*sp.diff(J,s)*K-1), s)
        constraints = [sp.factor(sp.together(
            residual.coeff_monomial(s**n)).as_numer_denom()[0]
            .subs(substitution)) for n in range(degree+1, degree+4)]
        basis = sp.groebner(constraints+[w*p2-1], w, p2, p1, p0)
        univariate = [sp.factor(g) for g in basis
                      if g.free_symbols <= {p0}]
        if degree in (0, 4):
            real_branch = "p0=p1=0; p2 nonzero"
        else:
            real_branch = "none: strictly positive polynomial in p0 is forced to zero"
        reports.append({
            "K_degree": degree, "lambda_over_mu": str(ratio),
            "constraints": [str(c) for c in constraints],
            "groebner_basis": [str(g) for g in basis],
            "univariate_witnesses": [str(g) for g in univariate],
            "real_branch": real_branch,
        })
    return reports


def normalized_saturated_wall_certificates():
    """The full five-function family on its strictly normalized wall branch.

    J=1+(Q+e*Q**2/sqrt(1+Q**2))**2, 0<abs(e)<1.
    Removal of the independent atan/asinh terms leaves the K below. The
    All five wall derivatives were retained before eliminating logs.  Without
    all strict wall normalizations the unconstrained fitted family is not
    excluded by this routine alone. The e=0 branch is not handled here.
    """
    q, R, e, h, r, D, C, T, w = sp.symbols("q R e h r D C T w")
    J = 1+q**2+2*e*q**3*R/(1+q**2)+e**2*q**4/(1+q**2)
    K = D+C*q+T*(q**3+e*((q**2-2)*R+2))/3
    derivative = lambda f: sp.diff(f,q)+q/R*sp.diff(f,R)
    residual = ((q+h)*derivative(K)*J+2*(r-1)*K*J
                -2*(q+h)*K*derivative(J)-1)
    numerator = sp.together(residual).as_numer_denom()[0]
    remainder = sp.rem(numerator, R**2-1-q**2, R)
    coefficients = sp.Poly(remainder,q,R).coeffs()
    reports = []
    for ratio, substitutions in ((sp.Integer(3), {T:0,C:0}),
                                  (sp.Rational(5,2), {T:0}),
                                  (sp.Rational(3,2), {})):
        constraints = [sp.factor(c.subs({r:ratio, **substitutions}))
                       for c in coefficients]
        constraints = [c for c in constraints if c != 0]
        basis = sp.groebner(constraints+[w*e*(1-e**2)-1], w,D,C,T,h,e)
        assert list(basis) == [1]
        reports.append({"lambda_over_mu":str(ratio),
                        "forced_tangent_coefficients":{str(k):str(v) for k,v in substitutions.items()},
                        "constraints":[str(c) for c in constraints],
                        "groebner_basis":[str(g) for g in basis]})
    return {"J":str(J), "mu_times_K":str(K), "density_shift":"h=v/mu",
            "scope":"full five-function family under strict Q=0,q=1,q_x=0,q_xx=-2 wall normalization, 0<abs(e)<1; all five height derivatives retained before exact log elimination; e=0 not excluded",
            "normalization_note":"If Q=0,q=1,q_x=0,q_xx=-2 hold at the same wall point and a>0 in the original full family, they force b=0,C=1,a=1 without a change of coordinates.",
            "branches":reports}


def cubic_forced_curvature():
    """Positive cubic wall data force a rational second normal derivative."""
    s, y = sp.symbols("s y")
    a, b, mu = sp.symbols("a b mu", positive=True)
    J0 = a+3*b*s**2
    c = 1/(4*mu*a)
    h2 = -3*b*(1+c**2)*s/(5*J0**2)
    X = a*s+b*s**3+c*y+h2*y**2
    E = full_inverse_residual(X, s, y, 3*mu, mu)["E"]
    coefficients = [sp.factor(sp.diff(E,y,n).subs(y,0)/sp.factorial(n))
                    for n in range(3)]
    assert coefficients[0] == coefficients[1] == 0
    assert coefficients[2] != 0
    return {"X": str(X), "X_yy_at_wall": str(2*h2),
            "E_y_coefficients_0_1_2": [str(c) for c in coefficients],
            "scope": "wall and first y coefficient vanish; full closure fails at y**2"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path(
        "runs/fresh_symbolic_closure_20260928"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    reports = {
        "full_coordinate_formula_tests": verify_full_formula(),
        "sos_quartic_wall": sos_wall_certificates(),
        "normalized_saturated_wall": normalized_saturated_wall_certificates(),
        "cubic_forced_curvature": cubic_forced_curvature(),
        "algebraic_mapping_tail_gate": {
            "assumptions":"f=(1+t**2)**(-p), G=t+k*(sqrt(1+t**2)-1), 0<p<1/2, abs(k)<1, A,W>0; all six height tangents m,A,W,c,k,p permitted; global transverse family",
            "first_layer_operator_power_support":["2*p-1-2*n", "-2*n", "1-2*p-2*n", "-1-2*n", "-2*p-2*n"],
            "n":"nonnegative integer; p derivative adds logs but no exponents",
            "forcing_leading_power":"-4*p",
            "necessary_p":["1/6", "1/4"],
            "previous_loose_gate":"p=1/3 is allowed only by the old integer-step superset, and fails the tighter even-step support gate",
            "p_derivative_condition":"(2*lambda-mu/(1-2*p))*p_prime=0; at exception do not apply the fixed-p primitive-independence proof"
        },
        "exact_positive_wall_slope": {
            "assumptions":"mu>0, lambda>mu, J(s)>0 smooth across s=0",
            "formula":"K(s)=J(s)**2/mu * Integral(t**(2*(lambda/mu-1)-1)/J(s*t)**3, (t,0,1))",
            "conclusions":"K>0; exact wall first-layer identity; unique solution bounded at s=0",
            "scope":"Not a full two-dimensional closure or a global Biot-Savart certificate"
        },
        "certification_scope": "Exact local differential identities only; centered global BS remains an independent obligation.",
    }
    (args.out/"symbolic_certificates.json").write_text(
        json.dumps(reports, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps({"output": str(args.out/"symbolic_certificates.json"),
                      "formula_checks": reports["full_coordinate_formula_tests"],
                      "wall_branches": [{k:r[k] for k in
                        ("lambda_over_mu", "K_degree", "real_branch")}
                        for r in reports["sos_quartic_wall"]]},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
