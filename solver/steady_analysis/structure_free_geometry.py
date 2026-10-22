"""Exact full-PDE reduction for curved square-root density contours.

No discrete profile correction is performed. The eight polynomial coefficients
are simultaneous ODE constraints, not eight terms of a wall Taylor truncation.
"""
import json
from pathlib import Path

import sympy as sp


def residual_polynomial():
    s, t, k, lam, mu = sp.symbols("s t k lambda mu", real=True)
    F, A, C = (sp.Function(name)(s) for name in ("F", "A", "C"))
    Y = C * (t**2 + 2 * t)
    X = F + A * t + k * Y
    J = sp.expand(sp.diff(X, s) * sp.diff(Y, t) - sp.diff(X, t) * sp.diff(Y, s))
    N = sp.integrate((2 * lam - mu) * J - mu * s * sp.diff(J, s), t)
    U = [sp.expand(mu * s * sp.diff(r, s) * J + N * sp.diff(r, t)) for r in (X, Y)]
    L = sp.expand(U[0] * sp.diff(X, t) + U[1] * sp.diff(Y, t))
    M = sp.expand(U[0] * sp.diff(X, s) + U[1] * sp.diff(Y, s))
    R = sp.Poly(sp.expand(
        J * sp.diff(L, s) - L * sp.diff(J, s)
        - J * sp.diff(M, t) + M * sp.diff(J, t)
        + sp.diff(Y, t) * J**2
    ), t)
    assert R.degree() == 7
    return s, t, k, lam, mu, F, A, C, X, Y, J, N, R


def jet_substitution(old, new, variable):
    return {sp.diff(old, variable, j): sp.diff(new, variable, j) for j in range(4)}


def checks_and_export():
    s, t, k, lam, mu, F, A, C, X, Y, J, N, R = residual_polynomial()
    v = sp.symbols("v", real=True)
    phi = (C * (lam * F + v - mu * s * sp.diff(F, s)) * (t**2 + 2 * t)
           + mu * s * (A * sp.diff(C, s) - C * sp.diff(A, s)) * t**2
           + (lam * A * C - 2 * mu * s * C * sp.diff(A, s)
              + mu * s * A * sp.diff(C, s)) * t**3 / 3)
    transport = sp.diff(phi, t) - lam * (X * sp.diff(Y, t) - Y * sp.diff(X, t)) - v * sp.diff(Y, t) + mu * s * J
    assert sp.expand(transport) == 0
    D, H = sp.Function("D")(s), sp.Function("H")(s)
    change = jet_substitution(A, D * sp.sqrt(C), s)
    top = sp.factor(sp.simplify(R.coeff_monomial(t**7).subs(change, simultaneous=True)))
    expected_top = 16 * mu * (1 + k**2) * C**5 / 3 * (
        -s * sp.diff(D, s) * sp.diff(D, s, 3)
        + s * sp.diff(D, s, 2)**2 - sp.diff(D, s) * sp.diff(D, s, 2)
    )
    assert sp.simplify(top - expected_top) == 0

    d = sp.symbols("d", nonzero=True, real=True)
    special_A = d * s * sp.sqrt(C)
    change = jet_substitution(A, special_A, s)
    change.update(jet_substitution(F, H + special_A, s))
    sixth = sp.factor(sp.simplify(R.coeff_monomial(t**6).subs(change, simultaneous=True)))
    bracket = ((1 + k**2) * (3 * mu * s * sp.diff(H, s, 3)
               + (4 * mu - 2 * lam) * sp.diff(H, s, 2))
               + 2 * d**2 * k * (lam - 2 * mu))
    assert sp.simplify(sixth + 8 * d * C**sp.Rational(9, 2) * bracket / 3) == 0

    d0, d1, c0, c1, c2 = sp.symbols("d0 d1 c0 c1 c2", real=True)
    special_C = c0 + c1 * s + c2 * s**2
    special_A = (d0 + d1 * s) * sp.sqrt(special_C)
    change = jet_substitution(A, special_A, s)
    change.update(jet_substitution(C, special_C, s))
    change.update(jet_substitution(F, H + special_A, s))
    sixth_n2 = sp.factor(sp.simplify(
        R.coeff_monomial(t**6).subs(lam, 2 * mu).subs(change, simultaneous=True)))
    assert sp.simplify(sixth_n2 + 8 * d1 * mu * s * (1 + k**2)
                       * special_C**sp.Rational(9, 2) * sp.diff(H, s, 3)) == 0

    # Exact wall test of the n=2 quadratic surface, after the wall resonance
    # 2 d1 c0 = d0 c1 and a harmless spatial scaling have been imposed.
    a, b, h0, h1, z = sp.symbols("a b h0 h1 z", real=True)
    C2 = 1 + 2 * a * s + b * s**2
    D2 = d * (1 + a * s)
    F2 = h0 * s + h1 * s**2 / 2 + D2 * sp.sqrt(C2)
    J2 = sp.diff(F2, s)
    K2 = k + D2 / (2 * sp.sqrt(C2))
    wall = sp.together(mu * s * sp.diff(K2, s) * J2 + 2 * mu * K2 * J2
                       - 2 * mu * s * K2 * sp.diff(J2, s) - 1)
    numerator, denominator = wall.as_numer_denom()
    numerator = sp.expand(numerator).subs(sp.sqrt(C2), z)
    reduced = sp.rem(sp.Poly(numerator, z), sp.Poly(z**2 - C2, z)).as_expr()
    rational_part = sp.factor(reduced.coeff(z, 0))
    radical_part = sp.factor(reduced.coeff(z, 1))
    q = sp.Poly(radical_part / (d * mu), s)
    assert sp.expand(q.coeff_monomial(1) - 2 * (h0 + 4 * a * k)) == 0
    assert sp.expand(q.coeff_monomial(s**2).subs(h0, -4 * a * k) - 8 * a * k * (b - a**2)) == 0
    assert sp.expand(q.coeff_monomial(s**3).subs({k: 0, h0: 0}) - h1 * (a**2 - b)) == 0
    p = sp.Poly(rational_part.subs({k: 0, h0: 0, h1: 0}), s)
    assert sp.expand(p.coeff_monomial(1) - (4 * a * d**2 * mu - 2)) == 0
    assert sp.expand(p.coeff_monomial(s**2).subs(mu, 1 / (2 * a * d**2)) - 2 * (a**2 - b)) == 0

    output = {
        "meaning": "Full PDE is equivalent to all eight coefficients vanishing where J != 0; no wall truncation.",
        "X": str(X), "Y": str(Y), "Jacobian": str(J), "transverse_numerator": str(N),
        "transport_streamfunction": str(phi), "transport_identity": "exact zero",
        "full_ODE_coefficients_t0_through_t7": [str(sp.factor(R.coeff_monomial(t**j))) for j in range(8)],
        "highest_coefficient_after_A_equals_D_sqrtC": str(top),
        "D_equals_d_s_sixth_coefficient": str(sixth),
        "n2_affine_D_quadratic_C_sixth_coefficient": str(sixth_n2),
        "n2_wall_denominator": str(denominator),
        "n2_wall_rational_part": str(rational_part),
        "n2_wall_radical_part": str(radical_part),
        "claims": {
            "finite_ODE_reduction": "exact",
            "general_variable_C_family_has_full_PDE_solution": "not established",
            "n2_nonconstant_affine_D_positive_definite_quadratic_C": "excluded within this ansatz",
            "identified_with_discrete_branch": False,
        },
    }
    folder = Path(__file__).resolve().parents[2] / "runs" / "structure_free_alpha_20260928"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "free_geometry_checks.json").write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    print("Exact full-PDE coefficient reduction and all symbolic assertions passed.")
    print(folder / "free_geometry_checks.json")


if __name__ == "__main__":
    checks_and_export()
