"""Finite real-power reduction for a free transverse exponent beta.

The ansatz is tested algebraically; no full steady profile is claimed.
"""
import json
from itertools import combinations
from pathlib import Path

import sympy as sp


def run():
    s, beta, k, lam, mu, v = sp.symbols("s beta k lambda mu v", real=True)
    H, D, C = (sp.Function(name)(s) for name in ("H", "D", "C"))
    L = lam * H + v - mu * s * sp.diff(H, s)
    A = (lam * (1 - beta) * D - mu * s * sp.diff(D, s)) / (beta + 1)
    B = D * (lam * C - mu * s * sp.diff(C, s))

    # A pair (a,b) denotes exponent a*beta+b, not a floating-point power.
    def add(*terms):
        out = {}
        for term in terms:
            for power, coefficient in term.items():
                out[power] = out.get(power, 0) + coefficient
        return {power: sp.expand(coefficient) for power, coefficient in out.items()
                if sp.expand(coefficient) != 0}

    def scale(term, factor):
        return {power: coefficient * factor for power, coefficient in term.items()
                if coefficient * factor != 0}

    def multiply(*terms):
        out = {(0, 0): sp.Integer(1)}
        for term in terms:
            new = {}
            for (a, b), coefficient in out.items():
                for (aa, bb), other in term.items():
                    power = a + aa, b + bb
                    new[power] = new.get(power, 0) + coefficient * other
            out = new
        return out

    def derivative(term, variable):
        out = {}
        for (a, b), coefficient in term.items():
            power = a * beta + b
            if variable == "s":
                out[a, b] = out.get((a, b), 0) + sp.diff(coefficient, s)
            factor = sp.diff(C, s) if variable == "s" else 1
            out[a, b - 1] = out.get((a, b - 1), 0) + power * coefficient * factor
        return {power: sp.expand(coefficient) for power, coefficient in out.items()
                if sp.expand(coefficient) != 0}

    J = {(0, 0): sp.diff(H, s), (1, 0): sp.diff(D, s),
         (1, -1): beta * D * sp.diff(C, s)}
    K = {(0, 0): k, (1, -1): beta * D}
    phi = {(0, 1): L, (1, 1): A, (1, 0): B,
           (0, 0): -L * C - A * C**(beta + 1) - B * C**beta}
    X = {(0, 0): H - k * C, (0, 1): k, (1, 0): D}
    height = {(0, 1): sp.Integer(1), (0, 0): -C}
    transport = add(derivative(phi, "y"), scale(X, -lam),
                    scale(multiply(height, K), lam), {(0, 0): -v}, scale(J, mu * s))
    assert all(sp.factor(coefficient) == 0 for coefficient in transport.values())
    one_plus_K2 = add({(0, 0): sp.Integer(1)}, multiply(K, K))
    residual = add(
        multiply(one_plus_K2, J, derivative(derivative(phi, "s"), "s")),
        scale(multiply(K, J, J, derivative(derivative(phi, "s"), "y")), -2),
        multiply(J, J, J, derivative(derivative(phi, "y"), "y")),
        multiply(add(scale(multiply(K, J, derivative(K, "s")), 2),
                     scale(multiply(derivative(K, "y"), J, J), -1),
                     scale(multiply(one_plus_K2, derivative(J, "s")), -1)),
                 derivative(phi, "s")),
        multiply(J, J),
    )
    assert len(residual) == 17
    expected_top = mu * (1 + k**2) / (beta + 1) * (
        -s * sp.diff(D, s) * sp.diff(D, s, 3)
        + s * sp.diff(D, s, 2)**2 - sp.diff(D, s) * sp.diff(D, s, 2))
    assert sp.factor(residual[2, 1] - expected_top) == 0

    d0, d1 = sp.symbols("d0 d1", real=True)
    change = {sp.diff(D, s, j): sp.diff(d0 + d1 * s, s, j) for j in range(4)}
    affine = {power: sp.factor(coefficient.subs(change, simultaneous=True)
                               .subs(lam, mu / (1 - beta)))
              for power, coefficient in residual.items()}
    expected_first = d1 * (1 + k**2) * (
        (mu / (1 - beta) - 2 * mu) * sp.diff(H, s, 2)
        - mu * s * sp.diff(H, s, 3))
    expected_second = mu * beta * d0 * d1**3
    assert sp.factor(affine[1, 1] - expected_first) == 0
    assert sp.factor(affine[4, -1] - expected_second) == 0
    assert affine[2, 1] == 0 and affine[3, 0] == 0

    collision = sp.factor((affine[1, 1] + affine[4, -1]).subs(beta, sp.Rational(2, 3)))
    h0, h1, h3 = sp.symbols("h0 h1 h3")
    trial = h0 + h1 * s - d0 * d1**2 * s**2 / (3 * (1 + k**2)) + h3 * s**3
    trial_change = {sp.diff(H, s, j): sp.diff(trial, s, j) for j in range(4)}
    assert sp.factor(collision.subs(trial_change, simultaneous=True)) == 0

    # One further exact level, keeping C arbitrary. Group equal powers first.
    beta23_grouped = {}
    for (a, b), coefficient in residual.items():
        power = sp.Rational(2, 3) * a + b
        specialized = (coefficient.subs(change, simultaneous=True)
                       .subs(trial_change, simultaneous=True)
                       .subs({beta: sp.Rational(2, 3), lam: 3 * mu}))
        beta23_grouped[power] = beta23_grouped.get(power, 0) + specialized
    for power in (sp.Rational(7, 3), sp.Integer(2), sp.Rational(5, 3)):
        assert sp.factor(beta23_grouped[power]) == 0
    kk = 1 + k**2
    c_operator = (s * (d0 + d1 * s) * sp.diff(C, s, 3)
                  + (d1 * s - 2 * d0) * sp.diff(C, s, 2) - 4 * d1 * sp.diff(C, s))
    c_rhs = (d1 / (mu * kk) - 4 * d0 * d1**2 * k * (2 * d0 + d1 * s) / (9 * kk**2)
             - 4 * d1 * h1 * k / (3 * kk))
    next_coefficient = sp.factor(beta23_grouped[sp.Rational(4, 3)])
    assert sp.factor(next_coefficient + d1 * mu * kk * (c_operator - c_rhs)) == 0
    c0 = sp.symbols("c0", real=True)
    c_linear = -1 / (4 * mu * kk) + 4 * d0**2 * d1 * k / (27 * kk**2) + h1 * k / (3 * kk)
    c_quadratic = 2 * d0 * d1**2 * k / (27 * kk**2)
    c_particular = c0 + c_linear * s + c_quadratic * s**2
    c_change = {sp.diff(C, s, j): sp.diff(c_particular, s, j) for j in range(4)}
    assert sp.factor((c_operator - c_rhs).subs(c_change, simultaneous=True)) == 0
    c_homogeneous = d1**2 * s**3 - d0 * d1 * s**2 + d0**2 * s
    c_change_hom = {sp.diff(C, s, j): sp.diff(c_homogeneous, s, j) for j in range(4)}
    assert sp.factor(c_operator.subs(c_change_hom, simultaneous=True)) == 0

    # Global analytic screening uses reduction of order, not a numerical solve.
    z = sp.symbols("z", real=True)
    w1 = 1 - 2 * z + 3 * z**2
    assert sp.expand(z * (1 + z) * sp.diff(w1, z, 2)
                     + (z - 2) * sp.diff(w1, z) - 4 * w1) == 0
    reduction_integrand = z**2 / ((1 + z)**3 * w1**2)
    assert sp.limit((z + 1)**3 * reduction_integrand, z, -1) == sp.Rational(1, 36)
    assert w1.subs(z, -1) == 6
    # Consequently the independent w2 has leading term -1/[12*(z+1)^2].
    exponent = sp.symbols("exponent")
    euler_indicial = sp.expand(exponent * (exponent - 1) * (exponent - 2)
                              + exponent * (exponent - 1) - 4 * exponent)
    assert sp.factor(euler_indicial) == exponent * (exponent - 3) * (exponent + 1)
    d0_zero_rhs = 1 / (mu * kk) - 4 * h1 * k / (3 * kk)
    d0_zero_particular = -d0_zero_rhs * s / 4
    assert sp.factor(s**2 * sp.diff(d0_zero_particular, s, 3)
                     + s * sp.diff(d0_zero_particular, s, 2)
                     - 4 * sp.diff(d0_zero_particular, s) - d0_zero_rhs) == 0

    collision_betas = set()
    for (a, b), (aa, bb) in combinations(residual, 2):
        if a != aa:
            value = sp.Rational(bb - b, a - aa)
            if 0 < value < 1:
                collision_betas.add(value)
    target_collision_betas = set()
    for a, b in residual:
        if (a, b) in ((4, -1), (2, 1), (3, 0)) or a == 4:
            continue
        value = sp.Rational(b + 1, 4 - a)
        if 0 < value < 1:
            target_collision_betas.add(value)
    assert target_collision_betas == {sp.Rational(1, 4), sp.Rational(1, 3),
                                     sp.Rational(1, 2), sp.Rational(2, 3)}
    output = {
        "ansatz": "X=H(s)+k*y+D(s)*(C(s)+y)**beta, 0<beta<1, C>0",
        "power_pairs_a_beta_plus_b": sorted(residual),
        "number_of_formal_power_terms": len(residual),
        "possible_exponent_collision_betas": [str(x) for x in sorted(collision_betas)],
        "necessary_betas_when_affine_D_and_d0_d1_nonzero": [str(x) for x in sorted(target_collision_betas)],
        "corresponding_necessary_alphas": [str(x) for x in sorted(1 - x for x in target_collision_betas)],
        "phi_coefficients": {"L": str(L), "A": str(A), "B": str(B)},
        "transport_identity": "exact zero",
        "general_highest_coefficient_u_2beta_plus_1": str(sp.factor(expected_top)),
        "affine_D_farfield_relation": "lambda=mu/(1-beta)",
        "coefficient_u_beta_plus_1": str(sp.factor(expected_first)),
        "coefficient_u_4beta_minus_1": str(expected_second),
        "beta_2_over_3_combined_equation": str(collision),
        "beta_2_over_3_H_trial": str(trial),
        "trial_satisfies_this_combined_equation": "exact zero",
        "beta_2_over_3_next_nonzero_power": "4/3",
        "beta_2_over_3_C_ODE_operator": str(c_operator),
        "beta_2_over_3_C_ODE_rhs": str(c_rhs),
        "next_coefficient_equals": "-d1*mu*(1+k**2)*(C_ODE_operator-C_ODE_rhs)",
        "C_ODE_quadratic_particular_solution": str(c_particular),
        "C_ODE_particular_solution_residual": "exact zero",
        "C_ODE_nonconstant_homogeneous_solution": str(c_homogeneous),
        "C_ODE_homogeneous_solution_residual": "exact zero",
        "d0_nonzero_reduction_of_order_integrand": str(reduction_integrand),
        "second_w_solution_near_z_minus_1": "-1/(12*(z+1)**2) + lower singular order",
        "d0_zero_C_homogeneous_exponents": [-1, 0, 3],
        "global_analytic_beta_2_over_3_affine_D_nonconstant": {
            "status": "excluded within this single-power height-shift ansatz",
            "assumptions": ["s ranges over all real numbers", "C is globally real analytic and strictly positive",
                            "F'=d(H+D*C**beta)/ds is strictly positive", "regular wall branch at s=0",
                            "D=d0+d1*s with d1 != 0", "beta=2/3 and lambda=3*mu>0"],
            "reason": "The next-power ODE and global C positivity force a quadratic C for d0!=0, or constant C for d0=0; the wall slope then changes sign, contradicting the positive integral wall formula.",
        },
        "remaining_power_coefficients_solved": False,
        "identified_with_discrete_profile": False,
    }
    folder = Path(__file__).resolve().parents[2] / "runs" / "structure_free_alpha_20260928"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "free_power_checks.json").write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    print("Free-power reconstruction and leading-coefficient assertions passed; remaining equations are unsolved.")


if __name__ == "__main__":
    run()
