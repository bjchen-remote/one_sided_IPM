"""Exact checks of algebraic parametrizations, without a Newton correction.

Run with sympy available. The output does not claim to identify the numerical
profile; its exclusions have the scope stated in COMPLEX.md.
"""

import json
from pathlib import Path

import sympy as sp


def cubic_parametrization():
    s, t, k, b = sp.symbols("s t k b", real=True)
    mu = 1 / (4 * k)
    h = 1 + k**2
    y = t * (1 + s**2) + k * s * t**2 + b * t**3
    x = s + s**3 / 3 + k * y - h * s * t**2 / 5 - 8 * k * h * t**3 / 105
    jac = sp.expand(sp.diff(x, s) * sp.diff(y, t) - sp.diff(x, t) * sp.diff(y, s))
    transverse_numerator = mu * sp.integrate(5 * jac - s * sp.diff(jac, s), t)
    wxnum = sp.expand(mu * s * sp.diff(x, s) * jac + transverse_numerator * sp.diff(x, t))
    wynum = sp.expand(mu * s * sp.diff(y, s) * jac + transverse_numerator * sp.diff(y, t))
    oneform_t = sp.expand(wxnum * sp.diff(x, t) + wynum * sp.diff(y, t))
    oneform_s = sp.expand(wxnum * sp.diff(x, s) + wynum * sp.diff(y, s))
    numerator = sp.expand(
        sp.diff(oneform_t, s) * jac - oneform_t * sp.diff(jac, s)
        - sp.diff(oneform_s, t) * jac + oneform_s * sp.diff(jac, t)
        + sp.diff(y, t) * jac**2
    )
    polynomial = sp.Poly(numerator, t)
    assert all(polynomial.coeff_monomial(t**j) == 0 for j in range(3))
    cubic = sp.factor(polynomial.coeff_monomial(t**3))
    expected = 2 * s * h * (1 + s**2)**2 / (105 * k) * (
        (-420 * b + 127 * k**2 + 196) * s**2
        - 420 * b + 106 * k**2 - 140
    )
    assert sp.expand(cubic - expected) == 0
    coefficient_difference = sp.expand((127 * k**2 + 196) - (106 * k**2 - 140))
    assert coefficient_difference == 21 * (k**2 + 16)
    # Even the purely complex value passing this first obstruction fails next.
    complex_next = sp.factor(polynomial.coeff_monomial(t**4).subs({k: 4 * sp.I, b: -sp.Rational(153, 35)}))
    assert complex_next == 36 * (s**2 + 1)**2 * (111 * s**2 - 1) / 245
    return {
        "X": str(x), "Y": str(y), "Jacobian": str(sp.factor(jac)),
        "curl_numerator_orders_0_1_2": "exact zero",
        "curl_numerator_order_3": str(cubic),
        "necessary_parameter_equation": "k**2 = -16",
        "complex_continuation_order_4": str(complex_next),
    }


def straight_levels():
    s, y, lam, mu, v = sp.symbols("s y lambda mu v")
    f, a = sp.Function("F")(s), sp.Function("a")(s)
    x = f + a * y
    jac = sp.diff(x, s)
    phi = (lam * f + v - mu * s * sp.diff(f, s)) * y - mu * s * sp.diff(a, s) * y**2 / 2
    dx = lambda expr: sp.diff(expr, s) / jac
    dy = lambda expr: sp.diff(expr, y) - a * sp.diff(expr, s) / jac
    numerator = sp.cancel((-dx(dx(phi)) - dy(dy(phi)) - 1 / jac) * jac**3)
    top = sp.factor(sp.Poly(numerator, y).coeff_monomial(y**3))
    theta_prime = sp.diff(a, s) / (1 + a**2)
    theta_log_derivative = sp.diff(theta_prime, s) / theta_prime
    angular_ode = sp.diff(s * theta_log_derivative, s) + 4 * s * theta_prime**2
    reconstructed = mu * (1 + a**2) * sp.diff(a, s)**2 * angular_ode / 2
    assert sp.cancel(top - reconstructed) == 0
    c, p = sp.symbols("C p", real=True, nonzero=True)
    positive_s = sp.Symbol("positive_s", positive=True)
    theta_trial_prime = p * c * positive_s**(p - 1) / (1 + c**2 * positive_s**(2 * p))
    test = sp.diff(positive_s * sp.diff(theta_trial_prime, positive_s) / theta_trial_prime, positive_s) + 4 * positive_s * theta_trial_prime**2
    assert sp.simplify(test) == 0
    return {
        "highest_y_coefficient": str(top),
        "angular_ode": "(s theta''/theta')' + 4 s (theta')**2 = 0",
        "analytic_nonconstant_branch": "theta(s) = theta_0 + atan(C s**p), p positive integer",
        "branch_substitution": "exact zero",
    }


def asymmetric_abelian_wall():
    z, eps = sp.symbols("z epsilon", real=True)
    aa, bb = 1 + eps, 1 - eps
    divisor = (z**2 + 1) * (aa * z**2 + bb)
    rational_differential = (aa * z**4 + bb) * (z**2 + 1) / z**4
    logarithmic_derivative = sp.Rational(2, 3) * z * (1 / (z**2 + 1) - aa / (aa * z**2 + bb))
    primitive_trial = (z**2 - 1) * divisor / (3 * z**3)
    obstruction = sp.factor(sp.cancel(
        (sp.diff(primitive_trial, z) + logarithmic_derivative * primitive_trial - rational_differential) * z**4
    ))
    expected = -2 * z**2 * ((8 * eps + 3) * z**2 - 8 * eps + 3) / 9
    assert sp.expand(obstruction - expected) == 0
    return {
        "curve": "w**3 = (z**2 + 1) / ((1 + epsilon)*z**2 + 1 - epsilon)",
        "genus_for_0_lt_abs_epsilon_lt_1": 2,
        "primitive_rational_factor_obstruction": str(obstruction),
        "scope": "Parent agent's C=2 asymmetric wall family only; no two-dimensional solution is asserted.",
    }


if __name__ == "__main__":
    result = {"cubic_parametrization": cubic_parametrization(), "straight_levels": straight_levels(),
              "asymmetric_abelian_wall": asymmetric_abelian_wall()}
    destination = Path(__file__).resolve().parents[2] / "runs" / "structure_20260928" / "complex_checks.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
