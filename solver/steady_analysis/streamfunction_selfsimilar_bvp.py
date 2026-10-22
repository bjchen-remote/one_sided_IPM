"""Exact log-radius closure and a numerical well-posedness audit.

Derive the angular-profile system in the drift-centred outer patch and
demonstrate why it must be solved as a two-point boundary-value problem,
not by Cauchy marching in radius.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline

from fixed_point import fd_matrix
from streamfunction_sources import load_native


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/streamfunction_selfsimilar_closure_20260929"


def symbolic_checks() -> dict:
    t, theta = sp.symbols("t theta", real=True)
    lam, mu = sp.symbols("lambda mu", positive=True)
    alpha = mu / lam
    power = 1 + alpha
    A = sp.Function("A")(t, theta)
    B = sp.Function("B")(t, theta)
    radius = sp.exp(t)
    density = sp.exp(alpha * t) * A
    psi = sp.exp(power * t) * B

    density_r = sp.exp(-t) * sp.diff(density, t)
    density_theta_over_r = sp.exp(-t) * sp.diff(density, theta)
    velocity_r = -sp.exp(-t) * sp.diff(psi, theta)
    velocity_theta = sp.exp(-t) * sp.diff(psi, t)
    transport = ((lam * radius + velocity_r) * density_r
                 + velocity_theta * density_theta_over_r - mu * density)
    epsilon = sp.exp((alpha - 1) * t)
    reduced_transport = (lam * sp.diff(A, t)
        + epsilon * (-sp.diff(B, theta) * (sp.diff(A, t) + alpha * A)
                     + (sp.diff(B, t) + power * B) * sp.diff(A, theta)))
    transport_check = sp.simplify(sp.exp(-alpha * t) * transport
                                  - reduced_transport)
    assert transport_check == 0

    laplacian = sp.exp(-2 * t) * (sp.diff(psi, t, 2)
                                  + sp.diff(psi, theta, 2))
    density_x = sp.exp(-t) * (sp.cos(theta) * sp.diff(density, t)
                              - sp.sin(theta) * sp.diff(density, theta))
    reduced_poisson = (
        -(sp.diff(B, t, 2) + 2 * power * sp.diff(B, t)
          + power**2 * B + sp.diff(B, theta, 2))
        - sp.cos(theta) * (sp.diff(A, t) + alpha * A)
        + sp.sin(theta) * sp.diff(A, theta))
    poisson_check = sp.simplify(
        sp.exp((1 - alpha) * t) * (-laplacian - density_x)
        - reduced_poisson)
    assert poisson_check == 0

    C = sp.symbols("C", real=True)
    At = sp.symbols("A_t", real=True)
    Aval, Atheta, Bval, Btheta = sp.symbols(
        "A A_theta B B_theta", real=True)
    algebraic_transport = (
        lam * At
        + epsilon * (-Btheta * (At + alpha * Aval)
                     + (C + power * Bval) * Atheta))
    solved = sp.solve(algebraic_transport, At)[0]
    expected = (
        epsilon
        * (alpha * Btheta * Aval - (C + power * Bval) * Atheta)
        / (lam - epsilon * Btheta))
    assert sp.simplify(solved - expected) == 0

    return dict(
        transportIdentity=str(transport_check),
        poissonIdentity=str(poisson_check),
        alpha="mu/lambda",
        power="1+mu/lambda",
        epsilon="exp(-(1-alpha)t)",
        firstOrderSystem=dict(
            B_t="C",
            A_t=("epsilon[alpha B_theta A-(C+pB)A_theta]"
                 "/[lambda-epsilon B_theta]"),
            C_t=("-2pC-p^2B-B_theta_theta"
                 "-cos(theta)(A_t+alpha A)+sin(theta)A_theta")),
        characteristicDenominator="D=lambda-epsilon B_theta=V_r/r",
        zeroDenominatorCompatibility=(
            "alpha B_theta A-(C+pB)A_theta=0"))


def source_splines(profile: dict) -> dict:
    return {
        key: RectBivariateSpline(profile["y"], profile["x"], profile[key])
        for key in ("Q", "psi", "u1", "u2")
    }


def characteristic_scan(profile: dict, splines: dict) -> list[dict]:
    lam, _, v = map(float, profile["rates"][:3])
    theta = (np.arange(801) + 0.5) * np.pi / 801
    rows = []
    radii = (0.0625, 0.125, 0.25, 0.5, 1, 2, 4,
             8, 16, 64, 256, 1024, 4096)
    for native_radius in radii:
        radius = native_radius / profile["ell"]
        x = -v / lam + radius * np.cos(theta)
        y = radius * np.sin(theta)
        u1 = splines["u1"].ev(y, x)
        u2 = splines["u2"].ev(y, x)
        radial_velocity = u1 * np.cos(theta) + u2 * np.sin(theta)
        denominator = lam + radial_velocity / radius
        rows.append(dict(
            nativeRadius=float(native_radius),
            minimum=float(denominator.min()),
            maximum=float(denominator.max()),
            negativeFraction=float(np.mean(denominator < 0))))
    return rows


def cauchy_pilot(profile: dict, splines: dict) -> list[dict]:
    lam, mu, v = map(float, profile["rates"][:3])
    alpha = mu / lam
    power = 1 + alpha
    kappa = 1 - alpha
    start_native, target_native = 1.0, 16.0
    records = []

    for count in (17, 25, 33, 41):
        theta = np.linspace(0, np.pi, count)
        D1 = fd_matrix(theta, 1, min(7, count)).toarray()
        D2 = fd_matrix(theta, 2, min(7, count)).toarray()

        def observed(native_radius: float):
            radius = native_radius / profile["ell"]
            x = -v / lam + radius * np.cos(theta)
            y = radius * np.sin(theta)
            Q = splines["Q"].ev(y, x)
            psi = splines["psi"].ev(y, x)
            u1 = splines["u1"].ev(y, x)
            u2 = splines["u2"].ev(y, x)
            psi_r = np.cos(theta) * u2 - np.sin(theta) * u1
            A = (Q + v / mu) / radius**alpha
            B = psi / radius**power
            C = (radius * psi_r - power * psi) / radius**power
            B[[0, -1]] = 0
            C[[0, -1]] = 0
            return A, B, C

        A0, B0, C0 = observed(start_native)
        initial = np.r_[A0, B0[1:-1], C0[1:-1]]

        def unpack(state):
            A = state[:count]
            B = np.zeros(count)
            C = np.zeros(count)
            B[1:-1] = state[count:2 * count - 2]
            C[1:-1] = state[2 * count - 2:]
            return A, B, C

        def rhs(t, state):
            A, B, C = unpack(state)
            Btheta = D1 @ B
            Atheta = D1 @ A
            epsilon = np.exp(-kappa * t)
            denominator = lam - epsilon * Btheta
            At = (epsilon
                  * (alpha * Btheta * A - (C + power * B) * Atheta)
                  / denominator)
            Ct = (-2 * power * C - power**2 * B - D2 @ B
                  - np.cos(theta) * (At + alpha * A)
                  + np.sin(theta) * Atheta)
            return np.r_[At, C[1:-1], Ct[1:-1]]

        t0 = np.log(start_native / profile["ell"])
        t1 = np.log(target_native / profile["ell"])
        solution = solve_ivp(
            rhs, (t0, t1), initial, method="BDF",
            rtol=1e-6, atol=1e-8, max_step=0.02)
        A, B, _ = unpack(solution.y[:, -1])
        epsilon = np.exp(-kappa * solution.t[-1])
        denominator = lam - epsilon * (D1 @ B)
        records.append(dict(
            angularNodes=count,
            succeeded=bool(solution.success),
            requestedNativeRadius=target_native,
            reachedNativeRadius=float(
                np.exp(solution.t[-1]) * profile["ell"]),
            minimumDenominator=float(denominator.min()),
            maximumAbsA=float(np.max(abs(A))),
            maximumAbsB=float(np.max(abs(B))),
            functionEvaluations=int(solution.nfev),
            message=solution.message))
    return records


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    profile = load_native(20)
    splines = source_splines(profile)
    checks = symbolic_checks()
    scan = characteristic_scan(profile, splines)
    pilot = cauchy_pilot(profile, splines)
    result = dict(
        source=profile["source"],
        rates=dict(lambda_=float(profile["rates"][0]),
                   mu=float(profile["rates"][1]),
                   v=float(profile["rates"][2])),
        symbolicChecks=checks,
        characteristicScan=scan,
        cauchyPilot=pilot,
        interpretation=(
            "The transformed PDE is locally closed where D is nonzero. "
            "The global construction is a two-point boundary-value problem; "
            "radial Cauchy marching is numerically ill posed."))
    (OUT / "bvp_certificate.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(dict(
        symbolicChecks=checks,
        firstPositiveFullCircle=next(
            row for row in scan if row["minimum"] > 0),
        cauchyPilot=pilot,
        output=str(OUT / "bvp_certificate.json")), indent=2))


if __name__ == "__main__":
    main()
