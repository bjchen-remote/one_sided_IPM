"""Exact Fourier/Mellin tail algebra for the half-plane IPM profile.

No saved discrete profile is loaded. These are formal asymptotic identities,
not existence, convergence, or selection results. Run with sympy available.
"""
from collections import defaultdict
import json
from pathlib import Path

import sympy as sp


def add_sine(target, mode, amplitude):
    if mode > 0:
        target[mode] += amplitude
    elif mode < 0:
        target[-mode] -= amplitude


def horizontal_derivative(coefficients, beta):
    """r^beta sum f_n sin(n theta) -> r^(beta-1) sum a_n sin(n theta)."""
    result = defaultdict(lambda: sp.S.Zero)
    for mode, amplitude in coefficients.items():
        add_sine(result, mode + 1, (beta - mode) * amplitude / 2)
        add_sine(result, mode - 1, (beta + mode) * amplitude / 2)
    return {n: sp.factor(a) for n, a in result.items() if a != 0}


def poisson_sine(coefficients, beta):
    source = horizontal_derivative(coefficients, beta)
    result = {}
    for mode, amplitude in source.items():
        denominator = mode**2 - (beta + 1)**2
        if denominator == 0:
            raise ValueError(f"Poisson resonance: beta={beta}, mode={mode}, amplitude={amplitude}")
        result[mode] = sp.factor(amplitude / denominator)
    return result


def advection_sine(stream, density, beta_i, beta_j):
    """Angular coefficients of u_i dot grad P_j."""
    result = defaultdict(lambda: sp.S.Zero)
    for m, g in stream.items():
        for n, f in density.items():
            add_sine(result, n + m, -beta_j * m * g * f / 2)
            add_sine(result, n - m, -beta_j * m * g * f / 2)
            add_sine(result, m + n, (beta_i + 1) * n * g * f / 2)
            add_sine(result, m - n, (beta_i + 1) * n * g * f / 2)
    return result


def symbolic_identities():
    theta, beta = sp.symbols('theta beta', real=True)
    beta_i, beta_j = sp.symbols('beta_i beta_j', real=True)
    r = sp.symbols('r', positive=True)
    F = sp.Function('F')(theta)
    G = sp.Function('G')(theta)
    P = r**beta * F
    dx = sp.cos(theta) * sp.diff(P, r) - sp.sin(theta) / r * sp.diff(P, theta)
    A = beta * sp.cos(theta) * F - sp.sin(theta) * sp.diff(F, theta)
    assert sp.simplify(dx - r**(beta - 1) * A) == 0
    psi = r**(beta + 1) * G
    lap = sp.diff(psi, r, 2) + sp.diff(psi, r) / r + sp.diff(psi, theta, 2) / r**2
    assert sp.simplify(lap - r**(beta - 1) * ((beta + 1)**2 * G + sp.diff(G, theta, 2))) == 0
    P_j = r**beta_j * F
    psi_i = r**(beta_i + 1) * G
    advection = -sp.diff(psi_i, theta) / r * sp.diff(P_j, r) + sp.diff(psi_i, r) / r * sp.diff(P_j, theta)
    bracket = -beta_j * sp.diff(G, theta) * F + (beta_i + 1) * G * sp.diff(F, theta)
    assert sp.simplify(advection - r**(beta_i + beta_j - 1) * bracket) == 0
    # Integration by parts against sin(2 theta); endpoint terms vanish.
    strain_weight = (beta + 1) * sp.cos(theta) * sp.sin(2 * theta) + 2 * sp.sin(theta) * sp.cos(2 * theta)
    strain_expected = (beta - 1) / 2 * sp.sin(theta) + (beta + 3) / 2 * sp.sin(3 * theta)
    assert sp.trigsimp(sp.expand_trig(strain_weight - strain_expected)) == 0
    for k in range(1, 7):
        positive = k * (sp.cos(theta) * sp.sin(k * theta) + sp.sin(theta) * sp.cos(k * theta))
        negative = k * (-sp.cos(theta) * sp.sin(k * theta) + sp.sin(theta) * sp.cos(k * theta))
        assert sp.trigsimp(sp.expand_trig(positive - k * sp.sin((k + 1) * theta))) == 0
        assert sp.trigsimp(sp.expand_trig(negative + k * sp.sin((k - 1) * theta))) == 0
    # Fourier coefficient at beta=-1: no pole for any positive sine mode.
    gamma_values = {}
    for n in range(1, 7):
        value = 2**sp.S.One * sp.pi * (-sp.I)**n * sp.gamma(sp.Rational(n + 1, 2)) / sp.gamma(sp.Rational(n + 1, 2))
        assert sp.simplify(value - 2 * sp.pi * (-sp.I)**n) == 0
        gamma_values[str(n)] = str(sp.simplify(value))
    # The beta=-3, sin(2 theta) resonant stream correction.
    test_psi = r**(-2) * sp.log(r) * sp.sin(2 * theta)
    minus_lap = -(sp.diff(test_psi, r, 2) + sp.diff(test_psi, r) / r + sp.diff(test_psi, theta, 2) / r**2)
    assert sp.simplify(minus_lap - 4 * r**(-4) * sp.sin(2 * theta)) == 0
    return dict(
        horizontalDerivative='exact symbolic zero',
        polarPoisson='exact symbolic zero',
        polarAdvection='exact symbolic zero',
        strainMomentIdentity='exact symbolic zero',
        resonanceProjectionIdentitiesChecked=list(range(1, 7)),
        fourierCoefficientAtBetaMinus1=gamma_values,
        betaMinus1DirichletEigenvalues='n^2, n=1,2,...; strictly positive',
        resonantLogStreamMinusLaplacian='4*r**(-4)*sin(2*theta)',
    )


def formal_example():
    """An arbitrary exact angular seed, not a candidate or discrete fit."""
    alpha = sp.Rational(1, 3)
    lam = sp.S.One
    beta = {n: n * alpha - (n - 1) for n in range(1, 7)}
    density = {1: {1: sp.S.One}}
    stream = {}
    checks = []
    for n in range(1, 7):
        if n > 1:
            source = defaultdict(lambda: sp.S.Zero)
            for i in range(1, n):
                for mode, amplitude in advection_sine(stream[i], density[n - i], beta[i], beta[n - i]).items():
                    source[mode] += amplitude
            density[n] = {
                mode: sp.factor(amplitude / ((n - 1) * lam * (1 - alpha)))
                for mode, amplitude in source.items() if amplitude != 0
            }
        if n < 6:
            stream[n] = poisson_sine(density[n], beta[n])
        checks.append(dict(n=n, beta=str(beta[n]), densitySine1=str(density[n].get(1, 0)),
                           densitySineModes={str(k): str(v) for k, v in sorted(density[n].items())}))
    resonant_source = horizontal_derivative(density[6], beta[6])[2]
    assert resonant_source == -2 * density[6][1]
    assert resonant_source != 0
    return dict(
        purpose='Illustrate a genuine later resonance, not solve or fit the profile.',
        alpha='1/3', lambdaValue='1', leadingDensityAngularFunction='sin(theta)',
        coefficients=checks,
        firstResonantPoissonSource=dict(order=6, beta='-3', sineMode=2, amplitude=str(resonant_source)),
        necessaryLogStreamCoefficient=str(resonant_source / 4),
        warning='No wall normalization, positive Q_x, inner matching, or convergence is claimed.',
    )


def main():
    root = Path(__file__).resolve().parents[2]
    output = root / 'runs' / 'fourier_20260928' / 'tail_probe.json'
    result = dict(
        conventions=dict(fourier='integral exp(-i*k.x)*P(x) dx',
                         oddExtension='p(x,y)=sign(y)*P(x,abs(y))',
                         velocityMultiplier='(kx*ky/|k|^2, -kx^2/|k|^2)',
                         gammaCoefficient='2^(beta+2)*pi*(-i)^n*Gamma((n+beta+2)/2)/Gamma((n-beta)/2), n>=1'),
        symbolic=symbolic_identities(),
        arbitraryFormalSeed=formal_example(),
        selectionConclusion='Neither homogeneity nor the beta=-1 third term selects alpha=1/3.',
        purePowerCaveat='A constant effective drift also generates degrees beta_n-m; remove by an unknown horizontal shift or retain these terms.',
        claimsExcluded=['existence', 'convergence', 'global profile', 'alpha selection', 'discrete data validation'],
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(output=str(output), symbolic=result['symbolic'],
                         firstResonance=result['arbitraryFormalSeed']['firstResonantPoissonSource']), indent=2))


if __name__ == '__main__':
    main()
