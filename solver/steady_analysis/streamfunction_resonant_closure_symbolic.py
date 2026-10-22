"""Exact algebra for the N=2, m=1 centered Mellin resonance."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    delta, radius = sp.symbols("delta radius", real=True, positive=True)
    time, reference_time = sp.symbols("time reference_time", real=True)
    alpha = sp.Rational(1, 2) + delta
    beta2 = sp.simplify(2 * alpha - 1)
    kappa2 = sp.simplify(beta2 + 1)
    denominator = sp.factor(1 - kappa2**2)

    # Centering subtracts the m=1 harmonic r sin(theta).  Writing r=exp(t)
    # isolates the radial factor and makes the resonance limit elementary.
    centered = sp.simplify(radius * (
        sp.exp(2 * delta * time)
        - sp.exp(2 * delta * reference_time)) / denominator)
    resonant_limit = sp.simplify(sp.limit(centered, delta, 0))
    expansion = sp.series(centered / radius, delta, 0, 3).removeO()

    assert beta2 == 2 * delta
    assert kappa2 == 1 + 2 * delta
    assert sp.simplify(denominator + 4 * delta * (1 + delta)) == 0
    assert sp.simplify(
        resonant_limit + radius * (time - reference_time) / 2) == 0
    assert sp.simplify(
        expansion.coeff(delta, 0)
        + (time - reference_time) / 2) == 0
    assert sp.simplify(
        expansion.coeff(delta, 1)
        - ((time - reference_time)
           - (time**2 - reference_time**2)) / 2) == 0

    certificate = {
        "alpha": str(alpha),
        "beta2": str(beta2),
        "kappa2": str(kappa2),
        "m1Denominator": str(denominator),
        "centeredRadialFactor": str(centered),
        "resonantLimit": str(resonant_limit),
        "deltaSeriesThroughSecondOrder": str(expansion),
        "referenceScale": "reference_time = log(r_*) is fixed by Green centering",
        "leadingDetuningError": "delta*(time-reference_time)**2",
        "finiteBoxBalance": "delta = O((log(L/r_*))**-2)",
        "allAssertionsPassed": True,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as stream:
        json.dump(certificate, stream, indent=2)
        stream.write("\n")
    print(json.dumps(certificate, indent=2))


if __name__ == "__main__":
    main()
