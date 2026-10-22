"""Wall-only normal reconstruction; no interior fitting or steady-state claim."""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_jacobi
from construct_profile import Wall, load_profile, inverse_rows, flux_audit
from normal_lift_algebra import curl_coefficients


class NormalLift:
    def __init__(self, wall, count=5, degree=64, quadrature=64):
        self.wall = wall
        self.degree = degree
        self.domain = [min(-1.4+wall.d, -.15), .7+wall.d]
        self.coeff = []
        self.quadrature = quadrature
        for n in range(1, count+1):
            gamma = n*(2*wall.lam/wall.mu-1)-1
            nodes, weights = roots_jacobi(quadrature, 0, gamma-1)
            t = (nodes+1)/2
            weight = weights/2**gamma
            def generate(s):
                s = np.asarray(s)
                st = s[..., None]*t
                jt = wall.jacobian(st-wall.d)[0]
                source = self.source(n, st)
                j = wall.jacobian(s-wall.d)[0]
                return j**(n+1)/wall.mu*np.sum(weight*source/jt**(n+1), axis=-1)
            polynomial = Chebyshev.interpolate(generate, degree, domain=self.domain)
            # Suppress sub-roundoff coefficient tails before taking derivatives.
            polynomial = polynomial.trim(tol=2e-14*max(1., abs(polynomial.coef).max()))
            self.coeff.append(polynomial)

    def source(self, n, s):
        w = self.wall
        J, Jp, Jpp = w.jacobian(s-w.d)
        j = [J]+[a.deriv()(s) for a in self.coeff]
        js = [Jp]+[a.deriv(2)(s) for a in self.coeff]
        jss = [Jpp]+[a.deriv(3)(s) for a in self.coeff]
        k = [(i+1)*a(s) for i, a in enumerate(self.coeff)]
        ks = [(i+1)*a.deriv()(s) for i, a in enumerate(self.coeff)]
        ky = [(i+1)*(i+2)*a(s) for i, a in enumerate(self.coeff[1:])]
        return curl_coefficients(j, js, jss, k, ks, ky, s, w.lam, w.mu, n-1)[-1]/(n*J**3)

    def values(self, r):
        return np.array([a(np.asarray(r)+self.wall.d) for a in self.coeff])

    def evaluate(self, r, y, kind):
        """Broadcast heights against the last (density) axis."""
        r = np.asarray(r)
        h = np.asarray(y)[..., None]
        F = self.wall.F(r)
        a = self.values(r)
        if kind.startswith('taylor'):
            order = int(kind.removeprefix('taylor'))
            return F+sum(a[n-1]*h**n for n in range(1, order+1))
        if kind == 'pade22':
            det = a[1]**2-a[0]*a[2]
            d1 = (a[0]*a[3]-a[1]*a[2])/det
            d2 = (a[2]**2-a[1]*a[3])/det
            return F+(a[0]*h+(a[1]+d1*a[0])*h*h)/(1+d1*h+d2*h*h)
        raise ValueError(kind)

    def integral(self, r, y, kind):
        nodes, weights = leggauss(48)
        t = (nodes+1)/2
        # shape (height, quadrature, density)
        return np.asarray(y)[:, None]*np.sum(
            self.evaluate(r, np.asarray(y)[:, None]*t, kind)*weights[None, :, None]/2, axis=1)


def run(out, names, degrees):
    root = Path(__file__).resolve().parents[2]
    source = root/'runs/steady_20260928_tau20'
    old = root/'runs/construction_20260928'
    results = []
    for name in names:
        saved = np.load(old/f'{name}.npz')
        wall = Wall(saved['rates'], saved['parameters'])
        profile = load_profile(source, name)
        lifts = [NormalLift(wall, 5, n) for n in degrees]
        lift = lifts[-1]
        r = np.linspace(-1.25, .5, 161)
        y = np.linspace(0, 2, 161)
        truth = inverse_rows(profile, r, y)
        a = lift.values(r)
        refinement = []
        for i in range(1, len(lifts)):
            change = lifts[i].values(r)-lifts[i-1].values(r)
            refinement.append(dict(degrees=degrees[i-1:i+1],
                maxChangeByOrder=abs(change).max(axis=1).tolist()))
        det = a[1]**2-a[0]*a[2]
        d1 = (a[0]*a[3]-a[1]*a[2])/det
        d2 = (a[2]**2-a[1]*a[3])/det
        pade_a5 = -d1*a[3]-d2*a[2]
        positive_poles = []
        for i, (b, c) in enumerate(zip(d1, d2)):
            for z in np.roots([c, b, 1.]):
                if abs(z.imag)<1e-10 and z.real>0:
                    positive_poles.append([float(r[i]), float(z.real)])
        trials = []
        arrays = dict(r=r, y=y, observed=truth, coefficients=a, d1=d1, d2=d2)
        for kind in ['taylor1', 'taylor2', 'taylor3', 'taylor4', 'taylor5', 'pade22']:
            pred = lift.evaluate(r, y, kind)
            arrays[kind] = pred
            errors = []
            for lo, hi in [(0, .1), (.1, .25), (.25, .5), (.5, 1.), (1., 2.)]:
                mask = (y>lo)&(y<=hi)
                delta = pred[mask]-truth[mask]
                errors.append(dict(heightBand=[lo,hi], rms=float(np.sqrt(np.mean(delta**2))),
                    relativeDisplacementL2=float(np.linalg.norm(delta)/np.linalg.norm(truth[mask]-truth[0]))))
            audits = []
            for ymax in [.05, .1, .2, .5, 1., 2.]:
                yy = np.linspace(0, ymax, 81)
                xx = lift.evaluate(r, yy, kind)
                integral = lift.integral(r, yy, kind)
                audit = flux_audit(wall, r, yy, xx, integral)
                audit['heightMax'] = ymax
                if kind == 'pade22':
                    audit['hasSampledPoleInsideStrip'] = any(z[1]<=ymax for z in positive_poles)
                    audit['certificationStatus'] = 'rejected: determinant crossing destroys a uniform smooth wall lift; raw finite-difference diagnostics only'
                else:
                    audit['certificationStatus'] = 'local finite-difference diagnostic; invalid as a profile wherever J is nonpositive'
                audits.append(audit)
            trials.append(dict(kind=kind, errors=errors, fullPDEAudits=audits))
        independent_a1 = wall.first_jet(r)[0]
        independent_a2 = wall.second_jet(r)
        result = dict(name=name, rates=saved['rates'].tolist(),
            sourceKind='finite-time trajectory' if name=='tau20' else 'conditional finite-box root',
            soleFittedInput='previous wall-only Jacobian fit; saved rates fixed; no interior fit',
            coefficientMaxima=abs(a).max(axis=1).tolist(),
            interpolationDegrees=[len(p.coef)-1 for p in lift.coeff],
            refinement=refinement,
            comparisonWithIndependentLowOrderFormulas=[float(abs(a[0]-independent_a1).max()),float(abs(a[1]-independent_a2).max())],
            pade22=dict(fifthCoefficientMismatchLinf=float(abs(pade_a5-a[4]).max()),
                fifthCoefficientMismatchRelativeL2=float(np.linalg.norm(pade_a5-a[4])/np.linalg.norm(a[4])),
                positivePoles=positive_poles, minimumDeterminantMagnitude=float(abs(det).min())),
            trials=trials,
            caveats=['Finite-order reconstruction is not a convergent or global skeleton.',
                'Wall fitting and fixed rates introduce errors not measured by spectral refinement.',
                'Finite differences audit the full local PDE, not the centered infinite-domain Biot-Savart law.'])
        results.append(result)
        (out/'numeric.json').write_text(json.dumps(results, indent=2)+'\n')
        np.savez_compressed(out/f'{name}.npz', **arrays)
        print(name, 'coefficient maxima', result['coefficientMaxima'],
              'refinement', refinement, 'pade5 mismatch', result['pade22']['fifthCoefficientMismatchRelativeL2'], flush=True)
    return results


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--names', nargs='+', default=['L4_257x129','tau20','L8_257x129_continued'])
    parser.add_argument('--degrees', nargs='+', type=int, default=[48,64,80])
    args = parser.parse_args()
    out = Path(__file__).resolve().parents[2]/'runs/normal_lift_20260928'
    out.mkdir(exist_ok=True)
    run(out, args.names, args.degrees)
