"""Independent whole-space Fourier / half-space wall-strain identities.

Convention: fhat(k,l) = int exp(-i*(k*x+l*y))*f(x,y) dx dy.
For -Delta psi=q, u=(-psi_y,psi_x), odd-y continuation of q gives
    (u1_x)hat = k*l/(k*k+l*l) * qhat.
At the wall, equivalently, with A(k)=int_0^infty exp(-k*y) qhat_x(k,y)dy,
    u1_x(0,0) = (1/pi) int_0^infty k*Im(A(k))dk
              = -(2/pi) int_{y>0} x*y/(x*x+y*y)**2 q(x,y) dxdy.

The saved profiles are OPTIONAL ancillary finite-box source diagnostics.
Their exterior source is unknown; none of these checks verifies a global root.
No saved state is modified. Set BLAS threads to two before invoking Python.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.interpolate import RectBivariateSpline
from scipy.special import erfcx


def gauss_interval(a, b, n):
    z, w = leggauss(n)
    return (a+b)/2+(b-a)*z/2, w*(b-a)/2


def cell_gauss(grid, n):
    z, w = leggauss(n)
    widths = np.diff(grid)
    points = ((grid[1:]+grid[:-1])[:, None]/2
              + widths[:, None]*z[None, :]/2)
    weights = widths[:, None]*w[None, :]/2
    return points.ravel(), weights.ravel()


def polar_strain(qdiff, lx, ly, order=80, shell_edges=None):
    """Integrate paired qdiff=q(x,y)-q(-x,y), never sampling singular origin."""
    if shell_edges is None:
        shell_edges = [0, .25, .5, 1, 2, 4, 6, 8, np.hypot(lx, ly)]
    edges = sorted(set(float(r) for r in shell_edges if 0 <= r <= np.hypot(lx, ly)))
    if edges[-1] < np.hypot(lx, ly):
        edges.append(float(np.hypot(lx, ly)))
    cuts = [0, np.arctan2(ly, lx), np.pi/2]
    for r in edges[1:]:
        if lx < r:
            cuts.append(np.arccos(lx/r))
        if ly < r:
            cuts.append(np.arcsin(ly/r))
    cuts = sorted(set(cuts))
    contributions = np.zeros(len(edges)-1)
    z, w = leggauss(order)
    for ta, tb in zip(cuts[:-1], cuts[1:]):
        theta, wt = gauss_interval(ta, tb, order)
        c, s = np.cos(theta), np.sin(theta)
        rmax = np.minimum(lx/c, ly/s)
        for j, (ra, rb) in enumerate(zip(edges[:-1], edges[1:])):
            lo, hi = np.minimum(ra, rmax), np.minimum(rb, rmax)
            radius = (lo[:, None]+hi[:, None])/2+(hi-lo)[:, None]*z/2
            wr = (hi-lo)[:, None]*w/2
            xx, yy = radius*c[:, None], radius*s[:, None]
            values = qdiff(xx.ravel(), yy.ravel()).reshape(radius.shape)
            kernel = -(2/np.pi)*c[:, None]*s[:, None]/radius
            contributions[j] += np.sum(wt[:, None]*wr*kernel*values)
    return dict(total=float(contributions.sum()), shells=[
        dict(rmin=a, rmax=b, contribution=float(v))
        for a, b, v in zip(edges[:-1], edges[1:], contributions)])


def fourier_laplace(qdiff, xgrid, ygrid, cell_order=4,
                    cutoffs=(8., 16., 32., 64.)):
    """Continuous x-transform then Laplace transform; no FFT periodic images.

    Composite Gauss on every interpolation cell resolves x/y variation.
    Composite Gauss in k has panels no wider than pi/Lx. This diagnostic
    reports cutoff dependence, not an extrapolated infinite-frequency value.
    """
    xx, wx = cell_gauss(xgrid, cell_order)
    yy, wy = cell_gauss(ygrid, cell_order)
    xmesh, ymesh = np.meshgrid(xx, yy)
    source = qdiff(xmesh.ravel(), ymesh.ravel()).reshape(ymesh.shape)
    weighted_source = source*wy[:, None]*wx[None, :]
    previous = 0.
    cumulative = 0.
    result = []
    for cutoff in cutoffs:
        count = int(np.ceil((cutoff-previous)*xgrid[-1]/np.pi))
        kk, wk = cell_gauss(np.linspace(previous, cutoff, count+1), 8)
        integral = 0.
        for first in range(0, len(kk), 64):
            k = kk[first:first+64]
            # Im A(k)=-int sin(k*x)exp(-k*y)*qdiff(x,y) dxdy.
            imag_a = -np.sum((np.exp(-k[:, None]*yy)@weighted_source)
                            *np.sin(k[:, None]*xx), axis=1)
            integral += float(np.dot(wk[first:first+64], k*imag_a/np.pi))
        cumulative += integral
        result.append(dict(kmax=cutoff, strain=cumulative,
                           last_band_contribution=integral))
        previous = cutoff
    return result


def manufactured_checks():
    """q=x*y*exp(-x*x-y*y), odd in y; all-plane exact strain=-1/8."""
    def laplace_first_moment(k):
        return .5-np.sqrt(np.pi)*k*erfcx(k/2)/4

    def spectral_integrand(k):
        return -k*k*np.exp(-k*k/4)*laplace_first_moment(k)/(2*np.sqrt(np.pi))

    # Both Fourier coordinates, using polar integration of the multiplier.
    angle, ea = quad(lambda t: np.cos(t)**2*np.sin(t)**2, 0, 2*np.pi,
                     epsabs=1e-13, epsrel=1e-12)
    radius, er = quad(lambda r: r**3*np.exp(-r*r/4), 0, np.inf,
                      epsabs=1e-13, epsrel=1e-12)
    whole_space = -angle*radius/(16*np.pi)
    horizontal, eh = quad(spectral_integrand, 0, np.inf, epsabs=1e-13, epsrel=1e-12)
    paired = lambda x, y: 2*x*y*np.exp(-x*x-y*y)
    physical = polar_strain(paired, 6, 6, order=64)['total']
    direct_transforms = fourier_laplace(paired, np.linspace(0, 6, 49),
                                      np.linspace(0, 6, 49), cell_order=6,
                                      cutoffs=(4., 8., 16., 32.))
    exact = -.125
    errors = dict(whole_space=abs(whole_space-exact),
                  horizontal_fourier_laplace=abs(horizontal-exact),
                  physical_kernel=abs(physical-exact),
                  numerical_source_transforms=abs(direct_transforms[-1]['strain']-exact))
    if max(errors.values()) > 1e-8:
        raise AssertionError(f'Manufactured identity check failed: {errors}')
    return dict(source='q(x,y)=x*y*exp(-x*x-y*y), y>=0',
                odd_extension_fourier='-pi*k*ell/4*exp(-(k*k+ell*ell)/4)',
                exact_strain=exact, whole_space_2d_fourier=whole_space,
                horizontal_fourier_laplace=horizontal,
                physical_kernel_box_L6=physical,
                numerical_source_transforms=direct_transforms,
                absolute_errors=errors,
                quadrature_reported_error=dict(angle=ea, radius=er, horizontal=eh))


def profile_probe(path):
    with np.load(path) as d:
        x, y, q = d['x'].copy(), d['y'].copy(), d['qx'].copy()
        rates = d['rates'].copy()
    if not np.allclose(x, -x[::-1]):
        raise ValueError('Only symmetric x boxes are supported')
    spline = RectBivariateSpline(y, x, q, kx=3, ky=3, s=0)
    paired = lambda xx, yy: spline.ev(yy, xx)-spline.ev(yy, -xx)
    physical = [polar_strain(paired, x[-1], y[-1], order=n) for n in (48, 96)]
    spectra = {str(n): fourier_laplace(paired, x[x >= 0], y, cell_order=n)
               for n in (4, 6)}
    desired = float(rates[1]-rates[0])
    interior = physical[-1]['total']
    return dict(name=path.parent.name, source_file=str(path),
                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                bounds=[float(x[0]), float(x[-1]), float(y[0]), float(y[-1])],
                rates=dict(lambda_=float(rates[0]), mu=float(rates[1]),
                           v=float(rates[2]), b=float(rates[3])),
                target_wall_strain_mu_minus_lambda=desired,
                retained_source_strain=interior,
                required_exterior_strain=desired-interior,
                exterior_fraction_of_required_strain=(desired-interior)/desired,
                polar_quadrature_order48=physical[0]['total'],
                polar_quadrature_order96=interior,
                polar_order_change=abs(physical[1]['total']-physical[0]['total']),
                radial_shell_contributions=physical[-1]['shells'],
                horizontal_fourier_laplace_by_cell_order=spectra,
                limitations=[
                    'Uses cubic interpolation of the saved finite-difference source qx.',
                    'Source is set to zero outside the saved rectangle for this diagnostic.',
                    'Required exterior strain is a necessary matching datum, not a computed tail.',
                    'Fourier cutoff dependence and interpolation error are separate uncertainties.',
                    'Polar order change measures quadrature sensitivity, not discretization error.',
                    'No global steady state or identification of a root is verified.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--states', type=Path, default=Path('runs/steady_20260928_tau20'))
    parser.add_argument('--out', type=Path, default=Path('runs/fourier_20260928'))
    parser.add_argument('--manufactured-only', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    result = dict(convention='fhat(k,ell)=int exp(-i*(k*x+ell*y))*f(x,y) dxdy',
                  status='identity validation and finite-source diagnostic only',
                  manufactured=manufactured_checks(), profiles=[])
    print(json.dumps(result['manufactured'], indent=2), flush=True)
    if not args.manufactured_only:
        for name in ('L4_129x65', 'L4_257x129', 'L6_193x97', 'L8_257x129_continued'):
            entry = profile_probe(args.states/name/'state.npz')
            result['profiles'].append(entry)
            print(json.dumps({k: entry[k] for k in (
                'name', 'target_wall_strain_mu_minus_lambda', 'retained_source_strain',
                'required_exterior_strain', 'polar_order_change')}), flush=True)
    output = args.out/'wall_strain_probe.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(str(output.resolve()), flush=True)


if __name__ == '__main__':
    main()
