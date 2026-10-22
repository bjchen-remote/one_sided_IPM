"""Exact operator audit and strict finite-skeleton falsification gates."""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from construct_profile import Wall, flux_audit
from normal_lift_algebra import curl_coefficients
from normal_lift_numeric import NormalLift

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/normal_lift_20260928'


def exact_gates():
    s = sp.symbols('s', real=True)
    J = sp.Function('J')(s)
    a = sp.Function('a')(s)
    lam, mu = sp.symbols('lambda mu')
    checks = []
    for n in range(1, 8):
        j = [J]+[0]*(n-1)+[sp.diff(a, s)]
        js = [sp.diff(v, s) for v in j]
        jss = [sp.diff(v, s) for v in js]
        k = [0]*(n-1)+[n*a]
        ks = [sp.diff(v, s) for v in k]
        ky = [(i+1)*k[i+1] for i in range(len(k)-1)]
        e = curl_coefficients(j, js, jss, k, ks, ky, s, lam, mu, n-1, sp.expand)[-1]
        e0 = curl_coefficients([J], [sp.diff(J,s)], [sp.diff(J,s,2)],
                              [], [], [], s, lam, mu, n-1, sp.expand)[-1]
        expected = -n*J**3*(mu*s*sp.diff(a,s)+(n*(2*lam-mu)-mu-(n+1)*mu*s*sp.diff(J,s)/J)*a)
        assert sp.simplify(e-e0-expected)==0
        checks.append(n)
    saved = json.loads((OUT/'cubic_symbolic_k.json').read_text())
    k = sp.symbols('k', positive=True)
    jets = [sp.sympify(row['a_n'], locals={'s':s,'k':k}) for row in saved['coefficients']]
    a1,a2,a3,a4,a5 = jets
    det = sp.factor(a2*a2-a1*a3)
    B = sp.factor(a1*a4-a2*a3)
    C = sp.factor(a3*a3-a2*a4)
    # Cleared-denominator gate, valid independently of rational evaluation.
    pade_gate = sp.factor(det*a5+B*a4+C*a3)
    root_gate = sp.factor(24*a2*a2*a5-60*a2*a3*a4+35*a3**3)
    p = sp.Rational(1,3)
    j = 1+s*s
    b = sp.diff(a1,s)/(p*j)
    c = sp.diff(a2,s)/(p*j)-(p-1)*b*b/2
    predicted_j3 = j*(p*(p-1)*b*c+p*(p-1)*(p-2)*b**3/6)
    jgate = sp.factor(sp.diff(a3,s)-predicted_j3)
    result = dict(triangularOperatorChecks=checks,
        exactPade22FifthCoefficientGate=str(pade_gate),
        exactCurvatureQuadraticPowerGate=str(root_gate),
        exactJacobianQuadraticPowerThirdCoefficientGate=str(jgate),
        jacobianQuadraticCoefficientAtSzero=str(sp.factor(c.subs(s,0))),
        scope='Exact cubic-wall benchmark; rejects a universal resummation of this wall, not every possible wall mother function.')
    (OUT/'exact_gates.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def numerical_gates():
    results = []
    for name in ['L4_257x129','tau20','L8_257x129_continued']:
        saved = np.load(ROOT/'runs/construction_20260928'/f'{name}.npz')
        wall = Wall(saved['rates'],saved['parameters'])
        lift = NormalLift(wall,5,80)
        r = np.linspace(-1.25,.5,2001)
        s = r+wall.d
        a = lift.values(r)
        p = (1-wall.mu/wall.lam)/2
        J = wall.jacobian(r)[0]
        j1,j2,j3 = [lift.coeff[i].deriv()(s) for i in range(3)]
        b = j1/(p*J)
        c = j2/(p*J)-(p-1)*b*b/2
        guess3 = J*(p*(p-1)*b*c+p*(p-1)*(p-2)*b**3/6)
        roots = []
        for rr,bb,cc in zip(r,b,c):
            for root in np.roots([cc,bb,1]):
                if abs(root.imag)<1e-10 and root.real>0:
                    roots.append([float(rr),float(root.real)])
        def determinant(rr):
            a1,a2,a3 = lift.values(rr)[:3]
            return a2*a2-a1*a3
        d = determinant(r)
        crossings = []
        for i in np.flatnonzero(d[:-1]*d[1:]<0):
            rr = brentq(determinant,r[i],r[i+1],xtol=1e-14)
            a1,a2,a3,a4 = lift.values(rr)[:4]
            crossings.append(dict(r=float(rr),B=float(a1*a4-a2*a3)))
        refinements = []
        if name=='L4_257x129':
            for nr,ny in [(161,81),(321,161)]:
                rr = np.linspace(-1.25,.5,nr)
                yy = np.linspace(0,.1,ny)
                for kind in ['taylor1','taylor5']:
                    audit = flux_audit(wall,rr,yy,lift.evaluate(rr,yy,kind),lift.integral(rr,yy,kind))
                    audit.update(kind=kind,nr=nr,ny=ny,heightMax=.1)
                    refinements.append(audit)
        results.append(dict(name=name,
            jacobianQuadraticPower=dict(exponent=p,
                thirdCoefficientRelativeMismatch=float(np.linalg.norm(guess3-j3)/np.linalg.norm(j3)),
                negativeQuadraticCoefficientFraction=float(np.mean(c<0)),
                sampledPositiveRoots=len(roots),earliestSampledPositiveRoot=min(x[1] for x in roots) if roots else None),
            padeDeterminantCrossings=crossings,fullPDERefinement=refinements))
    (OUT/'numerical_gates.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


if __name__=='__main__':
    print(json.dumps(exact_gates(),indent=2),flush=True)
    print(json.dumps(numerical_gates(),indent=2),flush=True)
