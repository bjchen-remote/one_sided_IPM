"""Exact symbolic wall-closure test for the minimal alpha=2/5 family.

This is a test of a particular wall trace, not an exclusion of its exponent.
The first normal equation has an explicit rational-plus-logarithmic solution.
No saved field is modified, and no evolution or full steady solve is run.
"""
from pathlib import Path
import json
import sympy as sp
import mpmath as mp

OUT = Path(__file__).resolve().parents[2] / 'runs/skeleton_closure_20260928'


def symbolic_gate():
    z, k, c = sp.symbols('z k c', real=True)
    a, b = (2 + 5*k)/8, (5*k - 2)/8
    den = (1+k)*z**10 + 1-k
    num = sp.expand(z**2*(a*z**4+c*z**2+b)**2*(a*z**4-b)**4)
    residue_test = sp.expand(z**2*sp.diff(num,z,2)-27*z*sp.diff(num,z)+171*num)
    remainder = sp.Poly(sp.rem(residue_test, den, z), z)
    powers = [0,2,4,6,8]
    polynomials = [sp.fraction(sp.factor(remainder.coeff_monomial(z**i)))[0]
                   for i in powers]
    # These factors cannot vanish in the admissible interval |k|<2/5.
    polynomials[0] = sp.cancel(polynomials[0]/(3*(k-1)*(5*k-2)*(5*k+2)))
    polynomials[4] = sp.cancel(polynomials[4]/(3*(5*k-2)*(5*k+2)))
    basis = list(sp.groebner(polynomials, c, k))
    assert basis == [1]
    return dict(
        exponent='2/5', gamma=3, admissible='-2/5 < k < 2/5; c arbitrary real',
        normalizedR=str(a*z**2+c+b/z**2), Xz='D/(2*z**6)',
        D=str(den), N=sp.sstr(num, order='lex'),
        integrand='128*N/D**3',
        residueCriterion='D divides z**2*N.diff(z,2)-27*z*N.diff(z)+171*N',
        remainderCoefficients={str(i):str(sp.factor(remainder.coeff_monomial(z**i))) for i in powers},
        nonzeroFactorsRemoved=['3*(k-1)*(5*k-2)*(5*k+2)', '3*(5*k-2)*(5*k+2)'],
        reducedPolynomialConstraints=[str(t) for t in polynomials],
        groebnerVariables=['c','k'], groebnerBasis=[str(t) for t in basis],
        conclusion='At least one logarithmic residue is nonzero for every admissible minimal wall trace. Its regular first normal coefficient is not rational in z. This does not exclude alpha=2/5 or more general wall traces.'), (z,k,c,num,den,residue_test)


def explicit_solution_check(symbols):
    mp.mp.dps=60
    z,k,c,num,den,residue_test=symbols
    params={k:sp.Rational(1,5), c:0}
    n=sp.lambdify(z,num.subs(params),'mpmath')
    n1=sp.lambdify(z,sp.diff(num,z).subs(params),'mpmath')
    n2=sp.lambdify(z,sp.diff(num,z,2).subs(params),'mpmath')
    d=sp.lambdify(z,den.subs(params),'mpmath')
    d1=sp.lambdify(z,sp.diff(den,z).subs(params),'mpmath')
    d2=sp.lambdify(z,sp.diff(den,z,2).subs(params),'mpmath')
    d3=sp.lambdify(z,sp.diff(den,z,3).subs(params),'mpmath')
    kk=mp.mpf(1)/5
    aa,bb=(2+5*kk)/8,(5*kk-2)/8
    rr=lambda w:aa*w**2+bb/w**2
    rz=lambda w:2*aa*w-2*bb/w**3
    xz=lambda w:d(w)/(2*w**6)
    qq=lambda w:rz(w)/xz(w)
    zs=(-bb/aa)**mp.mpf('.25')
    radius=((1-kk)/(1+kk))**mp.mpf('.1')
    poles=[radius*mp.exp(mp.j*mp.pi*(2*j+1)/10) for j in range(10)]
    coeff=[]
    for r in poles:
        a3=128*n(r)/d1(r)**3
        a2=128*(n1(r)/d1(r)**3-3*n(r)*d2(r)/(2*d1(r)**4))
        a1=64*(n2(r)/d1(r)**3-3*n1(r)*d2(r)/d1(r)**4
               +n(r)*(3*d2(r)**2-d1(r)*d3(r))/d1(r)**5)
        coeff.append((r,a1,a2,a3))
    def integral(w):
        return sum(-a3/2*((w-r)**-2-(zs-r)**-2)
                   -a2*((w-r)**-1-(zs-r)**-1)
                   +a1*mp.log((w-r)/(zs-r)) for r,a1,a2,a3 in coeff)
    def khat(w):
        return integral(w)/(qq(w)**2*rr(w)**3)
    rows=[]
    for text_z in ['0.45','0.6','0.9','1','1.3','1.7']:
        w=mp.mpf(text_z)
        val=khat(w)
        direct=mp.quad(lambda t:128*n(t)/d(t)**3,[zs,w])/(qq(w)**2*rr(w)**3)
        residual=rr(w)*mp.diff(khat,w)+(3*rz(w)+2*rr(w)*mp.diff(qq,w)/qq(w))*val-rz(w)*qq(w)
        rows.append(dict(z=float(w),muTimesPhysicalXy=float(mp.re(val)),
                         quadratureDifference=float(abs(val-direct)),
                         normalEquationResidual=float(abs(residual)),
                         imaginaryCancellation=float(abs(mp.im(val)))))
    return dict(k='1/5',constantTermC=0,
                note='Representative admissible wall trace only. Closed-form first normal solution, not full IPM solution.',
                stagnationZ=float(zs),precisionDigits=mp.mp.dps,
                maxLogResidue=float(max(abs(a1) for _,a1,_,_ in coeff)),
                checks=rows,
                maxQuadratureDifference=max(r['quadratureDifference'] for r in rows),
                maxNormalEquationResidual=max(r['normalEquationResidual'] for r in rows))


if __name__=='__main__':
    gate,symbols=symbolic_gate()
    result=dict(symbolicGate=gate,explicitFirstNormalSolution=explicit_solution_check(symbols),
                sympyVersion=sp.__version__)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'symbolic_closure_gate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'groebnerBasis':gate['groebnerBasis'],
                      'explicitSolution':result['explicitFirstNormalSolution']},indent=2))
