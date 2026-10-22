"""Held-out fourth/fifth ridge orders test finite symbolic resummations.

Normal coefficients are computed by the existing regular integral recurrence;
spectral degree refinement is mandatory. Candidate ridge expressions are
fixed by lower orders, never fitted to the orders used for checking them.
"""
import json
from pathlib import Path
import signal
import time
import numpy as np
import sympy as sp
from sympy.polys.fields import field
from math import factorial
from ridge_generated_jets import WallLift,ridge_series
from normal_lift_numeric import NormalLift

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_generated_20260928'


class WallAdapter(WallLift):
    def jacobian(self,r):
        j=self.J(r,2)
        return j[0],j[1],2*j[2]


def power_series(p,n):
    a=[1.]
    for k in range(1,n+1):a.append(a[-1]*(p-k+1)/k)
    a[0]=0.
    return np.array(a)


def candidates(mt):
    n=len(mt)-1;d1,d2,d3=mt[1:4]
    # Coalescing powers: t^3 and t^3 log(t).
    logarithm=np.zeros(n+1)
    for k in range(1,n+1):logarithm[k]=(-1)**(k+1)/k
    base=power_series(3.,n);unit=base.copy();unit[0]=1
    logbasis=np.convolve(unit,logarithm)[:n+1]
    beta=(d2-d1)/1.5;alpha=(d1-beta)/3
    logpred=alpha*base+beta*logbasis
    # Two distinct powers with p=3. q and both amplitudes use only d1,d2,d3.
    q=(d1-3*d3)/(d1-d2)
    second=power_series(q,n)
    co=np.linalg.solve(np.array([[3.,q],[3.,q*(q-1)/2]]),np.array([d1,d2]))
    powerpred=co[0]*base+co[1]*second
    assert np.max(abs(powerpred[:4]-mt[:4]))<1e-10
    result=[dict(name='t3_plus_t3_log_t',fittedOrders=[1,2],parameters=dict(alpha=alpha,beta=beta),
                 coefficients=logpred.tolist(),heldOutAbsoluteErrors=(logpred[3:]-mt[3:]).tolist()),
            dict(name='t3_plus_tq',fittedOrders=[1,2,3],parameters=dict(q=q,alpha=co[0],beta=co[1]),
                 coefficients=powerpred.tolist(),heldOutAbsoluteErrors=(powerpred[4:]-mt[4:]).tolist())]
    for factor in [0.,3.]:
        invfactor=power_series(-factor,n);invfactor[0]=1
        reduced=np.convolve(mt,invfactor)[:n+1]
        d1,d2,d3,d4=reduced[1:5]
        b1,b2=np.linalg.solve(np.array([[d2,d1],[d3,d2]]),[-d3,-d4])
        numer=np.array([0.,d1,d2+b1*d1]);denom=np.array([1.,b1,b2])
        pred=np.zeros(n+1)
        for i in range(1,n+1):pred[i]=(numer[i] if i<=2 else 0)-b1*pred[i-1]-(b2*pred[i-2] if i>=2 else 0)
        outer=power_series(factor,n);outer[0]=1
        pred=np.convolve(pred,outer)[:n+1]
        roots=np.roots(denom[::-1]);positive=[float(z.real) for z in roots if abs(z.imag)<1e-10 and z.real>0]
        result.append(dict(name=f't_power_{int(factor)}_times_pade22',fittedOrders=[1,2,3,4],
            numerator=numer.tolist(),denominator=denom.tolist(),coefficients=pred.tolist(),
            heldOutAbsoluteErrors=(pred[5:]-mt[5:]).tolist(),positiveZPoles=positive,
            farPower=int(factor)))
    return result


def exact_fourth_control():
    r,L=sp.symbols('r L');F,rr,ll=field((r,L),sp.QQ)
    saved=json.loads((ROOT/'runs/ridge_density_closure_20260928/fourth_elementary_gate.json').read_text())
    a=F.from_expr(sp.sympify(saved['coefficient'],locals={'r':r,'L':L}))
    direct=WallLift(0.,3.,1.,1.,20).fourth(0.,2)
    errors=[];reference=[]
    for n in range(3):
        exact=a.as_expr().subs(r,0).subs(L,sp.pi/4).evalf(50)
        error=abs(float(exact)-factorial(n)*direct[n])
        assert error<1e-10
        errors.append(float(error));reference.append(str(exact))
        a=a.diff(rr)+a.diff(ll)/(1+rr*rr)
    return dict(scope='Independent exact symmetric fourth normal coefficient and first two density derivatives at Q=0.',
                exactReference=reference,absoluteErrors=errors)


def main():
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('45 second local bound')))
    signal.alarm(45);start=time.monotonic()
    fits=json.loads((ROOT/'runs/ridge_density_closure_20260928/asymmetric_wall_data.json').read_text())
    fit=next(r for r in fits if r['name']=='L4_257x129')
    lam,mu,v,_=fit['rates'];b=fit['fit']['b'];results=[]
    for resonance in [False,True]:
        wall=WallAdapter(b,3*mu if resonance else lam,mu,v)
        direct=wall.third_ridge()
        rows=[]
        for degree in [48,64,80,96]:
            lift=NormalLift(wall,count=5,degree=degree,quadrature=96)
            j=wall.J(0.,6);a0=np.r_[0.,j/np.arange(1,8)]
            coefficients=[a0]+[np.array([a.deriv(k)(wall.d)/factorial(k) for k in range(8-n)])
                                for n,a in enumerate(lift.coeff,1)]
            ridge=ridge_series(coefficients,5)
            error=float(np.max(abs(np.array(ridge['m_tSeries'][:4])-direct['m_tSeries'])))
            row=dict(degree=degree,lowOrderIndependentError=error,**ridge)
            rows.append(row)
            print('ridge fifth gate',resonance,degree,ridge['m_tSeries'],'independent',error,flush=True)
        assert max(row['lowOrderIndependentError'] for row in rows)<2e-6
        assert np.ptp([row['m_tSeries'][5] for row in rows])<1e-3
        checks=candidates(np.array(rows[-1]['m_tSeries']))
        independent=[]
        if resonance:
            for nodes in [12,16,20]:
                direct4=WallLift(b,3*mu,mu,v,nodes).fourth_ridge()
                independent.append(dict(nodes=nodes,**direct4))
                print('independent fourth',nodes,direct4['m_tSeries'],flush=True)
            assert np.max(abs(np.array(independent[-1]['m_tSeries'])-independent[-2]['m_tSeries']))<1e-7
            assert abs(independent[-1]['m_tSeries'][4]-rows[-1]['m_tSeries'][4])<2e-6
        results.append(dict(enforcedResonance=resonance,refinement=rows,candidates=checks,
                            independentFourthOrder=independent,
                            scope='Numerical finite-order necessary-condition checks; not an exact symbolic no-go theorem.'))
    result=dict(name=fit['name'],b=b,rates=fit['rates'],results=results,
                exactFourthControl=exact_fourth_control(),seconds=time.monotonic()-start)
    (OUT/'ridge_resummation_gate.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0)
    print('seconds',result['seconds'],flush=True)


if __name__=='__main__':main()
