"""Exact structural tests for X_s(s,y)^2=P(s)+y*R(s).
No finite-difference fitting or Taylor-error correction is performed.
"""
from pathlib import Path
import json
import sympy as sp

s,y,mu,lam=sp.symbols('s y mu lam',real=True,nonzero=True)
J=sp.Function('J')(s);K=sp.Function('K')(s)

def wall_source(J,K):
    return (-mu*(1+K*K)*(J*sp.diff(J,s)+s*(J*sp.diff(J,s,2)-sp.diff(J,s)**2))/(2*J**3)
        -(lam-2*mu)*K*sp.diff(K,s)/J
        +mu*s*K*sp.diff(K,s)*sp.diff(J,s)/J**2
        +mu*s*K*sp.diff(K,s,2)/J
        -sp.Rational(3,2)*mu*s*sp.diff(K,s)**2/J+sp.diff(K,s)/J**2)


def source_check():
    Ap=(lam-mu)*J-mu*s*sp.diff(J,s)
    App=sp.diff(Ap,s)
    B=-mu*s*sp.diff(K,s)/2
    # y^1 coefficient of J^3 Delta(Phi), with X=F+K*y.
    Nlap1=(1+K*K)*(J*App-sp.diff(J,s)*Ap)-2*K*J*sp.diff(K,s)*Ap
    Nlap1-=4*K*J**2*sp.diff(B,s)
    Nlap1+=6*B*J**2*sp.diff(K,s)
    expected=Nlap1/(2*J**3)+sp.diff(K,s)/J**2
    assert sp.simplify(expected-wall_source(J,K))==0
    return 'exact symbolic equality with cleared Poisson coefficient'


def symmetric_obstruction():
    j0,j2,k0,k2=sp.symbols('j0 j2 k0 k2',nonzero=True,real=True)
    jj=j0+j2*s*s/2;kk=k0+k2*s*s/2
    wall=mu*s*sp.diff(kk,s)*jj+2*(lam-mu)*kk*jj-2*mu*s*kk*sp.diff(jj,s)-1
    wall2=sp.factor(sp.diff(wall,s,2).subs(s,0))
    assert sp.factor(wall2-2*(lam*j0*k2+(lam-3*mu)*k0*j2))==0
    srcprime=sp.factor(sp.diff(wall_source(jj,kk),s).subs(s,0))
    srcprime=sp.factor(srcprime.subs(k2,-(lam-3*mu)*k0*j2/(lam*j0)))
    srcprime=sp.factor(srcprime.subs(j0,1/(2*(lam-mu)*k0)))
    target=-mu*j2*(2*(lam-mu)*k0)**2*(1+(lam/mu-1-3*mu/lam)*k0*k0)
    assert sp.factor(srcprime-target)==0
    return {'factor':'1+(lambda/mu-1-3*mu/lambda)*K(0)^2',
        'excludedLambdaOverMuLowerBound':float((1+sp.sqrt(13))/2),
        'hypotheses':'J_prime(0)=K_prime(0)=0, J_doubleprime(0)>0, positive J and K; exact m=2 law',
        'sourcePrime':str(srcprime)}


def exact_riccati_check():
    j=sp.Function('j')(s,y);k=sp.Function('k')(s,y);D=sp.Function('D')(s,y)
    R=sp.Function('R')(s)
    V1=mu*s*j+k*D
    # Physical x derivative is d_s/J; physical y derivative d_y-K/J*d_s.
    curl_times_J=sp.diff(D,s)-j*sp.diff(V1,y)+k*sp.diff(V1,s)
    curl_times_J=curl_times_J.subs(sp.diff(j,y),sp.diff(k,s)).subs(sp.diff(k,s),R/(2*j))
    S=j*D;aa=sp.diff(D,s)
    bb=-j*sp.diff(D,y)+R*D/(2*j)+mu*(j+s*sp.diff(j,s))
    cc=sp.diff(D,s)-mu*s*R/2+1
    assert sp.simplify(curl_times_J+1-(-S*sp.diff(k,y)+aa*k*k+bb*k+cc))==0
    return 'exact symbolic zero for the full physical curl/Poisson equation'


def main():
    out={'wallSource':source_check(),'symmetricObstruction':symmetric_obstruction(),
         'fullRiccatiReduction':exact_riccati_check(),
         'farFieldPowerConstraint':'m*(1-alpha)=1, assuming a common fixed-density chart and C1 homogeneous far field at a simple angular zero',
         'scope':'Necessary structural constraints and exact reduced system; not a full classification of asymmetric solutions.'}
    dest=Path(__file__).resolve().parents[2]/'runs/structure_free_alpha_20260928/free_jacobian_checks.json'
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
