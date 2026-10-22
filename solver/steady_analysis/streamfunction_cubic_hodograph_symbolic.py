"""Exact wall elimination for a cubic inverse regular streamfunction.

No fitted modes, floating point, or numerical integration.
"""
from pathlib import Path
import json
import sympy as s

OUT=Path(__file__).resolve().parents[2]/'runs/streamfunction_symbolic_20260929'
OUT.mkdir(parents=True,exist_ok=True)
h,a,b,c,k,l,m,n,lam,mu,v=s.symbols('h a b c k l m n lambda mu v',real=True)
X=a*h+b*h*h+c*h**3
p=s.diff(X,h);w=k+l*h+m*h*h+n*h**3
rho=2*k*h+l*h*h+s.Rational(2,3)*m*h**3+n*h**4/2
q=2*w/p
E=s.expand(2*(lam*X+v-h)*w-p*(mu*rho+v))
co=s.Poly(E,h)
checks={}
def check(name,expr):
    out=s.factor(expr)
    assert out==0,(name,out)
    checks[name]='exact zero'
check('rho_h_equals_2w',s.diff(rho,h)-2*w)
check('highest_h6',co.nth(6)-c*n*(4*lam-3*mu)/2)
check('next_h5_after_n0',co.nth(5).subs(n,0)-2*c*m*(lam-mu))
q0=q.subs(h,0)
qx=s.diff(q,h)/p
qxx=s.diff(qx,h)/p
peak={k:a/2,l:b,m:0,n:0}
check('peak_value',(q0-1).subs(peak))
check('peak_first',qx.subs(peak).subs(h,0))
check('peak_curvature',qxx.subs(peak).subs(h,0)+6*c/a**3)
normalized={**peak,c:a**3/3}
check('normalized_curvature',qxx.subs(normalized).subs(h,0)+2)
ep=s.Poly(E.subs(normalized),h)
check('h1_slope',ep.nth(1)-a*(a*(lam-mu)-1))
base={**normalized,lam:mu+1/a,b:v*a**3}
eb=s.Poly(s.expand(E.subs(normalized).subs({lam:mu+1/a,b:v*a**3})),h)
check('h2_drift',s.factor(ep.nth(2).subs(lam,mu+1/a))-(b-a**3*v))
check('h3_remaining',eb.nth(3)-a**3*(6*a**2*v**2-2*a*mu+1)/3)
check('h4_remaining',eb.nth(4)+a**5*v*(a*mu-2)/3)
branch={**normalized,b:v*a**3,lam:3/a,mu:2/a}
epbranch=s.Poly(s.expand(E.subs(normalized).subs({b:v*a**3,lam:3/a,mu:2/a})),h)
check('nonzero_drift_only_obstruction',epbranch.as_expr()-a**3*h**3*(2*a*a*v*v-1))
qbranch=s.factor(q.subs(normalized).subs({b:v*a**3,lam:3/a,mu:2/a}))
check('wall_q_formula',qbranch-(1+2*v*a*a*h)/(1+2*v*a*a*h+a*a*h*h))
# At h=-1/(v*a^2), numerator=-1. Under 2a^2 v^2=1,
# denominator=1, so q=-1 exactly; no limiting argument is needed.
hbad=-1/(v*a*a)
check('explicit_negative_point',(qbranch.subs(h,hbad)+1).subs(v*v,1/(2*a*a)))
check('quadratic_discriminant',(2*v*a*a)**2-4*a*a-(-2*a*a)-2*a*a*(2*a*a*v*v-1))
vzero={**normalized,b:0,v:0,lam:3*mu,a:1/(2*mu)}
# Sequential substitution is intentional for dependencies c=a^3/3, k=a/2.
zero_expr=E.subs(normalized).subs({b:0,v:0,lam:3*mu}).subs(a,1/(2*mu))
check('zero_drift_wall_only',zero_expr)
payload=dict(
    ansatz='x=Z(y)+a(y)h+b(y)h^2+c(y)h^3, h=psi/y',
    assumptions=['lambda>mu>0','v!=0','a(0)>0,c(0)>0','global monotone map of h in R to x in R',
                 'rho(0,0)=0,q(0,0)=1,q_x(0,0)=0,q_xx(0,0)=-2','smooth wall traces'],
    result='No globally positive wall q in this cubic inverse-streamfunction class.',
    coefficientPolynomials={str(i):str(s.factor(co.nth(i))) for i in range(7)},
    elimination=['n=0','m=0','k=a/2','l=b','c=a^3/3','lambda-mu=1/a','b=v*a^3',
                 'mu=2/a, lambda=3/a','2*a^2*v^2=1'],
    wallQ=str(qbranch),negativePoint='h=-1/(v*a^2) gives q=-1',checks=checks,
    count=len(checks),zeroDriftScope='Only a wall-compatible branch is exhibited; no full two-dimensional closure is asserted.')
(OUT/'cubic_hodograph_certificate.json').write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps(dict(result=payload['result'],count=len(checks),checks=checks),indent=2))
