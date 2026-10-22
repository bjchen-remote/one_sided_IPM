"""Finite Laurent transport closure for a root geometry, before Poisson."""
import json
from pathlib import Path
import sympy as s

t,z,e,lam,delta=s.symbols('t z epsilon lambda delta',real=True)
k,a,b,g,h,P,R,A,H,N=[s.Function(n)(t) for n in ('k','a','b','g','h','P','R','A','H','N')]
x=g+a*z**3+e*z
y=h+b*z**3
rho=t*z-e/z+delta
psi=P*z**4+A*z**3+N*z**2+R*z+H
dt=lambda f:s.diff(f,t)
dz=lambda f:s.diff(f,z)
det=dt(x)*dz(y)-dz(x)*dt(y)
T=dt(psi)*dz(rho)-dz(psi)*dt(rho)
T+=lam*((x*dz(y)-y*dz(x))*dt(rho)+(y*dt(x)-x*dt(y))*dz(rho))
T-=lam*rho*det/3
subs={dt(g):k*dt(h),dt(a):1+k*dt(b),a:t/3+k*b}
T=s.expand(T.subs(subs,simultaneous=True))
q_certificate=s.expand((dt(rho)*dz(y)-dz(rho)*dt(y))-det/z**2)
assert s.factor(q_certificate.subs(subs,simultaneous=True))==0
terms={str(n[0]-2):str(s.factor(c)) for n,c in s.Poly(s.expand(T*z**2),z).terms()}
out=Path(__file__).resolve().parents[2]/'runs/normal_flow_closure_20260928'
out.mkdir(exist_ok=True)
(out/'root_transport.json').write_text(json.dumps(terms,indent=2)+'\n')
print(json.dumps(terms,indent=2))

# Leading angular Poisson equation forced by the degree-four stream function.
# Transport for the finite set of all nonnegative root powers forces b''=0.
B0,B1,kv=s.symbols('B0 B1 k_value',real=True)
bb=B0+B1*t
aa=t/3+kv*bb
pp=e*lam*(B0+s.Rational(8,9)*B1*t)
ds=lambda f:s.diff(f,t)+s.Rational(2,3)/bb*s.diff(f,kv)
xx=aa*z**3
yy=bb*z**3
ps=pp*z**4
dd=s.expand(ds(xx)*s.diff(yy,z)-s.diff(xx,z)*ds(yy))
dx=lambda f:s.cancel((s.diff(yy,z)*ds(f)-ds(yy)*s.diff(f,z))/dd)
dy=lambda f:s.cancel((-s.diff(xx,z)*ds(f)+ds(xx)*s.diff(f,z))/dd)
poisson=s.factor(-dx(dx(ps))-dy(dy(ps))-z**-2)
num,den=s.fraction(poisson)
angular=s.Poly(num,kv)
compat=s.factor(s.resultant(num,s.cancel(3*bb*ds(num)),kv))
assert s.factor(s.Poly(compat,t).LC()-1728*2187*B1**12*e**2*lam**2)==0
angular_result={'residual':str(poisson),'numerator_coefficients_in_k':{
    str(n[0]):str(s.factor(co)) for n,co in angular.terms()},
    'b':str(bb),'K_prime':'2/(3*b)',
    'P4':str(pp),'differential_resultant':str(compat),
    'scope':'Necessary leading Poisson equation for finite nonnegative-root stream function.'}
(out/'root_angular.json').write_text(json.dumps(angular_result,indent=2)+'\n')
print(json.dumps(angular_result,indent=2))

# Independent physical-coordinate substitution of the remaining affine-in-x branch.
X,Y,bb0,kk0=s.symbols('X Y b k0',positive=True)
ss=bb0**s.Rational(2,3)*(X-kk0*Y)*Y**(-s.Rational(2,3))
ss+=s.Rational(9,2)*bb0**s.Rational(4,3)/lam*Y**(-s.Rational(1,3))
pss=-s.Rational(9,4)*bb0**s.Rational(2,3)*Y**s.Rational(4,3)
ux=-s.diff(pss,Y);uy=s.diff(pss,X)
poisson_exact=s.simplify(-s.diff(pss,X,2)-s.diff(pss,Y,2)-s.diff(ss,X))
transport_exact=s.simplify((ux+lam*X)*s.diff(ss,X)+(uy+lam*Y)*s.diff(ss,Y)-lam*ss/3)
assert poisson_exact==0 and transport_exact==0
(out/'explicit_rejected_branch.json').write_text(json.dumps({
    's':str(ss),'psi':str(pss),'poisson_residual':str(poisson_exact),
    'transport_residual':str(transport_exact),
    'q':str(s.diff(ss,X)),
    'rejection':'q is independent of X and diverges at Y=0; no localized smooth wall peak.'},indent=2)+'\n')
