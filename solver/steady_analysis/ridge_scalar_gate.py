"""Symbolic checks for the curvature gate and scalar skeleton compatibility."""
import json
from pathlib import Path
import sympy as sp

u,a,f,A,Ap,b,k,mp,h,T0=sp.symbols('U a f A Ap beta k mp h T0')
r0,r1,r0p,r1p,p0,p1,p0p,p1p=sp.symbols('r0 r1 r0p r1p p0 p1 p0p p1p')
R=r1*u+r0
Rp=(r1p+a*r1)*u+r0p+f*r1
Pp=(p1p+a*p1)*u+p0p+f*p1
px=(mp*Pp-Rp)/(1+mp*mp)
rx=(Pp+mp*Rp)/(1+mp*mp)
assert sp.simplify(mp*rx-px-Rp)==0
assert sp.simplify(rx+mp*px-Pp)==0
physical=(T0-b*u)*h+2*R*Ap+A*(px-k*rx)
E1=-b*h+2*Ap*r1+A*(b*(p1p+a*p1)-(1+k*mp)*(r1p+a*r1))/(1+mp*mp)
E0=T0*h+2*Ap*r0+A*(b*(p0p+f*p1)-(1+k*mp)*(r0p+f*r1))/(1+mp*mp)
assert sp.factor((physical-E1*u-E0).subs(b,mp-k))==0
e0,e1,e0p,e1p=sp.symbols('E0 E1 E0p E1p',nonzero=True)
quotient_derivative=(-e0p*e1+e0*e1p)/e1**2
cleared=sp.expand(e1**2*(quotient_derivative-a*(-e0/e1)-f))
assert sp.expand(cleared-(-e0p*e1+e0*e1p+a*e0*e1-f*e1**2))==0
out=Path(__file__).resolve().parents[2]/'runs/ridge_parameter_20260928'
result=dict(hessianReconstructionCheck='0',curvatureAffineInVelocityCheck='0',
    eliminatedScalarEquation='-E0_prime*E1+E0*E1_prime+a*E0*E1-f*E1^2=0',
    amplitudeParameterEquation='-dot(E0)*E1+E0*dot(E1)+Ydot*a*E0*E1-Ydot*f*E1^2=0',
    velocityRecovery='U=-E0/E1, only on nondegenerate intervals E1 != 0',
    orderBound='At most fourth order in Y(t) if A,m,C,curvatureWidth are prescribed functions of amplitude parameter t.',
    rationality='Laurent-polynomial A,m,C,curvatureWidth imply rational differential coefficients in t,Ydot,Yddot,Ydddot,Yddddot.',
    scope='Exact necessary ridge and curvature conditions. No claim of complete 2D symbolic closure, global positivity, or Biot-Savart selection. E1=0 requires the original coupled equations.')
(out/'scalar_gate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
