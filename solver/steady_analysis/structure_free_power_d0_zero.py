"""Exact Poisson obstruction for X=H0+h*s+k*y+d*s*(C(s)+y)^beta.
Only the affine-H, d!=0, C(0)>0 branch is considered.
"""
from pathlib import Path
import json
import sympy as sp

beta,lam,c0,u=sp.symbols('beta lam c0 u',positive=True)
h,k,d,c1=sp.symbols('h k d c1',real=True,nonzero=True)
mu=lam*(1-beta);nu=lam*(1+beta)
q=d*u**beta;q0=d*c0**beta
J=h+q;Jy=d*beta*u**(beta-1);Js=2*c1*Jy
T=h*(u-c0)+d*(u**(beta+1)-c0**(beta+1))/(beta+1)
Ts=2*c1*(q-q0)
D=nu*T/J
Ds=(nu-mu)*Ts/J-nu*T*Js/J**2
Dy=nu*(J*J-T*Jy)/J**2
# At s=0, K=k and K_y=0. The full curl condition times J is zero.
curl=(1+k*k)*Ds-k*J*Dy+k*Jy*D+mu*k*J+1
actual=sp.simplify(curl*J**2/(2*lam*beta))
a=(1+k*k)*c1;g=1/(2*lam*beta)
E=c0*((beta+1)*h+q0)
target=2*a*(q-q0)*(h+q)+((beta+1)*h*q+q*q-d*E*u**(beta-1))*(k*(h+q)-a)-k*(h+q)**3+g*(h+q)**2
assert sp.simplify(sp.expand(actual-target))==0

# Extract the polynomial-in-q sector without any numerical exponent matching.
Q,Q0,aa,gg=sp.symbols('Q Q0 a g',real=True)
polynomial=sp.Poly(sp.expand(2*aa*(Q-Q0)*(h+Q)+((beta+1)*h*Q+Q*Q)*(k*(h+Q)-aa)-k*(h+Q)**3+gg*(h+Q)**2),Q)
assert polynomial.degree()==2
coefficient2=sp.factor(polynomial.nth(2).subs(gg,k*(h+Q0)))
a_selected=-k*(beta*h+Q0)
assert sp.factor(coefficient2.subs(aa,a_selected))==0
coefficient_negative=sp.factor((-d*c0*((beta+1)*h+Q0)*(k*h-aa)).subs(aa,a_selected))
assert sp.factor(coefficient_negative+d*c0*k*((beta+1)*h+Q0)**2)==0
coefficient1=sp.factor(polynomial.nth(1).subs({gg:k*(h+Q0),aa:a_selected}).subs(Q0,-(beta+1)*h))
assert coefficient1==k*h*h

out={'fullPoissonIdentity':'exact symbolic zero after restriction to s=0',
     'coefficient_u_2beta_before_d_squared':str(coefficient2),
     'selected_a':str(a_selected),
     'coefficient_u_beta_minus_1':str(coefficient_negative),
     'remaining_coefficient_u_beta_before_d':str(coefficient1),
     'excludedBranch':'0<beta<1, beta=1-mu/lambda, d!=0, C(0)>0, H affine, finite positive wall J; smooth local Poisson/transport solution',
     'notExcluded':['non-affine H','nonlinear D','d0!=0 branch','singular C(0)=0 branch','general free-beta ansatz'],
     'degenerate_d_zero':'affine-density solution remains possible; it has no wall-gradient peak'}
path=Path(__file__).resolve().parents[2]/'runs/structure_free_alpha_20260928/free_power_d0_zero_checks.json'
path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
