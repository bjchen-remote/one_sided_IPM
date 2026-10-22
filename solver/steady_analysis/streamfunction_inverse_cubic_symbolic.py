"""Exact transport coefficients for psi=y F(x-Z(y)), F'=1/(a+bF+cF^2)."""
import json
from pathlib import Path
import sympy as s

OUT=Path(__file__).resolve().parents[2]/'runs/streamfunction_symbolic_20260929'
OUT.mkdir(parents=True,exist_ok=True)
f=s.symbols('F')
a,b,c,y,lam,mu,D,A,B,Ap,Bp,C,Cp,v=s.symbols('a b c y lambda mu D A B Ap Bp C Cp v',real=True)
Q=a+b*f+c*f*f
z=a*f+b*f*f/2+c*f**3/3
P=1/Q
rho=A*f-B*P+C
q=A*P-B*P*s.diff(P,f)
R=(lam*z+D-f)*q+y*(lam+P)*(Ap*f-Bp*P+Cp)-mu*rho-v
poly=s.Poly(s.cancel(R*Q**3),f)
coeff={str(i):s.factor(poly.nth(i)) for i in range(poly.degree()+1)}
print('Degree:',poly.degree())
for i in reversed(range(poly.degree()+1)):
    print(i,coeff[str(i)])
(OUT/'inverse_cubic_coefficients.json').write_text(json.dumps(dict(
    ansatz='psi=y F(x-Z(y)); z=a F+b F^2/2+c F^3/3; F_prime=1/(a+bF+cF^2)',
    definitions='A=2 Z_prime+y Z_second; B=y(1+Z_prime^2); Ap=dA/dy; Bp=dB/dy; D=lambda(Z-yZ_prime)+v; C=density addition',
    coefficientPolynomials={k:str(vv) for k,vv in coeff.items()}),indent=2))
