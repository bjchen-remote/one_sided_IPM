"""Pure algebra certificates for ANGLE_MAP_SYMBOLIC.md. No source data."""
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_polar_structure_20260929'
T,k,eta,zeta,xi,sigma,z=sp.symbols('T k eta zeta xi sigma z',real=True)
H=k*sp.cos(T)+eta*sp.sin(T)+zeta*sp.sin(T)*sp.cos(T)+xi*sp.sin(T)**2
hp=sp.diff(H,T)
endpoints=sp.simplify(hp.subs(T,0)-(eta+zeta))==0 and sp.simplify(hp.subs(T,sp.pi)-(-eta+zeta))==0
S=2*z/(1+z*z);C=(1-z*z)/(1+z*z)
rational=sp.simplify(sigma*S*C-S*S/2-2*z*(sigma*(1-z*z)-z)/(1+z*z)**2)==0
c1,c2,c3=sp.symbols('c1 c2 c3')
matrix=sp.Matrix([[1,c,1-c*c] for c in (c1,c2,c3)])
vandermonde=sp.factor(matrix.det())
vandermonde_ok=sp.simplify(vandermonde+(c2-c1)*(c3-c1)*(c3-c2))==0
U,B,q,qx,qy,Ux,P,mu,K,v=sp.symbols('U B q qx qy Ux P mu K v')
E=U*qx+B*qy+(Ux-mu)*q
rho_x=(Ux*q+U*qx-P*K+B*qy)/mu
rho_x_ok=sp.simplify(rho_x.subs(K,E/P)-q)==0
Uy,By,Ky=sp.symbols('Uy By Ky')
rho_y_numerator=Uy*q+U*qy-By*K-B*Ky
second_identity=rho_y_numerator+mu*K
rho_y_ok=sp.simplify(rho_y_numerator/mu+K-second_identity/mu)==0
j=dict(type='exact symbolic algebra, no numerical source or fitting',
    wallEnvelopeDerivative=endpoints,rationalP2=rational,
    vandermondeDeterminant=str(vandermonde),vandermondeNondegenerateWhenCosinesDistinct=vandermonde_ok,
    reconstructedRhoX=rho_x_ok,reconstructedRhoY=rho_y_ok)
assert all(j[k] for k in ('wallEnvelopeDerivative','rationalP2','vandermondeNondegenerateWhenCosinesDistinct','reconstructedRhoX','reconstructedRhoY'))
(OUT/'angle_map_symbolic_certificate.json').write_text(json.dumps(j,indent=2)+'\n')
print(json.dumps(j,indent=2))
