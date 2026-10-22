"""Exact local checks for the new angular-mode radial DAE.

No old peak/halfwidth/mother-profile representation is used.
"""
import json
from pathlib import Path
import sympy as s


def zero(expression):
    value = s.trigsimp(s.simplify(expression))
    assert value == 0, value
    return "exact zero"


def run():
    r=s.symbols('r', positive=True)
    t=s.symbols('theta', real=True)
    lam,mu,v=s.symbols('lambda mu v', nonzero=True, real=True)
    a=[s.Function(f'a{i}')(r) for i in range(2)]
    b=[s.Function(f'b{i}')(r) for i in range(2)]
    F=[s.Function(f'F{i}')(t) for i in range(2)]
    G=[s.Function(f'G{i}')(t) for i in range(2)]
    psi=sum(ai*fi for ai,fi in zip(a,F))
    rho=sum(bi*gi for bi,gi in zip(b,G))
    Dx=lambda f:s.cos(t)*s.diff(f,r)-s.sin(t)*s.diff(f,t)/r
    Dy=lambda f:s.sin(t)*s.diff(f,r)+s.cos(t)*s.diff(f,t)/r
    lap=s.diff(psi,r,2)+s.diff(psi,r)/r+s.diff(psi,t,2)/r**2
    ux,uy=-Dy(psi),Dx(psi)
    vr=lam*r+v*s.cos(t)-s.diff(psi,t)/r
    vt=-v*s.sin(t)+s.diff(psi,r)
    P=lap+Dx(rho)
    T=vr*s.diff(rho,r)+vt*s.diff(rho,t)/r-mu*rho-v
    checks={}
    checks['cartesian_to_polar_laplacian']=zero(Dx(Dx(psi))+Dy(Dy(psi))-lap)
    checks['radial_velocity']=zero((ux+v)*s.cos(t)+uy*s.sin(t)+lam*r-vr)
    checks['angular_velocity']=zero(-(ux+v)*s.sin(t)+uy*s.cos(t)-vt)
    physical=(ux+lam*r*s.cos(t)+v)*Dx(rho)+(uy+lam*r*s.sin(t))*Dy(rho)-mu*rho-v
    checks['transport']=zero(physical-T)
    for i in range(2):
        checks[f'P_a_second_{i}']=zero(s.diff(P,s.diff(a[i],r,2))-F[i])
        checks[f'T_a_second_{i}']=zero(s.diff(T,s.diff(a[i],r,2)))
        checks[f'P_b_first_{i}']=zero(s.diff(P,s.diff(b[i],r))-s.cos(t)*G[i])
        checks[f'T_b_first_{i}']=zero(s.diff(T,s.diff(b[i],r))-vr*G[i])
    # Boundary DAE after using the Poisson equation to eliminate rho_r.
    Vr,L2,R0,R1=s.symbols('Vr L2 R0 R1')
    P0=L2/r**2+R1; T0=Vr*R1-mu*R0-v
    A0=-Vr*L2/r**2-mu*R0-v
    Ppi=L2/r**2-R1; Api=Vr*L2/r**2-mu*R0-v
    checks['wall_theta_zero']=zero(T0-A0-Vr*P0)
    checks['wall_theta_pi']=zero(T0-Api+Vr*Ppi)
    Vx,Vy,Vx0,Vy0,q0,qx,I,Iy,C=s.symbols('Vx Vy Vx0 Vy0 q0 qx I Iy C')
    Cp=(mu*C+v-Vx0*q0)/Vy0
    direct=Vx*qx+Vy*(Iy+Cp)-mu*(I+C)-v
    reduced=Vx*qx+Vy*Iy-mu*I+(Vy/Vy0-1)*(mu*C+v)-(Vy/Vy0)*Vx0*q0
    checks['streamfunction_only_density_offset']=zero(direct-reduced)
    # A finite angular exact solution disproves any blanket finite-mode no-go.
    x,y,q=s.symbols('x y q', real=True)
    strain=lam-mu
    psia=strain*x*y-q*y*y/2
    rhoa=q*x-q*q*y/(2*strain)+v*(q-1)/mu
    va=(lam*x+v-s.diff(psia,y),lam*y+s.diff(psia,x))
    checks['affine_example_poisson']=zero(-s.diff(psia,x,2)-s.diff(psia,y,2)-s.diff(rhoa,x))
    checks['affine_example_transport']=zero(va[0]*s.diff(rhoa,x)+va[1]*s.diff(rhoa,y)-mu*rhoa-v)
    # Same-power Mellin pair, with the density constant removed exactly.
    alpha=s.symbols('alpha', real=True)
    f=s.Function('F')(t);g=s.Function('G')(t)
    psim=r**(alpha+1)*f
    rhom=r**alpha*g-v/mu
    Q=alpha*s.cos(t)*g-s.sin(t)*s.diff(g,t)
    N=(alpha+1)*f*s.diff(g,t)-alpha*s.diff(f,t)*g
    vrm=lam*r+v*s.cos(t)-s.diff(psim,t)/r
    vtm=-v*s.sin(t)+s.diff(psim,r)
    tm=vrm*s.diff(rhom,r)+vtm*s.diff(rhom,t)/r-mu*rhom-v
    expected=r**(2*alpha-1)*N+(lam*alpha-mu)*r**alpha*g+v*r**(alpha-1)*Q
    checks['Mellin_density_derivative']=zero(Dx(rhom)-r**(alpha-1)*Q)
    checks['Mellin_transport_powers']=zero(s.expand_power_base(tm-expected,force=True))
    return {'checks':checks,'count':len(checks),
            'scope':'Exact local algebra only; no numerical modal closure or global half-plane solution is asserted.'}


if __name__=='__main__':
    result=run()
    root=Path(__file__).resolve().parents[2]
    path=root/'runs/streamfunction_polar_closure_20260928/exact_reduction_checks.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
