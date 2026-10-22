"""Exact jet certificate for unrestricted rank-one psi=A(r)F(theta).

F is not chosen in advance: smooth origin regularity forces its quadratic
angular jet. The polynomial used below is that forced Cartesian jet, with
free higher radial coefficients and unrestricted density offset jets.
"""
import json
from pathlib import Path
import sympy as s


def run():
    x,y=s.symbols('x y',real=True)
    sigma,lam,mu,v=s.symbols('sigma lambda mu v',real=True)
    u,w,t=s.symbols('u w t',real=True)
    C1,C2,C3=s.symbols('C1 C2 C3',real=True)
    radius2=x*x+y*y
    radial=1+u*radius2+w*radius2**2+t*radius2**3
    psi=(sigma*x*y-y*y/2)*radial
    q=s.expand(-s.diff(psi,x,2)-s.diff(psi,y,2))
    rho=s.integrate(q,x)+C1*y+C2*y*y+C3*y**3
    vx=lam*x+v-s.diff(psi,y)
    vy=lam*y+s.diff(psi,x)
    residual=s.expand(vx*s.diff(rho,x)+vy*s.diff(rho,y)-mu*rho-v)
    wall=s.expand(residual.subs(y,0))
    normal=s.expand(s.diff(residual,y).subs(y,0))
    checks={}

    def exact(name,expr):
        value=s.factor(expr)
        assert value==0,(name,value)
        checks[name]='exact polynomial identity'

    f=1+u*x*x+w*x**4+t*x**6
    exact('wall_q',q.subs(y,0)-f)
    exact('wall_horizontal_velocity',-s.diff(psi,y).subs(y,0)+sigma*x*f)
    exact('wall_density_offset_retained',s.diff(rho,y).subs(y,0)-(C1-6*sigma*u*x*x-8*sigma*w*x**4-10*sigma*t*x**6))
    exact('wall_even_transport',(wall+wall.subs(x,-x))/2-v*(f-1))
    exact('curvature_normalization',s.diff(q,x,2).subs({x:0,y:0})-2*u)
    exact('wall_rate_relation',wall.coeff(x,1)-(lam-mu-sigma))
    wall_mu=s.expand(wall.subs(mu,lam-sigma))
    coeff3=s.factor(wall_mu.coeff(x,3))
    coeff5=s.factor(wall_mu.coeff(x,5))
    exact('wall_cubic_coefficient',coeff3-u*(2*lam-5*sigma)/3)
    exact('wall_quintic_coefficient',coeff5-(w*(4*lam-9*sigma)/5-sigma*u*u))
    exact('normal_constant',normal.coeff(x,0).subs(mu,lam-sigma)-(1+2*sigma*C1))
    forced={lam:s.Rational(5,2)*sigma,mu:s.Rational(3,2)*sigma,
            C1:-1/(2*sigma),w:5*u*u,v:0}
    normal4=s.factor(normal.coeff(x,4).subs(forced))
    exact('normal_fourth_obstruction',normal4+u*u*(652*sigma*sigma+3)/2)
    assert not normal4.has(t,C2,C3)
    checks['higher_jets_cannot_change_obstruction']='exact absence of t,C2,C3'

    # General commutator identity, with no modal assumption for density.
    p=s.Function('psi')(x,y); r=s.Function('rho')(x,y)
    Vx,Vy=lam*x+v-s.diff(p,y),lam*y+s.diff(p,x)
    T=Vx*s.diff(r,x)+Vy*s.diff(r,y)-mu*r-v
    qx=s.diff(r,x)
    gradient=Vx*s.diff(qx,x)+Vy*s.diff(qx,y)+(s.diff(Vx,x)-mu)*qx+s.diff(p,x,2)*s.diff(r,y)
    exact('commutator_density_elimination',s.diff(T,x)-gradient)

    return dict(count=len(checks),checks=checks,
                certificates=dict(wall_even=str(v*(f-1)),wall_x3=str(coeff3),
                                  wall_x5=str(coeff5),normal_constant='1+2*sigma*C1',
                                  normal_x4_after_forced_relations=str(normal4)),
                conclusion='No sufficiently smooth global rank-one psi=A(r)F(theta) satisfies the centered wall equations and q(0,0)=1,q_x(0,0)=0,q_xx(0,0)=-2. Both v!=0 and v=0 are excluded; rho was not assumed separable.',
                scope='Rank one must extend to the regular origin and both wall rays. Does not exclude rank two or rank-one descriptions restricted to annuli/angular bands.')


if __name__=='__main__':
    result=run()
    root=Path(__file__).resolve().parents[2]
    path=root/'runs/streamfunction_separation_20260929/rank_one_certificate.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
