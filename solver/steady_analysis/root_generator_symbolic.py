"""Exact certificate excluding the quadratic-Q, gamma=2/3 generator.

All arithmetic is rational.  No sampling, fitting, numerical root search,
or large Groebner computation is performed.  Work in z^3=Q; differentiate
using z_x=Q_x/(3z^2), z_y=Q_y/(3z^2).
"""
import json
from pathlib import Path
import sympy as s


def run():
    x,y,z=s.symbols('x y z')
    a,b,c,d,e,A=s.symbols('a b c d e A')
    lam,mu,v,h,hp=s.symbols('lambda mu v h h_prime')
    Q=1+d*x+e*y+a*x*x+2*b*x*y+c*y*y
    Qx,Qy=s.diff(Q,x),s.diff(Q,y)
    Dx=lambda f:s.diff(f,x)+Qx*s.diff(f,z)/(3*z*z)
    Dy=lambda f:s.diff(f,y)+Qy*s.diff(f,z)/(3*z*z)
    relation=z**3-Q
    checks={}

    def exact(name,expr):
        num,den=s.fraction(s.cancel(expr))
        rem=s.rem(num,relation,z)
        rem=s.expand(rem)
        assert rem==0,(name,s.factor(rem))
        checks[name]='exact zero modulo z^3-Q'

    t=2*A/s.Integer(9)
    kappa=2*A*d/s.Integer(3)
    Phi=A*y*z*z-kappa*x*y
    P=s.expand(3*Q*(y*(s.diff(Q,x,2)+s.diff(Q,y,2))+2*Qy)-y*(Qx**2+Qy**2))
    U=s.expand(3*Q*(Qx+y*s.diff(Q,x,y))-y*Qx*Qy)
    W=s.expand(y*(3*Q*s.diff(Q,x,2)-Qx**2))
    Nx=s.expand(3*Q*s.diff(P,x)-4*P*Qx)
    Ny=s.expand(3*Q*s.diff(P,y)-4*P*Qy)
    exact('laplacian_generator',Dx(Dx(Phi))+Dy(Dy(Phi))-t*P/z**4)
    exact('mixed_derivative',Dx(Dy(Phi))-(t*U/z**4-kappa))
    exact('second_x_derivative',Dx(Dx(Phi))-t*W/z**4)
    exact('laplacian_x_derivative',Dx(t*P/z**4)-t*Nx/(3*z**7))
    exact('laplacian_y_derivative',Dy(t*P/z**4)-t*Ny/(3*z**7))

    E0=lam*y*hp-mu*h+v
    E1=s.expand(-U*Nx+W*Ny)
    E2=(lam*x+v+kappa)*Nx+lam*y*Ny-3*mu*P*Q+3*W*Q*hp
    vx=lam*x+v+kappa-t*U/z**4
    vy=lam*y+t*W/z**4
    residual=vx*t*Nx/(3*z**7)+vy*(t*Ny/(3*z**7)+hp)-mu*(t*P/z**4+h)+v
    expected=3*Q**4*E0+t*t*z*E1+t*Q*z*z*E2
    exact('three_root_classes',3*z**12*residual-expected)

    f=1+d*x+a*x*x
    wall_gate=(6*b-d*e)+(4*b*d-2*a*e)*x+2*a*b*x*x
    wall_expected=-18*f*f*(d+2*a*x)*wall_gate
    assert s.expand(E1.subs(y,0)-wall_expected)==0
    checks['wall_factorization']='exact polynomial identity'

    branch=s.Poly(E1.subs({b:0,e:0}),x,y)
    coef_x6y=s.factor(branch.coeff_monomial(x**6*y))
    coef_y7=s.factor(branch.coeff_monomial(y**7))
    assert s.expand(coef_x6y-36*a**4*(a+9*c))==0
    assert s.expand(coef_y7-12*a*c**3*(3*a+7*c))==0
    checks['branch_b_e_zero_x6_y']='exact coefficient identity'
    checks['branch_b_e_zero_y7']='exact coefficient identity'
    # On a != 0 the first equation fixes c=-a/9.  The second is nonzero.
    substituted=s.factor(coef_y7.subs(c,-a/9))
    assert s.expand(substituted+s.Rational(80,2187)*a**5)==0
    resultant=s.factor(s.resultant(a+9*c,c**3*(3*a+7*c),c))
    assert resultant==-20*a**4
    checks['a_nonzero_branch_contradiction']='exact substitution and resultant'

    linear_branch=s.Poly(s.expand(E1.subs({a:0,b:0,e:0})),y)
    linear_y=s.factor(linear_branch.coeff_monomial(y))
    assert s.expand(linear_y+9*d**4*(1+d*x))==0
    checks['a_b_e_zero_y']='exact coefficient identity'

    off_wall_branch=s.Poly(s.expand(E1.subs({a:0,d:0})),y)
    off_wall_y=s.factor(off_wall_branch.coeff_monomial(y))
    assert off_wall_y==-432*b*b
    checks['a_d_zero_y']='exact coefficient identity'
    trivial=s.expand(E1.subs({a:0,b:0,d:0}))
    assert trivial==0
    assert s.expand(Nx.subs({a:0,b:0,d:0}))==0
    checks['remaining_branch_has_zero_q']='exact polynomial identity'

    return dict(
        candidate='Phi=A*y*Q^(2/3)-(2*A*d/3)*x*y+chi(y), Q=1+d*x+e*y+a*x^2+2*b*x*y+c*y^2',
        definitions={key:str(s.factor(expr)) for key,expr in [('P',P),('U',U),('W',W),('Nx',Nx),('Ny',Ny)]},
        checks=checks,count=len(checks),
        certificates=dict(wall_factor=str(s.factor(wall_expected)),wall_gate=str(wall_gate),
                          b_e_zero_x6_y=str(coef_x6y),b_e_zero_y7=str(coef_y7),
                          c_minus_a_over_9_y7=str(substituted),resultant=str(resultant),
                          a_b_e_zero_y=str(linear_y),a_d_zero_y=str(off_wall_y)),
        conclusion='For A!=0 and x-dependent Q, the independent root class requires E1=0; its complete coefficient case split forces a=b=d=0, a contradiction. Thus no nontrivial q solution in this entire generator family, even with arbitrary chi(y).',
        scope='Q>0 on any nonempty open physical patch; chi sufficiently differentiable and y-only; gamma exactly 2/3; lambda,mu,v arbitrary constants. No claim about other generators or exponents.')


if __name__=='__main__':
    result=run()
    root=Path(__file__).resolve().parents[2]
    out=root/'runs/streamfunction_symbolic_20260929/root_generator_certificate.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'count':result['count'],'checks':result['checks'],
                      'certificates':result['certificates'],'conclusion':result['conclusion']},indent=2))
