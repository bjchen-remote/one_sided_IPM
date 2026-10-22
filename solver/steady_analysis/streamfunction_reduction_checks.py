"""Small exact checks for the streamfunction skeleton reduction.

No numerical trajectory, fitting, projection, or global existence claim.
The polynomial G is a local symbolic test chart; the derivation is general.
"""
import json
from pathlib import Path
import sympy as s


def check_zero(expr):
    ans = s.factor(s.cancel(expr))
    assert ans == 0, ans
    return "exact zero"


def run():
    x, y, xi, lam, mu, v = s.symbols("x y xi lambda mu v", real=True)
    n, W, k, B, c, chi = [s.Function(z)(y) for z in ("n", "W", "k", "B", "c", "chi")]
    G = xi + k*xi**2/2
    g = 1+xi**2
    J = xi+xi**3/3+k*(xi**2/2+xi**4/4)
    K = s.integrate(s.expand(J*s.diff(G, xi)), xi)
    X = n+W*G
    L, T = s.diff(X, xi), s.diff(X, y)
    H = c+B*W*J
    Dx = lambda f: s.diff(f, xi)/L
    Dy = lambda f: s.diff(f, y)-T*s.diff(f, xi)/L
    S = Dy(H)
    hd = s.diff(H, y)
    psi = y*H
    V1, V2 = lam*X+v-Dy(psi), lam*y+Dx(psi)
    a = lam+B*g
    U = lam*(X-y*T)+v-H-y*hd
    checks = {}
    checks["J_xi"] = check_zero(s.diff(J, xi)-g*s.diff(G, xi))
    checks["H_x"] = check_zero(Dx(H)-B*g)
    checks["H_y"] = check_zero(S-(hd-B*g*T))
    checks["V2"] = check_zero(V2-y*a)
    checks["moving_advection_cancellation"] = check_zero(V1-T*V2-U)
    q = -Dx(Dx(psi))-Dy(Dy(psi))
    qmix = -2*S-y*(s.diff(S, y)-T*s.diff(S, xi)/L+B*s.diff(g, xi)/L)
    checks["poisson_mixed"] = check_zero(q-qmix)
    # Tangent-space expression uses fixed xi partials in J_k and G_k.
    b = [-B*g, B*(J-g*G), B*W*(s.diff(J,k)-g*s.diff(G,k)), W*J, s.Integer(1)]
    theta = [n,W,k,B,c]
    checks["S_tangent_range"] = check_zero(S-sum(bb*s.diff(th,y) for bb,th in zip(b,theta)))
    checks["wall_normal_jet_xi"] = check_zero(s.diff(S,xi)-(W*s.diff(B,y)*g*s.diff(G,xi)-B*T*s.diff(g,xi)))
    wall_U = lam*X+v-H
    wall_q = -2*S
    wall_gradient_transport = wall_U*Dx(wall_q)+(lam-mu-B*g)*wall_q
    mother_ode = wall_U*(s.diff(B,y)*g-B*T*s.diff(g,xi)/L)+(lam-mu-B*g)*S
    checks["wall_mother_ode"] = check_zero(wall_gradient_transport+2*mother_ode)
    Phi = chi+y*(c*W*G+B*W**2*K)
    checks["potential_primitive"] = check_zero(Dx(Phi)-psi)
    checks["potential_poisson"] = check_zero(Dx(-Dx(Dx(Phi))-Dy(Dy(Phi)))-q)
    expected = [
        -y*H,
        y*B*W*(2*K-J*G),
        y*B*W**2*(s.diff(K,k)-J*s.diff(G,k)),
        y*W**2*K,
        y*W*G,
    ]
    delta_phi = [s.diff(Phi,th)-Dx(Phi)*s.diff(X,th) for th in theta]
    for th,actual,want in zip(theta,delta_phi,expected):
        checks[f"fixed_x_potential_variation_{str(th.func)}"] = check_zero(actual-want)
    # The coefficient of theta'' in Delta Phi is Phi_theta at fixed x.
    lap_phi = Dx(Dx(Phi))+Dy(Dy(Phi))
    for th,want in zip(theta,expected):
        checks[f"highest_second_height_{str(th.func)}"] = check_zero(s.diff(lap_phi,s.diff(th,y,2))-want)
    checks["highest_second_height_chi"] = check_zero(s.diff(lap_phi,s.diff(chi,y,2))-1)
    # Since V has only first height derivatives, the theta''' coefficient
    # in +V.grad(Delta Phi) is +V2 times the above coefficient.
    f = s.Function("Phi")(x,y)
    lap = s.diff(f,x,2)+s.diff(f,y,2)
    rho = -lap
    vx,vy = lam*x+v-s.diff(f,x,y), lam*y+s.diff(f,x,2)
    original = vx*s.diff(rho,x)+vy*s.diff(rho,y)-mu*rho-v
    potential = vx*s.diff(lap,x)+vy*s.diff(lap,y)-mu*lap+v
    checks["scalar_potential_transport_sign"] = check_zero(original+potential)
    # Divergence/gradient identities behind the streamfunction-only test.
    Vx,Vy,Vxx,Vyx,q0,qx,qy,r,rx = s.symbols("Vx Vy Vxx Vyx q qx qy r rx")
    differentiated_transport = Vx*qx+Vy*qy+(Vxx-mu)*q0+Vyx*r
    Kx_minus_muq = Vxx*q0+Vx*qx+Vyx*r+Vy*rx-mu*q0
    checks["density_compatibility_identity"] = check_zero(Kx_minus_muq-differentiated_transport-Vy*(rx-qy))
    return dict(checks=checks, count=len(checks),
                scope="Local exact algebra only; no projected closure, analytic existence, or global Biot-Savart claim.")


if __name__ == "__main__":
    out = Path(__file__).resolve().parents[2]/"runs/streamfunction_reduction_20260928/independent_checks.json"
    result = run()
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
