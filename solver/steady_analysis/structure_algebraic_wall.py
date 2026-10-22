"""A globally positive, asymmetric radical wall family at lambda=3*mu.

The family solves the FIRST wall-normal compatibility equation exactly.
It is not asserted to solve the two-dimensional IPM problem. Parameters are
identified from the stagnation offset alone (C=2), and optionally from that
offset plus the peak normal slope. Full wall data are held out for comparison.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.optimize import brentq, least_squares
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.integrate import cumulative_simpson

z, eps, cc = sp.symbols('z eps C', real=True)
H=1+eps*z/sp.sqrt(1+z*z)
f=(1+2*eps*z*sp.sqrt(1+z*z)+cc*z*z)*H**(-sp.Rational(1,3))
fun=sp.lambdify((z,eps,cc),f,'numpy')
der=sp.lambdify((z,eps,cc),sp.diff(f,z),'numpy')
der2=sp.lambdify((z,eps,cc),sp.diff(f,z,2),'numpy')


def normalized(e,c):
    hi=1.
    while der(hi,e,c)<0 and hi<1e4:
        hi*=2
    zp=brentq(lambda u:der(u,e,c),0,hi)
    fp=fun(zp,e,c)
    ell=np.sqrt(der2(zp,e,c)/(2*fp))
    hp=1+e*zp/np.sqrt(1+zp*zp)
    return dict(eps=e,C=c,zpeak=zp,L=ell,d=ell*zp,
                fpeak=fp,muKpeak=fp/(4*hp**(2/3)))


def exact_checks():
    # Dimensionless s=z (L=1); general L follows by substitution z=s/L.
    A,mu=sp.symbols('A mu',positive=True)
    K=A*H**(-sp.Rational(2,3))
    J=f/(4*mu*A)
    residual=mu*z*sp.diff(K,z)*J+4*mu*K*J-2*mu*z*K*sp.diff(J,z)-1
    zero=sp.simplify(sp.expand_power_base(residual,force=True))
    assert zero==0
    assert sp.simplify(sp.diff(H,z,2).subs(z,0))==0
    return dict(wall_compatibility='exact symbolic zero',
        resonance_obstruction='H_second_derivative_at_zero = 0',
        parameter_domain='-1 < eps < 0, C > 2*abs(eps), L,A,mu > 0',
        positivity='H > 0 and N >= 1+2*abs(eps)*z^2-2*abs(eps)*abs(z)*sqrt(1+z^2) > 0',
        status='Exact wall trace/first normal jet family; full PDE unsolved.')


def compare(d,model):
    x,Q=d['x'],d['Q'][0]
    mu=d['rates'][1]
    xq=PchipInterpolator(Q,x)
    q_spline=CubicSpline(x,d['qx'][0])
    # Fixed common density interval; all four saved wall traces cover it.
    qq=np.linspace(-1.1,1.1,2001)
    e,c,L,shift,fp=[model[k] for k in ['eps','C','L','d','fpeak']]
    zz=(qq+shift)/L
    jj=fun(zz,e,c)/fp
    xp=cumulative_simpson(jj,x=qq,initial=0)
    xp-=xp[len(qq)//2]
    true_x=xq(qq)
    true_q=q_spline(true_x)
    predicted_K=fp/(4*mu)*(1+e*zz/np.sqrt(1+zz*zz))**(-2/3)
    # Restrict to a common geometric core as well as the density window.
    core=abs(true_x)<=2
    common=core & (true_x>=xp[0]) & (true_x<=xp[-1])
    same_x_q=PchipInterpolator(xp,1/jj)(true_x[common])
    return dict(wall_inverse_x_linf=float(np.max(abs(xp[core]-true_x[core]))),
        wall_inverse_x_rms=float(np.sqrt(np.mean((xp[core]-true_x[core])**2))),
        wall_q_same_density_linf=float(np.max(abs(1/jj[core]-true_q[core]))),
        wall_q_same_x_linf=float(np.max(abs(same_x_q-true_q[common]))),
        same_x_interval=[float(true_x[common][0]),float(true_x[common][-1])],
        compared_Q=[float(qq[core][0]),float(qq[core][-1])],
        predicted_peak_K=float(model['muKpeak']/mu)), (qq, true_x, xp, true_q, 1/jj, predicted_K)


def main():
    root=Path(__file__).resolve().parents[2]
    out=root/'runs/structure_20260928'
    out.mkdir(exist_ok=True)
    rows=[]; plot_data=[]
    from fixed_point import fd_matrix
    for name in ['L4_129x65','L4_257x129','L6_193x97','L8_257x129_continued']:
        d=np.load(root/'runs/steady_20260928_tau20'/name/'state.npz')
        lam,mu,v,b=d['rates']; offset=v/mu
        qy=(fd_matrix(d['y'],1,7)@d['Q'])[0]
        i=np.argmin(abs(d['x']))
        Kpeak=-qy[i]/d['qx'][0,i]
        e=brentq(lambda a:normalized(a,2)['d']-offset,-.99999,-.0001)
        one=normalized(e,2)
        one_metrics,pp=compare(d,one)
        # Optional two-scalar matching; C>2|eps| guarantees positive J.
        def unpack(u):
            ee=-1/(1+np.exp(-u[0]))
            return normalized(ee,2*abs(ee)+np.exp(u[1]))
        def residual(u):
            m=unpack(u)
            return [(m['d']-offset)/offset,
                    (m['muKpeak']-mu*Kpeak)/(mu*Kpeak)]
        initial=[np.log(-e/(1+e)),np.log(2-2*abs(e))]
        opt=least_squares(residual,initial,bounds=([-8,-12],[12,6]),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=200)
        two=unpack(opt.x)
        two_metrics,_=compare(d,two)
        rows.append(dict(name=name,measured_alpha=float(mu/lam),measured_d=float(offset),
            measured_Kpeak=float(Kpeak),
            C2_from_offset_only=dict(parameters=one,comparison=one_metrics),
            two_scalars=dict(parameters=two,matching_defect=residual(opt.x),
                optimizer_success=bool(opt.success),optimizer_status=int(opt.status),
                optimizer_message=str(opt.message),active_bounds=opt.active_mask.tolist(),
                exact_scalar_match=bool(np.max(abs(np.array(residual(opt.x))))<1e-7),
                comparison=two_metrics),
            held_out='Full wall shape; the C=2 model also holds out the normal slope.',
            caution='Finite-box data are not known to converge; model fixes alpha=1/3.'))
        plot_data.append((name,pp))
    result=dict(exact=exact_checks(),data=rows)
    (out/'algebraic_wall_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for a,(name,p) in zip(ax.flat,plot_data):
        qq,xx,xp,truth,pred,K=p
        a.plot(xx,truth,label='saved box candidate',lw=2)
        a.plot(xp,pred,'--',label='radical wall, C=2',lw=1.7)
        a.set(xlim=(-2,2),ylim=(.2,1.05),title=name,xlabel='x',ylabel='wall Q_x')
        a.grid(alpha=.2)
    ax[0,0].legend(fontsize=8)
    fig.suptitle('One-parameter wall structure: offset selected, wall shape held out')
    fig.savefig(out/'algebraic_wall_comparison.png',dpi=160)
    plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
