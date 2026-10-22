"""Free-exponent wall generator identities and normalized illustrative profiles.

No data are fitted here. The plotted profiles satisfy first-normal wall
compatibility, not the full two-dimensional equation.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from scipy.integrate import cumulative_simpson


def symbolic_checks():
    s,r,A,mu=sp.symbols('s r A mu',positive=True)
    B,H=sp.Function('B')(s),sp.Function('H')(s)
    K=A*H**(-sp.Rational(2,3))
    J=B/(2*mu*A*H**sp.Rational(1,3))
    wall=mu*s*sp.diff(K,s)*J+2*mu*r*K*J-2*mu*s*K*sp.diff(J,s)-1
    expected=(r*B-s*sp.diff(B,s)-H)/H
    assert sp.simplify(wall-expected)==0
    return dict(generator_identity='wall_residual = (r*B-s*B_prime-H)/H',
        alpha='1/(r+1)',canonical_transverse_rate='mu*(2*r+1)',
        integer_resonance='[s^n]H=(r-n)[s^n]B; vanishes automatically when r=n',
        scope='Wall-normal compatibility only; no full 2D or data identification claim.')


def normalized_examples():
    z,r=sp.symbols('z r',real=True)
    b,c,e,delta=2.,1.,-.55,-.35
    w=1+z*z;chi=z/sp.sqrt(w);p=r/2
    B=b*(1+delta*chi)+c*(w**p-p*w**(p-1)
        +e*z*w**((r-1)/2)-e*(r-1)/2*z*w**((r-3)/2))
    H=r*b+b*delta*chi*(r-1/w)-c/2*w**(p-2)*(r*(r-2)+e*(r-1)*(r-3)*chi)
    f=B/H**sp.Rational(1,3)
    fun=sp.lambdify((z,r),f,'numpy')
    f1=sp.lambdify((z,r),sp.diff(f,z),'numpy')
    f2=sp.lambdify((z,r),sp.diff(f,z,2),'numpy')
    hfun=sp.lambdify((z,r),H,'numpy')
    from structure_free_wall import counterterm_positivity_bound
    rows=[];curves=[]
    for rr in [1.5,1.8,2.,2.2,2.6]:
        sufficient_bound=counterterm_positivity_bound(rr,c,e)/(1-abs(delta))
        assert b>sufficient_bound
        # These examples have one ordinary minimum. Find a sign change rather
        # than assuming it occurs at the stagnation density s=0.
        sample=np.linspace(-10,10,2001)
        dd=f1(sample,rr)
        ids=np.flatnonzero((dd[:-1]<0)&(dd[1:]>=0))
        assert len(ids)==1
        j=ids[0];zp=brentq(lambda u:f1(u,rr),sample[j],sample[j+1])
        fp=fun(zp,rr);L=np.sqrt(f2(zp,rr)/(2*fp));d=L*zp
        Q=np.linspace(-2.5,2.5,4001)
        zz=(Q+d)/L
        JJ=fun(zz,rr)/fp
        xx=cumulative_simpson(JJ,x=Q,initial=0);xx-=xx[len(xx)//2]
        muK=fp/(2*hfun(zz,rr)**(2/3))
        rows.append(dict(r=rr,alpha=1/(rr+1),b=b,c=c,eps=e,delta=delta,
            density_scale=L,density_at_peak=d,z_at_peak=zp,
            required_b_lower_bound=sufficient_bound,
            H_global_lower_bound=rr*(b*(1-abs(delta))-counterterm_positivity_bound(rr,c,e)),
            peak_J=1.,peak_J_prime=0.,peak_J_second=2.,
            K_right_over_left=((1+delta)/(1-delta))**(-2/3)))
        curves.append((rr,xx,1/JJ,Q,muK))
    return rows,curves


def main():
    root=Path(__file__).resolve().parents[2]
    out=root/'runs/structure_free_alpha_20260928'
    out.mkdir(exist_ok=True)
    examples,curves=normalized_examples()
    result=dict(symbolic=symbolic_checks(),examples=examples,
        note='Illustrative normalized wall family; no saved numerical profile was used to choose parameters.')
    (out/'free_generator_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    for rr,xx,qq,Q,muK in curves:
        label=f'alpha={1/(rr+1):.3f}'
        ax[0].plot(xx,qq,label=label)
        ax[1].plot(Q,muK,label=label)
    ax[0].set(xlim=(-3,3),ylim=(0,1.05),xlabel='x',ylabel='wall Q_x',title='Same peak normalization, continuous exponent')
    ax[1].set(xlim=(-2,2),xlabel='Q (peak at 0)',ylabel='mu K',title='Exact first-normal compatibility')
    for a in ax:a.grid(alpha=.2)
    ax[0].legend(fontsize=8)
    fig.suptitle('Analytic wall candidates; full two-dimensional closure remains open')
    fig.savefig(out/'free_alpha_wall_examples.png',dpi=160)
    plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
