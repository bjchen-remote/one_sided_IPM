"""Audit mapped geometry against velocity and Poisson on held-out radii.

Radial parameters are interpolated observations, not solutions of radial ODEs.
The common angular mother is frozen; interior q and velocities are not fitted.
"""
import os,json
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ipm_polar_mpl')
import numpy as np
from scipy.interpolate import CubicSpline,BSpline,RectBivariateSpline
from streamfunction_sources import load_native

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_polar_structure_20260929'

def model(archive,count,R,theta):
    train=np.arange(len(R))%2==0
    t=np.log(R)
    keys=['A','s','kappa','eta','zeta']
    values=np.stack([archive[f'tau20_{count}function_{key}'] for key in keys],axis=1)
    values[:,:2]=np.log(values[:,:2])
    cs=CubicSpline(t[train],values[train],axis=0)
    z,z1,z2=[cs(t,nu=n)[:,:,None] for n in range(3)]
    logA,k,kap,eta,zeta=[z[:,j] for j in range(5)]
    a,w,kap1,eta1,zeta1=[z1[:,j] for j in range(5)]
    a1,w1,kap2,eta2,zeta2=[z2[:,j] for j in range(5)]
    th=2*np.arctan2(np.exp(k)*np.sin(theta/2),np.cos(theta/2))
    S,C=np.sin(th),np.cos(th)
    H=kap*C+eta*S+zeta*S*C
    J=-kap*S+eta*C+zeta*(C*C-S*S)
    J1=-kap*C-eta*S-4*zeta*S*C
    Ht=kap1*C+eta1*S+zeta1*S*C
    Htt=kap2*C+eta2*S+zeta2*S*C
    Jt=-kap1*S+eta1*C+zeta1*(C*C-S*S)
    h=np.cosh(k)+np.sinh(k)*C
    h1=-np.sinh(k)*S
    P=a+Ht+w*S*J
    DP=a1+Htt+w1*S*J+2*w*S*Jt+w*w*S*C*J+w*w*S*S*J1
    L2=h*h+w*w*S*S
    L1=w1*S+w*w*S*C+2*P*w*S+h*h1+2*h*h*J
    L0=DP+P*P+h*h1*J+h*h*(J1+J*J)
    mother=BSpline(archive['referenceMotherSplineKnots'],
        archive['referenceMotherSplineCoefficients'],int(archive['referenceMotherSplineDegree']))
    F,F1,F2=[mother(th,nu=i) for i in range(3)]
    pref=np.exp(logA+H)
    rr=R[:,None]
    psi=pref*F
    ur=-pref*h*(F1+J*F)/rr
    ut=pref*(P*F+w*S*F1)/rr
    q=-pref*(L2*F2+L1*F1+L0*F)/(rr*rr)
    return dict(psi=psi,ur=ur,ut=ut,q=q)

def main():
    d=np.load(OUT/'angle_map_wall_envelope.npz')
    R=d['nativeRadii'];theta=np.linspace(0,np.pi,801)
    RR,TT=np.meshgrid(R,theta,indexing='ij')
    p=load_native(20);x=RR*np.cos(TT)/p['ell'];y=RR*np.sin(TT)/p['ell']
    fs={key:RectBivariateSpline(p['y'],p['x'],p[key]) for key in ('psi','u1','u2','q')}
    u1=p['P']*p['ell']*fs['u1'].ev(y,x)
    u2=p['P']*p['ell']*fs['u2'].ev(y,x)
    truth=dict(psi=p['P']*p['ell']**2*fs['psi'].ev(y,x),
        ur=np.cos(TT)*u1+np.sin(TT)*u2,ut=-np.sin(TT)*u1+np.cos(TT)*u2,
        q=p['P']*fs['q'].ev(y,x))
    baseline=-p['P']*(fs['psi'].ev(y,x,dx=2)+fs['psi'].ev(y,x,dy=2))
    allmodels={count:model(d,count,R,theta) for count in (3,5)}
    held=np.arange(len(R))%2==1
    records=[]
    for count,z in allmodels.items():
        rows=[]
        for i,r in enumerate(R):
            if not held[i]:continue
            rel=lambda a,b:float(np.linalg.norm(a-b)/np.linalg.norm(b))
            rows.append(dict(nativeR=float(r),psiRelativeL2=rel(z['psi'][i],truth['psi'][i]),
                velocityRelativeL2=rel(np.r_[z['ur'][i],z['ut'][i]],np.r_[truth['ur'][i],truth['ut'][i]]),
                qRelativeL2=rel(z['q'][i],truth['q'][i]),qSplineBaselineRelativeL2=rel(baseline[i],truth['q'][i]),
                predictedQMin=float(z['q'][i].min()),sourceQMin=float(truth['q'][i].min()),
                negativePredictedQAnglesDegrees=(theta[z['q'][i]<0]*180/np.pi).tolist(),
                wallQPrediction=z['q'][i,[0,-1]].tolist(),wallQSource=truth['q'][i,[0,-1]].tolist()))
        windows={}
        for name,lo,hi in [('middle',4,256),('outer',256,4096)]:
            rr=[z for z in rows if lo<=z['nativeR']<=hi]
            windows[name]={key:[min(r[key] for r in rr),max(r[key] for r in rr)] for key in
                ('psiRelativeL2','velocityRelativeL2','qRelativeL2','qSplineBaselineRelativeL2','predictedQMin','sourceQMin')}
        records.append(dict(functions=count,rows=rows,windows=windows))
    result=dict(protocol='Alternate native radii calibrate log A, log s, kappa, eta, zeta. Cubic interpolation in log R gives first/second radial derivatives. The other radii are held out. Frozen tau20 R32 mother; no interior velocity or q fitting and no radial ODE solution.',
        source=p['source'],tau=20,exactOperator='q=-A exp(H)/R^2*(L2 Fsecond+L1 Fprime+L0 F)',runs=records)
    (OUT/'mapped_operator_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(OUT/'mapped_operator_audit.npz',nativeRadii=R,theta=theta,heldout=held,
        **{f'source_{k}':v for k,v in truth.items()},**{f'model{n}_{k}':v for n,z in allmodels.items() for k,v in z.items()})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for row,target in enumerate((64.,512.)):
        indices=np.flatnonzero(held);i=indices[np.argmin(abs(np.log(R[indices]/target)))]
        for col,key in enumerate(('psi','ur','q')):
            scale=np.max(abs(truth[key][i]));ax=axs[row,col]
            ax.plot(theta*180/np.pi,truth[key][i]/scale,'k',label='saved field')
            for count,color in [(3,'#36658c'),(5,'#bd663c')]:
                ax.plot(theta*180/np.pi,allmodels[count][key][i]/scale,color=color,ls='--',label=f'{count} radial functions')
            ax.axhline(0,color='.7',lw=.6);ax.set(title=f'R={R[i]:.3g}, {key}',xlabel='Angle (degrees)',ylabel='Divided by source peak magnitude')
            ax.grid(alpha=.15)
    axs[0,0].legend(fontsize=8)
    fig.suptitle('Held-out radii | same mapped streamfunction generates velocity and Poisson source',fontsize=13)
    fig.savefig(OUT/'mapped_operator_audit.png',dpi=180);plt.close(fig)
    print(json.dumps([dict(functions=z['functions'],windows=z['windows']) for z in records],indent=2))

if __name__=='__main__':main()
