"""Observe the two zero branches of psi, without a peak/width template.

Exact factorization psi=y*(x-Z(y))*sigma(y)*A(x-Z(y),y).
Learn transverse modes of A on alternating heights, interpolate coefficients
to held-out heights, and audit psi, velocity and Poisson derivatives.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/ipm_streamfunction_mpl')
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from streamfunction_sources import load_native, load_conditional
from streamfunction_geometry import Fields

OUT = Path(__file__).resolve().parents[2] / 'runs/streamfunction_fresh_20260928'

def rel(a,b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))

def zero(f,y):
    xs=np.linspace(max(f.p['x'][0],-4),min(f.p['x'][-1],16),1001)
    vals=f.H(xs,y); ii=np.flatnonzero(vals[:-1]*vals[1:]<=0)
    roots=[brentq(lambda x:float(f.H(np.array(x),y)),xs[i],xs[i+1]) for i in ii]
    return min(roots,key=lambda z:abs(z-1.4*y))

def run(p,ymax,radius):
    f=Fields(p); y=np.linspace(0,ymax,81); s=np.linspace(-radius,radius,241)
    tr=np.arange(0,len(y),2); te=np.arange(1,len(y)-1,2)
    z=np.array([zero(f,float(h)) for h in y])
    assert np.min(z)-radius>=p['x'][0] and np.max(z)+radius<=p['x'][-1]
    sig=np.array([float(f.D_from_psi(np.array(zz),float(h))) for zz,h in zip(z,y)])
    hh=np.array([f.H(zz+s,float(h)) for zz,h in zip(z,y)])
    aa=np.empty_like(hh); mid=len(s)//2
    aa[:,:mid]=hh[:,:mid]/s[:mid]/sig[:,None]
    aa[:,mid+1:]=hh[:,mid+1:]/s[mid+1:]/sig[:,None]; aa[:,mid]=1.
    u,sv,vh=np.linalg.svd(aa[tr],full_matrices=False)
    zs=CubicSpline(y[tr],z[tr]); ss=CubicSpline(y[tr],sig[tr])
    hte=y[te]; zpred=zs(hte); sigpred=ss(hte)
    S,Y=np.meshgrid(s,hte); X=S+zpred[:,None]
    truthpsi=f.psi.ev(Y,X); truthu1=f.u1.ev(Y,X)
    truthu2=f.u2.ev(Y,X); truthq=f.q.ev(Y,X)
    results=[]; reconstructed={}
    for rank in (1,2,3,4):
        cc=(u[:,:rank]*sv[:rank])*sig[tr,None]
        cs=CubicSpline(y[tr],cc); modes=CubicSpline(s,vh[:rank].T)
        m0=modes(s).T; m1=modes(s,1).T; m2=modes(s,2).T
        a=cs(hte)@m0; ay=cs(hte,1)@m0; ayy=cs(hte,2)@m0
        ax=cs(hte)@m1; axy=cs(hte,1)@m1; axx=cs(hte)@m2
        H=S*a; Hs=a+S*ax; Hss=2*ax+S*axx
        Hy=S*ay; Hsy=ay+S*axy; Hyy=S*ayy
        zp=zs(hte,1)[:,None]; zpp=zs(hte,2)[:,None]
        hy=Hy-zp*Hs; hyy=Hyy-2*zp*Hsy+zp*zp*Hss-zpp*Hs
        psi=Y*H; u1=-H-Y*hy; u2=Y*Hs; q=-2*hy-Y*(Hss+hyy)
        interior=(hte>=.1*ymax)&(hte<=.9*ymax)
        sl=np.s_[interior,4:-4]
        results.append(dict(rank=rank,trainFactorRelativeL2=float(np.sqrt(np.sum(sv[rank:]**2)/np.sum(sv**2))),
            heldoutPsi=rel(psi,truthpsi),heldoutU1=rel(u1,truthu1),heldoutU2=rel(u2,truthu2),
            heldoutQInterior=rel(q[sl],truthq[sl]),minimumQInterior=float(q[sl].min()),
            sourceMinimumQInterior=float(truthq[sl].min())))
        reconstructed[rank]=psi
    # Independent source spline derivative floor, with the same physical audit points.
    qsrc=-f.psi.ev(Y,X,dx=2)-f.psi.ev(Y,X,dy=2)
    detail=dict(name=p['name'],kind=p['kind'],yRange=[0,ymax],shiftedXRange=[-radius,radius],
        normalization='A=H/[(x-Z)*sigma], sigma=Hx(Z,y), A(0,y)=1',
        sampling='81 heights, even 41 for training, odd 40 held out; coefficients, Z and sigma interpolated from training heights',
        singularValues=sv.tolist(),results=results,sourceSplineQRelativeL2=rel(qsrc[sl],truthq[sl]),
        minSigma=float(sig.min()),minFactor=float(aa.min()),maxFactor=float(aa.max()),
        zeroPredictionMaximumError=float(np.max(abs(zpred-z[te]))),
        y=y.tolist(),Z=z.tolist(),sigma=sig.tolist())
    return detail,(s,y,aa,sv,vh,z,sig)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    records=[]; plots=[]
    for tau in (10,15,20):
        r,pl=run(load_native(tau),2.,2.); records.append(r)
        if tau==20: plots=pl
    for name in ('L4_257x129','L8_257x129_continued'):
        r,_=run(load_conditional(name),1.,2.); records.append(r)
    r,_=run(load_native(20),4.,4.); records.append(r)
    (OUT/'zero_geometry.json').write_text(json.dumps(dict(records=records),indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    s,y,aa,sv,vh,z,sig=plots
    fig,ax=plt.subplots(1,3,figsize=(13.5,4),layout='constrained')
    for i in (0,10,20,40,80): ax[0].plot(s,aa[i],label=f'y={y[i]:g}')
    ax[0].set(xlabel='s = x - Z(y)',ylabel='H / [s sigma(y)]',title='Direct zero-line factor; no peak alignment');ax[0].legend()
    ax[1].semilogy(np.arange(1,9),sv[:8]/sv[0],'o-');ax[1].set(xlabel='mode',ylabel='singular value / first',title='Data-learned transverse modes')
    for r in records[:3]:
        ax[2].plot(r['y'],r['Z'],label=r['name'].replace('native_',''))
    ax[2].set(xlabel='inner y',ylabel='Z(y)',title='Internal zero streamline');ax[2].legend()
    fig.savefig(OUT/'zero_geometry.png',dpi=160);plt.close(fig)
    print(json.dumps([{k:r[k] for k in ('name','yRange','shiftedXRange','results','sourceSplineQRelativeL2')} for r in records],indent=2))

if __name__=='__main__': main()
