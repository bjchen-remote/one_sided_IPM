"""Wall traces of the data-learned polar streamfunction modes.

Use each learned angular mode's exact linear combination of source circles
to recover wall derivatives, instead of extrapolating midpoint splines.
This audits a necessary closure gate, not a full radial solution.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from streamfunction_sources import load_native

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'runs/streamfunction_fresh_20260928'
OUT=ROOT/'runs/streamfunction_polar_closure_20260928'

def score(pred,truth):
    e=np.linalg.norm(pred-truth)/np.linalg.norm(truth)
    return dict(relativeL2=float(e),maxRelative=float(np.max(abs((pred-truth)/truth))))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    p=load_native(20);lam,mu,v,off=p['rates']
    source={k:RectBivariateSpline(p['y'],p['x'],p[k]) for k in ('psi','u1','q','Q')}
    with np.load(OLD/'polar_geometry.npz') as z:
        R=z['nativeRadii'];r=R/p['ell'];train=z['trainingRadiusMask']
        psi=z['wall_saddle_tau20_psi']/(p['P']*p['ell']**2)
        amplitude=np.sqrt(np.mean(psi**2,axis=1));phi=psi/amplitude[:,None]
    X=r[:,None]*np.array([1.,-1.]);Y=np.zeros_like(X)
    qw=source['q'].ev(Y,X);uw=source['u1'].ev(Y,X);rho=source['Q'].ev(Y,X)
    # theta derivative at 0 is r*psi_y=-r*u1; at pi it is -r*psi_y=r*u1.
    first=np.stack([-r*uw[:,0],r*uw[:,1]],axis=1)
    second=-r[:,None]**2*qw
    second_spline=r[:,None]**2*source['psi'].ev(Y,X,dx=2)
    rows=[]
    for band,lo,hi in [('middle',4.,256.),('outer',256.,4096.)]:
        idx=np.flatnonzero((R>=lo)&(R<=hi));tr=idx[train[idx]];te=idx[~train[idx]]
        u,s,vh=np.linalg.svd(phi[tr],full_matrices=False)
        for rank in (1,2,3,4):
            basis=vh[:rank];comb=u[:,:rank].T/s[:rank,None]
            d1=comb@(first[tr]/amplitude[tr,None])
            d2=comb@(second[tr]/amplitude[tr,None])
            d2s=comb@(second_spline[tr]/amplitude[tr,None])
            # Best possible held-out coefficients: observations, not ODE predictions.
            aa=psi[te]@basis.T
            first_pred=aa@d1;second_pred=aa@d2
            up=first_pred/(r[te,None]*np.array([-1.,1.]))
            qp=-second_pred/r[te,None]**2
            # Compare source and reconstructed wall transport, using observed rho.
            xp=X[te]
            src=(lam*xp+v+uw[te])*qw[te]-mu*rho[te]-v
            pred=(lam*xp+v+up)*qp-mu*rho[te]-v
            denom=np.sqrt(np.sum(((lam*xp+v+uw[te])*qw[te])**2+(mu*rho[te]+v)**2))
            rows.append(dict(band=band,rank=rank,trainingNativeRadii=R[tr].tolist(),heldoutNativeRadii=R[te].tolist(),
                modeWallFirst=d1.tolist(),modeWallSecond=d2.tolist(),
                boundaryPsiExactlyZeroByLinearCombination=True,
                velocityWalls=score(up,uw[te]),qWalls=score(qp,qw[te]),
                leftQ=score(qp[:,1],qw[te,1]),rightQ=score(qp[:,0],qw[te,0]),
                minPredictedQ=float(qp.min()),minSavedQ=float(qw[te].min()),
                secondWallDerivativeSourceSensitivity=float(np.linalg.norm(d2-d2s)/np.linalg.norm(d2)),
                sourceWallTransportRelative=float(np.linalg.norm(src)/denom),
                predictedWallTransportWithObservedDensity=float(np.linalg.norm(pred)/denom),
                sourceTransportNormalization='L2 norm of the pair (Vx*q, mu*rho+v)',
                qObserved=qw[te].tolist(),qPredicted=qp.tolist()))
    output=dict(protocol=dict(source='native tau20, centered inner coordinates',
        rate=dict(lam=float(lam),mu=float(mu),v=float(v)),
        modes='Same training-radius normalized-psi SVD as fresh polar observation; m=1..4 diagnostics',
        coefficientTest='Heldout coefficients reprojected from full observed psi: conditional shape gate, not prediction',
        wallTrace='Modes are linear combinations of training profiles; extend each source circle to theta=0,pi with exact zero value and derivatives from saved u1 and q',
        limitation='Nonsteady finite-time data; original wall residual is separately reported'),rows=rows)
    (OUT/'wall_closure.json').write_text(json.dumps(output,indent=2))
    lines=['# 极坐标流函数模式的壁面闭合检查','','两模态在内部半圆上的速度拟合好，并不自动保证壁面密度梯度。这里不从角向中点样条外推，而利用每个 SVD 模式是训练半圆的线性组合，直接继承源场的壁面导数。', '',
        '每个源半圆的 psi 在两端为零，psi_theta(0)=-r*u1(r,0)、psi_theta(pi)=r*u1(-r,0)，psi_thetatheta=-r^2*q。因此有限数据模态可以有非零的二阶端点迹；有限正弦特征模态却不能。', '',
        '| 区间 | 模式数 | 壁面速度误差 | 壁面 q 误差 | 源壁面输运余量 | 重构后壁面输运余量 |', '|---|---:|---:|---:|---:|---:|']
    for d in rows:
        lines.append(f"| {d['band']} | {d['rank']} | {100*d['velocityWalls']['relativeL2']:.3f}% | {100*d['qWalls']['relativeL2']:.3f}% | {100*d['sourceWallTransportRelative']:.3f}% | {100*d['predictedWallTransportWithObservedDensity']:.3f}% |")
    lines.extend(['','这里留出半径的模式系数和密度都仍使用观测值；只审核模式的必要壁面能力，不是径向 ODE 的预测。源输运余量也非零，因为 tau20 不是严格稳态。全部角模式都在训练半径确定，没有在测试壁面重新拟合。', '',
        '对于有限模式的精确解，还存在两个壁面代数条件。记 f0=Phi_prime(0), h0=Phi_second(0)，fpi、hpi 类似，则', '',
        '−[lambda*r+v−a·f0/r]*(a·h0)/r^2 = mu*b·G(0)+v,', '',
        '[lambda*r−v−a·fpi/r]*(a·hpi)/r^2 = mu*b·G(pi)+v.', '',
        '它们来自壁面 Poisson 与输运联立；只有内部角向投影不能替代这些条件。Phi 的函数值在两端为零，并不要求 Phi_second 在两端为零。', '',
        '复现：solver/steady_analysis/fresh_polar_wall_closure.py；完整结果 wall_closure.json。'])
    (OUT/'WALL_CLOSURE.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps([{k:d[k] for k in ('band','rank','velocityWalls','qWalls','sourceWallTransportRelative','predictedWallTransportWithObservedDensity')} for d in rows],indent=2))

if __name__=='__main__': main()
