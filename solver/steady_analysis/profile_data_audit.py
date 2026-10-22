"""Saved-field geometry evidence audit; no PDE solve or symbolic closure search."""
from pathlib import Path
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline, PchipInterpolator, RectBivariateSpline

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/profile_geometry_20260928'
BETA=1.5
XI=np.linspace(-1.5,1.5,401)


def master(x):
    return 32/(x**6-3*x**5+3*x**4+3*x**3+28*x*x+32)


def metrics(a,b):
    e=np.asarray(a)-np.asarray(b)
    return dict(relative_l2_reference=float(np.linalg.norm(e)/np.linalg.norm(b)),
                relative_l2_observed=float(np.linalg.norm(e)/np.linalg.norm(a)),
                absolute_rms=float(np.sqrt(np.mean(e*e))),
                linf=float(np.max(abs(e))))


def sample_indices(n,stride,phase):
    return np.unique(np.r_[0,np.arange(phase,n,stride),n-1]).astype(int)


def measure(x,row,Qfield,h,mode='cubic',stride=1,phase=0):
    ids=sample_indices(len(x),stride,phase)
    xs,qs=x[ids],row[ids]
    cs=CubicSpline(xs,qs,extrapolate=False)
    if mode=='pchip':
        curve=PchipInterpolator(xs,qs,extrapolate=False)
        im=int(np.argmax(qs));m=float(xs[im]);A=float(qs[im])
        qxx=None
    else:
        curve=cs
        roots=cs.derivative().roots(extrapolate=False)
        pts=np.r_[xs[0],roots[(roots>=xs[0])&(roots<=xs[-1])],xs[-1]]
        m=float(pts[np.argmax(cs(pts))]);A=float(cs(m));qxx=float(cs(m,2))
        if mode=='quadratic5':
            i=int(np.searchsorted(xs,m));i=max(2,min(i,len(xs)-3))
            xx=xs[i-2:i+3];yy=qs[i-2:i+3]
            coeff=np.polyfit(xx-m,yy,2)
            if coeff[0]>=0: return dict(available=False,reason='local quadratic not concave')
            dm=-coeff[1]/(2*coeff[0])
            if not xx[0]<=m+dm<=xx[-1]: return dict(available=False,reason='local peak outside stencil')
            m+=float(dm);A=float(np.polyval(coeff,dm));qxx=float(2*coeff[0])
    half=(PchipInterpolator(xs,qs-A/2,extrapolate=False) if mode=='pchip'
          else CubicSpline(xs,qs-A/2,extrapolate=False)).roots(extrapolate=False)
    left=half[half<m];right=half[half>m]
    if not len(left) or not len(right):
        return dict(available=False,reason='half-height not inside saved x domain',m=m,A=A)
    wl=m-float(left[-1]);wr=float(right[0])-m
    W=(wl+wr)*BETA/(2*np.sinh(BETA))
    k=(wr-wl)/(wr+wl)/np.tanh(BETA/2)
    c=float(Qfield(h,m)[0,0]-Qfield(0,0)[0,0])
    wcurv=float(np.sqrt(-2*A/qxx)) if qxx is not None and qxx<0 else None
    ii=int(np.argmin(abs(x-m)))
    return dict(available=True,m=m,A=A,W=W,k=k,left=wl,right=wr,c=c,
                qxx=qxx,wcurv=wcurv,c_over_AW=c/(A*W),
                normalized_xi_curvature=None if qxx is None else qxx*W*W/A,
                local_dx=float(np.median(np.diff(x[max(0,ii-2):min(len(x),ii+3)]))),
                width_over_local_dx=float(W/np.median(np.diff(x[max(0,ii-2):min(len(x),ii+3)]))))


def delta(alt,base):
    if not alt['available'] or not base['available']:return None
    d=dict(peak_shift_over_W=(alt['m']-base['m'])/base['W'],
           A_relative=alt['A']/base['A']-1,W_relative=alt['W']/base['W']-1,
           k_absolute=alt['k']-base['k'],c_over_AW_absolute=alt['c_over_AW']-base['c_over_AW'])
    if alt['qxx'] is not None:
        d['qxx_relative']=alt['qxx']/base['qxx']-1
        d['curvature_width_relative']=alt['wcurv']/base['wcurv']-1
    return d


def audit_native_frames():
    heights=np.array([0,.25,1,4,16,64,256,1024.])
    allrows=[];scales=[];crosscheck={}
    for tau in [10,15,20]:
        d=loadmat(ROOT/f'runs/geometry_loop_20260928/frame_tau{tau}.mat',simplify_cells=True)
        a,P,ell=[float(d[k]) for k in ['a','P','ell']]
        x=a+ell*d['x'];y=ell*d['y'];q=P*d['q'];Q=P*ell*d['Q']
        qf=RectBivariateSpline(y,x,q);Qf=RectBivariateSpline(y,x,Q)
        scales.append(dict(tau=float(d['tau']),a=a,P=P,ell=ell,
                           native_height_1_in_inner=1/ell))
        methods=[('pchip',1,0),('quadratic5',1,0)]
        methods += [('cubic',2,j) for j in range(2)]
        methods += [('cubic',4,j) for j in range(4)]
        for h in heights:
            row=qf(h,x)[0];base=measure(x,row,Qf,h)
            alternatives={}
            for method,stride,phase in methods:
                alt=measure(x,row,Qf,h,method,stride,phase)
                alternatives[f'{method}_stride{stride}_phase{phase}']=dict(geometry=alt,delta=delta(alt,base))
            for phase in range(2):
                ids=sample_indices(len(y),2,phase)
                qfy=RectBivariateSpline(y[ids],x,q[ids]);Qfy=RectBivariateSpline(y[ids],x,Q[ids])
                alt=measure(x,qfy(h,x)[0],Qfy,h)
                alternatives[f'cubic_y_stride2_phase{phase}']=dict(geometry=alt,delta=delta(alt,base))
            allrows.append(dict(tau=float(d['tau']),native_y=float(h),baseline=base,alternatives=alternatives))
        if tau==20:
            raw=loadmat(ROOT/'runs/steady_20260928_tau20/run_analysis.mat',simplify_cells=True)
            direct=RectBivariateSpline(raw['yFull'],raw['xFull'],raw['rho'])
            qdiff=float(np.max(abs(q-raw['rhoX'])))
            dens=(Q-Qf(0,0)[0,0])-(raw['rho']-direct(0,0)[0,0])
            crosscheck=dict(restored_native_q_max_abs_difference=qdiff,
                            restored_centered_rho_max_abs_difference=float(np.max(abs(dens))),
                            restored_x_max_abs_difference=float(np.max(abs(x-raw['xFull']))),
                            restored_y_max_abs_difference=float(np.max(abs(y-raw['yFull']))))
    return scales,allrows,crosscheck


def audit_history():
    h=loadmat(ROOT/'runs/construction_20260928/geometry_history.mat',simplify_cells=True)
    saved=json.loads((ROOT/'runs/construction_20260928/geometry_transform.json').read_text())['timeConvergence']
    heights=np.geomspace(4,256,60);out=[];method_audit=[]
    for j,(tau,x,y,q,rho) in enumerate(zip(h['tau'],h['x'],h['y'],h['rhoX'],h['rho'])):
        qf=RectBivariateSpline(y,x,q);Qf=RectBivariateSpline(y,x,rho)
        rows=[];geos=[];support=[]
        for Y in heights:
            g=measure(x,qf(Y,x)[0],Qf,Y)
            xx=g['m']+g['W']*(np.sinh(BETA*XI)+g['k']*(np.cosh(BETA*XI)-1))/BETA
            assert xx.min()>=x[0] and xx.max()<=x[-1]
            rows.append(qf.ev(np.full_like(xx,Y),xx)/g['A']);geos.append(g)
            support.append([float(xx.min()),float(xx.max())])
        rows=np.array(rows);ref=np.broadcast_to(master(XI),rows.shape)
        result=metrics(rows,ref)
        excluded=np.min(abs(XI[:,None]-np.array([-1.,0.,1.])),axis=1)>.05
        per=[metrics(row,master(XI))['relative_l2_reference'] for row in rows]
        out.append(dict(tau=float(tau),step=int(h['steps'][j]),
                        samples=rows.size,**result,
                        saved_relative_l2=saved[j]['integerRationalRelativeL2'],
                        reproduction_difference=result['relative_l2_reference']-saved[j]['integerRationalRelativeL2'],
                        outside_anchor_neighborhoods=metrics(rows[:,excluded],ref[:,excluded]),
                        row_relative_l2_range=[float(min(per)),float(max(per))],
                        original_y_nodes_in_range=int(np.sum((y>=4)&(y<=256))),
                        mapped_x_extent=[float(np.min(support)),float(np.max(support))]))
        if j==len(h['tau'])-1:
            for mode,stride,phase in [('pchip',1,0),('quadratic5',1,0),('cubic',2,0),('cubic',2,1),('cubic',4,0),('cubic',4,1)]:
                geom_only=[];consistent=[]
                ids=sample_indices(len(x),stride,phase)
                for Y in heights:
                    row=qf(Y,x)[0];g=measure(x,row,Qf,Y,mode,stride,phase)
                    xx=g['m']+g['W']*(np.sinh(BETA*XI)+g['k']*(np.cosh(BETA*XI)-1))/BETA
                    assert xx.min()>=x[ids[0]] and xx.max()<=x[ids[-1]]
                    interp=(PchipInterpolator(x[ids],row[ids],extrapolate=False) if mode=='pchip'
                            else CubicSpline(x[ids],row[ids],extrapolate=False))
                    geom_only.append(qf.ev(np.full_like(xx,Y),xx)/g['A'])
                    consistent.append(interp(xx)/g['A'])
                method_audit.append(dict(method=f'{mode}_stride{stride}_phase{phase}',
                                         geometry_only=metrics(np.array(geom_only),ref),
                                         geometry_and_profile=metrics(np.array(consistent),ref)))
    return dict(heights=heights.tolist(),xi=XI.tolist(),height_count=60,xi_count=401,
                norm='Unweighted Frobenius norm over 60 log-spaced native heights x 401 uniform xi nodes; primary denominator is broadcast f0, not observed rows.',
                no_extrapolation=True,frames=out,tau20_method_sensitivity=method_audit)


def audit_candidates():
    names=['L4_129x65','L4_257x129','L6_193x97','L8_257x129_continued']
    data={};inventory=[]
    for name in names:
        d=np.load(ROOT/'runs/steady_20260928_tau20'/name/'state.npz')
        x,y,Q,q=[d[k] for k in ['x','y','Q','qx']]
        data[name]=(RectBivariateSpline(y,x,q),RectBivariateSpline(y,x,Q))
        inventory.append(dict(name=name,nx=len(x),ny=len(y),L=float(y[-1]),
                              dx=float(np.diff(x)[0]),dy=float(np.diff(y)[0]),rates=d['rates'].tolist()))
    x=np.linspace(-2,2,129);y=np.linspace(0,1,65)
    pairs=[(names[0],names[1],'same_box_mesh_change'),(names[0],names[2],'same_spacing_box_change'),
           (names[2],names[3],'same_spacing_box_change'),(names[1],names[3],'mesh_and_box_both_change')]
    comparisons=[]
    for a,b,kind in pairs:
        qa,Qa=data[a];qb,Qb=data[b]
        comparisons.append(dict(first=a,reference=b,kind=kind,q=metrics(qa(y,x),qb(y,x)),
                                Q=metrics(Qa(y,x),Qb(y,x)),
                                wall_q=metrics(qa(0,x),qb(0,x))))
    return dict(inventory=inventory,comparison_patch=dict(x=[-2,2],y=[0,1],nx=129,ny=65),comparisons=comparisons)


def summarize_sensitivity(rows):
    result=[]
    for lo,hi in [(0,1),(4,256),(1024,1024)]:
        part=[r for r in rows if lo<=r['native_y']<=hi]
        for category,match in [('pchip',lambda k:k.startswith('pchip')),
                               ('quadratic5',lambda k:k.startswith('quadratic5')),
                               ('x_stride2',lambda k:k.startswith('cubic_stride2')),
                               ('x_stride4',lambda k:k.startswith('cubic_stride4')),
                               ('y_stride2',lambda k:k.startswith('cubic_y_stride2'))]:
            ds=[a['delta'] for r in part for k,a in r['alternatives'].items() if match(k) and a['delta'] is not None]
            keys=sorted(set().union(*(d.keys() for d in ds)))
            result.append(dict(native_y_range=[lo,hi],method=category,comparisons=len(ds),
                               max_absolute_changes={k:max(abs(d[k]) for d in ds if k in d) for k in keys}))
    return result


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    scales,rows,cross=audit_native_frames()
    result=dict(scope='Independent saved-field audit only; interpolation sensitivity is not a rigorous continuum error bound.',
                scales=scales,native_restoration_check=cross,integer_master_history=audit_history(),
                candidate_comparisons=audit_candidates(),native_geometry_sensitivity=rows,
                sensitivity_summary=summarize_sensitivity(rows))
    (OUT/'data_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'history':[{k:r[k] for k in ['tau','relative_l2_reference','relative_l2_observed','linf','reproduction_difference']} for r in result['integer_master_history']['frames']],
                      'sensitivity':result['sensitivity_summary'],
                      'candidate_comparisons':result['candidate_comparisons'],
                      'restoration':cross},indent=2))


if __name__=='__main__':main()
