"""Measure a proposed single-invariant slope collapse in saved IPM fields."""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
from scipy.optimize import brentq
from construct_profile import load_profile

root=Path(__file__).resolve().parents[2]
folder=root/'runs/steady_20260928_tau20'
out=root/'runs/normal_flow_closure_20260928'
out.mkdir(exist_ok=True)
results=[]
for name in ['tau20','L4_129x65','L4_257x129','L6_193x97','L8_257x129_continued']:
    p=load_profile(folder,name)
    lam,mu,v,_=p['rates'];alpha=mu/lam;shift=v/mu
    Q=RectBivariateSpline(p['y'],p['x'],p['Q'])
    levels=np.linspace(-.1,1.15,251)
    heights=np.linspace(0,1,17)
    records=[]
    stationary=[]
    for height in heights:
        xs=np.array([brentq(lambda x:Q(height,x)[0,0]-(s-shift),p['x'][0],p['x'][-1])
                     for s in levels])
        ys=np.full_like(xs,height)
        q=Q.ev(ys,xs,dx=0,dy=1)
        assert np.all(q>0)
        K=-Q.ev(ys,xs,dx=1,dy=0)/q
        t=levels*q**(alpha/(1-alpha))
        if np.any(np.diff(t)<=0):
            raise ValueError(f'Invariant not injective in core: {name}, {height}')
        records.append((t,K,q))
        xx=brentq(lambda x:Q(height,x)[0,0]+shift,p['x'][0],p['x'][-1])
        qq=float(Q(height,xx,dx=0,dy=1)[0,0])
        KK=-float(Q(height,xx,dx=1,dy=0)[0,0])/qq
        qqs=float(Q(height,xx,dx=0,dy=2)[0,0])/qq
        curvature=-float(Q(height,xx,dx=2,dy=0)[0,0]
                         +2*KK*Q(height,xx,dx=1,dy=1)[0,0]
                         +KK*KK*Q(height,xx,dx=0,dy=2)[0,0])/qq
        stationary.append({'height':float(height),'x':float(xx),'q':qq,'K':KK,'q_s':qqs,
                           'X_yy_spline':curvature})
    lo=max(r[0][0] for r in records);hi=min(r[0][-1] for r in records)
    tt=np.linspace(lo,hi,201)
    profiles=np.array([PchipInterpolator(t,K)(tt) for t,K,q in records])
    qprofiles=np.array([PchipInterpolator(t,q)(tt) for t,K,q in records])
    errors=profiles[1:]-profiles[0]
    dq=qprofiles[1:]-qprofiles[0]
    train=heights[1:]<=.5
    hold=~train
    slope_a=float(np.sum(errors[train]*dq[train])/np.sum(dq[train]**2))
    corrected=errors-slope_a*dq
    wall=stationary[0]
    k_exact=wall['q']/(2*(lam-mu))
    predicted=wall['q_s']*(mu*(2*lam-mu)+(4*lam*lam-6*lam*mu-mu*mu)*k_exact*k_exact)/((2*lam-mu)*(4*lam-3*mu))
    results.append({'name':name,'alpha':float(alpha),'density_shift':float(shift),
                    'source_kind':'finite-time trajectory' if name=='tau20' else 'conditional finite-box root',
                    'height_interval':[0,1],'s_interval':[float(levels[0]),float(levels[-1])],
                    'shared_t_interval':[float(lo),float(hi)],
                    'slope_relative_l2_against_wall':float(np.linalg.norm(errors)/np.linalg.norm(profiles[1:])),
                    'slope_linf_against_wall':float(abs(errors).max()),
                    'wall_curvature_predicted_by_steady_PDE':float(predicted),
                    'wall_curvature_measured_by_spline':wall['X_yy_spline'],
                    'one_constant_q_correction':{
                        'ansatz':'K=kappa(s*q**(alpha/(1-alpha))) + a*q',
                        'a_fit':slope_a,'mu_times_a_fit':float(mu*slope_a),
                        'training_height_interval':[.0625,.5],
                        'heldout_height_interval':[.5625,1.],
                        'train_relative_l2':float(np.linalg.norm(corrected[train])/np.linalg.norm(profiles[1:][train])),
                        'heldout_relative_l2':float(np.linalg.norm(corrected[hold])/np.linalg.norm(profiles[1:][hold])),
                        'uncorrected_heldout_relative_l2':float(np.linalg.norm(errors[hold])/np.linalg.norm(profiles[1:][hold])),
                        'heldout_linf':float(abs(corrected[hold]).max()),
                        'one_third_stationary_necessary_a':float(-(5+17*k_exact*k_exact)/(18*k_exact*wall['q'])),
                        'scope':'One constant fit and an entire measured wall curve; a data-selected relation, not an analytic master curve or complete PDE closure.'},
                    'stationary':stationary,
                    'caveat':'Spline derivatives on saved grids; this is a geometry diagnostic, not a continuum proof.'})
(out/'data_probe.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps([{k:v for k,v in r.items() if k!='stationary'}|{'stationary_wall':r['stationary'][0],
                                                             'stationary_height1':r['stationary'][-1]}
                 for r in results],indent=2))
