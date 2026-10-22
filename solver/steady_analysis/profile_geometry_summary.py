"""Derived data checks and a geometry-parameter figure; no equation solving."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from profile_shape_atlas import f0

OUT=Path(__file__).resolve().parents[2]/'runs/profile_geometry_20260928'
data=json.loads((OUT/'geometry_detail.json').read_text())
beta=1.5
left=brentq(lambda z:f0(z)-.75,-1,0)
right=brentq(lambda z:f0(z)-.75,0,1)
checks=[]
slopes=[]
for frame in data['frames']:
    rows=frame['nativeYProfiles']
    for r in rows:
        k,W=r['k'],r['Wnative']
        G=lambda z:(np.sinh(beta*z)+k*(np.cosh(beta*z)-1))/beta
        lev=r['levels']['0.75']
        predicted=[-W*G(left),W*G(right)]
        observed=[lev['widthLeftNative'],lev['widthRightNative']]
        checks.append(dict(tau=frame['tau'],Y=r['Ynative'],
            leftRelativeDifference=float(predicted[0]/observed[0]-1),
            rightRelativeDifference=float(predicted[1]/observed[1]-1)))
    positive=[r for r in rows if r['Ynative']>0]
    for j in range(1,len(positive)-1):
        lo,center,hi=positive[j-1:j+2]
        logspan=np.log(hi['Ynative']/lo['Ynative'])
        def slope(key): return float(np.log(hi[key]/lo[key])/logspan)
        slopes.append(dict(tau=frame['tau'],Y=center['Ynative'],
            window=[lo['Ynative'],hi['Ynative']],
            mSlope=slope('mAbsoluteNative'),WSlope=slope('Wnative'),
            inferredAlphaFromA=1+slope('A_native'),
            AWslope=slope('AW_native'),cSlope=slope('c_origin_rhoMinusRho00')))
report=dict(
    note='Three-quarter widths are not calibration anchors; compared with fixed f0 using measured peak/half-width geometry.',
    xiAtThreeQuarter=[left,right],threeQuarterChecks=checks,
    localSlopes=slopes,
    slopeWarning='Centered secants on logarithmic native heights. These different quantities need not share a finite-height exponent; no asymptotic alpha is selected.')
(OUT/'derived_geometry_checks.json').write_text(json.dumps(report,indent=2)+'\n')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,3,figsize=(13,7.3),constrained_layout=True)
for frame,color in zip(data['frames'],['#8a6fba','#d38b33','#237fa2']):
    rows=[r for r in frame['nativeYProfiles'] if r['Ynative']>0]
    y=np.array([r['Ynative'] for r in rows])
    vectors=[
        [r['mOverYAbsolute'] for r in rows],
        [r['WnativeOverY'] for r in rows],
        [r['k'] for r in rows],
        [4*r['cOriginOverYA'] for r in rows],
        [r['A_native'] for r in rows],
        [r['Wnative']/r['curvatureWidthNative'] for r in rows]]
    for aa,values in zip(ax.ravel(),vectors):
        aa.semilogx(y,values,'o-',color=color,label=f"tau={frame['tau']:.0f}",markersize=4)
for aa,target in zip([ax[0,0],ax[0,1],ax[0,2],ax[1,0],ax[1,2]],
                     [1,.5,.75,1,np.sqrt(7/8)]):
    aa.axhline(target,color='black',ls='--',lw=1)
labels=[
    ('Peak-line ratio','m / Y'),('Width-to-height ratio','W / Y'),
    ('Geometric asymmetry','k'),('Quarter-density coincidence','4 c / (Y A)'),
    ('Peak gradient amplitude','A'),('Width / curvature width','W / w_curv')]
for aa,(title,ylabel) in zip(ax.ravel(),labels):
    aa.set(title=title,xlabel='native Y',ylabel=ylabel)
    aa.axvspan(4,256,color='gray',alpha=.08);aa.grid(alpha=.18)
ax[1,1].set_yscale('log')
ax[0,0].legend(fontsize=8)
fig.suptitle('Measured geometry evolves even when the normalized core shape nearly collapses')
fig.savefig(OUT/'profile_geometry_parameters.png',dpi=190)
plt.close(fig)
print(json.dumps({'threeQuarterXi':[left,right],
    'tau20ThreeQuarter':[r for r in checks if r['tau']>19 and r['Y'] in [4,16,64,256]],
    'tau20Slopes':[r for r in slopes if r['tau']>19 and r['Y'] in [4,16,64,256]]},indent=2))
