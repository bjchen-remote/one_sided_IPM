"""Display the position and height of the ridge as the primary skeleton."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parents[2]/'runs/ridge_amplitude_20260928'
report=json.loads((OUT/'amplitude_detail.json').read_text())
frames=report['timeFrames']
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
colors=['#8a6fba','#d38b33','#237fa2']

fig=plt.figure(figsize=(12,4.9),layout='constrained')
ax=fig.add_subplot(1,2,1,projection='3d')
aa=fig.add_subplot(1,2,2)
for f,col in zip(frames,colors):
    rows=f['profileRows']; lab=f"tau={f['tau']:.0f}"
    Y=np.array([r['Ynative'] for r in rows])
    A=np.array([r['A_native'] for r in rows])
    m=np.array([r['mAbsoluteNative'] for r in rows])
    near=Y<=16
    ax.plot(m[near],Y[near],A[near],color=col,lw=2,label=lab)
    ax.scatter(m[0],Y[0],A[0],color=col,s=24)
    pos=Y>0
    aa.loglog(Y[pos],A[pos]/A[0],color=col,lw=2,label=lab)
ax.set(xlabel='ridge position m(Y)',ylabel='native Y',zlabel='peak height A(Y)',
       title='The full ridge: (m(Y), Y, A(Y))')
ax.view_init(elev=24,azim=-64)
ax.legend(fontsize=9)
aa.set(xlabel='native Y',ylabel='A(Y) / A(0)',
       title='Vertical amplitude profile')
for level in [.5,.25,.125,.0625,.03125]:
    aa.axhline(level,color='gray',lw=.6,alpha=.35)
aa.grid(alpha=.18);aa.legend()
fig.savefig(OUT/'ridge_backbone.png',dpi=190);plt.close(fig)

fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
levels=[.75,.5,.25,.125,.0625,.03125]
names=['3/4','1/2','1/4','1/8','1/16','1/32']
for f,col in zip(frames,colors):
    rows=f['profileRows'];lab=f"tau={f['tau']:.0f}"
    Y=np.array([r['Ynative'] for r in rows])
    A=np.array([r['A_native'] for r in rows])
    m=np.array([r['mAbsoluteNative'] for r in rows])
    low=Y<=4
    ax[0,0].plot(Y[low],A[low]/A[0],color=col,lw=2,label=lab)
    slopes=np.array([r['pA']['5'] if r['pA']['5'] is not None else np.nan for r in rows])
    ax[0,1].semilogx(Y[1:],slopes[1:],color=col,lw=2,label=lab)
    cr=f['thresholdCrossings'];yy=[r['Ynative'] for r in cr]
    mm=[r['mAbsoluteNative'] for r in cr]
    ax[1,0].semilogy(names,yy,'o-',color=col,label=lab)
    ax[1,1].loglog(yy,mm,'o-',color=col,label=lab)
ax[0,0].set(xlabel='native Y',ylabel='A(Y) / A(0)',title='Near-wall amplitude loss')
ax[0,1].set(xlabel='native Y',ylabel='-d log A / d log Y',title='Local decay exponent (5-point windows)')
ax[0,1].axhline(.5,color='black',ls='--',lw=1)
ax[1,0].set(xlabel='fraction of wall peak',ylabel='native height Y at that fraction',
            title='Amplitude levels move toward the wall')
ax[1,1].set(xlabel='height of an amplitude level',ylabel='ridge position at that level',
            title='Position and amplitude measured together')
for aa in ax.ravel():aa.grid(alpha=.18);aa.legend(fontsize=9)
fig.savefig(OUT/'ridge_amplitude_scales.png',dpi=190);plt.close(fig)
relations=json.loads((OUT/'relations.json').read_text())['halfHeightTests']
z=np.array(relations['zLog'])
fixed=(1+3*z)**(-.5)
base=np.array(relations['curves'][0]['BLogZ'])
fig,ax=plt.subplots(1,3,figsize=(14,4.3),layout='constrained')
for f,col,c in zip(frames,colors,relations['curves']):
    lab=f"tau={f['tau']:.0f}"
    rows=f['profileRows'];Y=np.array([r['Ynative'] for r in rows])
    B=np.array([r['AOverA0Native'] for r in rows])
    pos=(Y>0)&(Y<=140)
    ax[0].loglog(Y[pos],B[pos],color=col,lw=2,label=lab)
    obs=np.array(c['BLogZ'])
    ax[1].loglog(z,obs,color=col,lw=2,label=lab)
    ax[2].semilogx(z,100*(obs/base-1),color=col,lw=2,label=lab)
ax[1].loglog(z,fixed,'--',color='black',lw=1.4,label='(1 + 3z)^(-1/2)')
ax[1].scatter([1],[.5],c='black',s=28,zorder=5)
ax[0].set(xlabel='native Y',ylabel='A(Y) / A(0)',title='Before vertical rescaling')
ax[1].set(xlabel='z = Y / H',ylabel='A(H z) / A(0)',title='One measured length H per time')
ax[2].set(xlabel='z = Y / H',ylabel='difference from tau=10 (%)',title='Residual shape change after rescaling')
ax[2].axhline(0,c='black',lw=.7)
for aa in ax:aa.grid(alpha=.18);aa.legend(fontsize=9)
fig.savefig(OUT/'ridge_amplitude_collapse.png',dpi=190);plt.close(fig)
print('Saved three ridge amplitude figures')
