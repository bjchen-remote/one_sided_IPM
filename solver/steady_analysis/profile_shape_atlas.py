"""Data-only profile alignment, residual structure, and low-rank diagnostics."""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline, RectBivariateSpline
from geometry_measure import trajectory_frame, interp_fields, peak_measure

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/profile_geometry_20260928'
OLD=ROOT/'runs/geometry_loop_20260928'
HEIGHTS=np.array([0,.25,.5,1,2,4,8,16,32,64,128,256,512,1024.])
XI=np.linspace(-1.5,1.5,601)
OUTER_XI=np.linspace(-3,3,601)
BET=1.5


def f0(x):
    return 32/(x**6-3*x**5+3*x**4+3*x**3+28*x*x+32)


def relative(pred,obs):
    return float(np.linalg.norm(pred-obs)/np.linalg.norm(obs))


def profile_data():
    result={}
    for tau in [10,15,20]:
        p=trajectory_frame(tau)
        meta=loadmat(OLD/f'frame_tau{tau}.mat',simplify_cells=True)
        field,density=interp_fields(p)
        ell=float(meta['ell']); a=float(meta['a']); amp=float(meta['P'])
        rows=[]; outer=[]; shapes={name:[] for name in ['affine_halfwidth','piecewise_halfwidth','sinh']}
        for Y in HEIGHTS:
            geo=peak_measure(p,float(Y/ell),field,density)
            m,A,W,k=geo['m'],geo['A'],geo['halfWidthGeometry']['W'],geo['halfWidthGeometry']['skew']
            left,right=geo['halfWidthGeometry']['left'],geo['halfWidthGeometry']['right']
            physical={
                'affine_halfwidth':m+(left+right)/2*XI,
                'piecewise_halfwidth':m+np.where(XI<0,left,right)*XI,
                'sinh':m+W*(np.sinh(BET*XI)+k*(np.cosh(BET*XI)-1))/BET}
            for name,xx in physical.items():
                if not np.all((xx>=p['x'][0])&(xx<=p['x'][-1])):
                    raise ValueError(f'Outside saved field: tau={tau},Y={Y},map={name}')
                shapes[name].append(field(float(Y/ell),xx)[0]/A)
            xx=m+W*(np.sinh(BET*OUTER_XI)+k*(np.cosh(BET*OUTER_XI)-1))/BET
            if not np.all((xx>=p['x'][0])&(xx<=p['x'][-1])):
                raise ValueError('Outside saved field in transverse outer diagnostic')
            outer.append(field(float(Y/ell),xx)[0]/A)
            rows.append(dict(Y=float(Y),m=float(a+ell*m),A=float(amp*A),
                W=float(ell*W),k=float(k),left=float(ell*left),right=float(ell*right),
                curvatureWidth=float(ell*geo['w']),WOverCurvatureWidth=float(W/geo['w'])))
        result[tau]=dict(rows=rows,shapes={k:np.asarray(v) for k,v in shapes.items()},
                         outer=np.asarray(outer),a=a,ell=ell,P=amp)
    return result


def analyze(data):
    middle=(HEIGHTS>=4)&(HEIGHTS<=256)
    selected=(HEIGHTS==4)|(HEIGHTS==16)|(HEIGHTS==64)|(HEIGHTS==256)
    hold_height=(HEIGHTS==8)|(HEIGHTS==32)|(HEIGHTS==128)
    names=list(data[20]['shapes'])
    alignment={}
    for name in names:
        mat=data[20]['shapes'][name][middle]
        master=mat.mean(axis=0)
        alignment[name]=dict(rowCount=len(mat),meanCurveSpread=relative(np.broadcast_to(master,mat.shape),mat),
             maxPointwiseStd=float(mat.std(axis=0).max()),
             note='Equal weighting of sampled heights and canonical coordinate; physical bands differ between maps.')
    allrows=np.concatenate([data[t]['shapes']['sinh'] for t in [10,15,20]])
    labels=[dict(tau=t,Y=float(y)) for t in [10,15,20] for y in HEIGHTS]
    scores=[]
    even_f=32/(XI**6+3*XI**4+28*XI**2+32)
    base=XI**2*(XI**2-1)
    rational_basis=np.column_stack([base,base*XI,base*XI**2])
    for label,F in zip(labels,allrows):
        odd_inverse=(1/F-1/F[::-1])/2
        odd_basis=XI**3*(1-XI**2)
        central=(abs(XI)>=.2)&(abs(XI)<=1.5)&(abs(abs(XI)-1)>.08)
        c=float(odd_basis[central]@odd_inverse[central]/(odd_basis[central]@odd_basis[central]))
        odd_res=relative(c*odd_basis[central],odd_inverse[central])
        fits={}
        for cutoff in [1.2,1.5]:
            mask=abs(XI)<=cutoff
            coef=np.linalg.lstsq(rational_basis[mask],(1/F-1-XI**2)[mask],rcond=None)[0]
            denominator=1+XI**2+rational_basis@coef
            prediction=1/denominator
            fits[str(cutoff)]=dict(abc=coef.tolist(),fitVariable='reciprocal profile',
                fullProfileRelativeL2=relative(prediction,F),
                outerHoldoutRelativeL2=relative(prediction[~mask],F[~mask]) if cutoff<1.5 else None,
                minimumCoreDenominator=float(denominator.min()))
        scores.append(dict(**label,integerF0RelativeL2=relative(f0(XI),F),
            evenIntegerRelativeL2=relative(even_f,F),lorentzRelativeL2=relative(1/(1+XI*XI),F),
            maximumAbsoluteF0Error=float(max(abs(f0(XI)-F))),
            oddReciprocalCoefficient=c,integerOddCoefficient=3/32,
            oddReciprocalRelativeResidual=odd_res,
            reciprocalThreeCoefficients=fits,
            values={str(x):float(np.interp(x,XI,F)) for x in [-1.5,-1.25,-.75,-.5,0,.5,.75,1.25,1.5]}))
    # Basis trained only on tau10/tau15 and alternating logarithmic heights.
    train=np.concatenate([data[t]['shapes']['sinh'][selected] for t in [10,15]])
    mean=train.mean(axis=0)
    _,s,vt=np.linalg.svd(train-mean,full_matrices=False)
    hold={}
    for label,mat in [
        ('tau10_15_unseen_heights',np.concatenate([data[t]['shapes']['sinh'][hold_height] for t in [10,15]])),
        ('tau20_all_middle_heights',data[20]['shapes']['sinh'][middle]),
        ('tau20_near_wall',data[20]['shapes']['sinh'][HEIGHTS<=1]),
        ('tau20_outer_diagnostic',data[20]['shapes']['sinh'][HEIGHTS>=512])]:
        delta=mat-mean
        hold[label]=dict(rows=len(mat),meanOnlyRelativeL2=relative(np.broadcast_to(mean,mat.shape),mat))
        for rank in [1,2]:
            pred=mean+(delta@vt[:rank].T)@vt[:rank]
            hold[label][f'meanPlus{rank}ModesRelativeL2']=relative(pred,mat)
            hold[label][f'fractionResidualEnergyCaptured{rank}']=float(1-np.linalg.norm(pred-mat)**2/np.linalg.norm(delta)**2)
    time={}
    regional_time={}
    regions={'left_outer':(OUTER_XI>=-3)&(OUTER_XI<=-1.5),
             'core':abs(OUTER_XI)<=1,
             'right_outer':(OUTER_XI>=1.5)&(OUTER_XI<=3)}
    for lo,hi in [(10,15),(15,20),(10,20)]:
        time[f'{lo}_to_{hi}']=[dict(Y=float(y),
            relativeL2=relative(data[hi]['shapes']['sinh'][i],data[lo]['shapes']['sinh'][i]))
            for i,y in enumerate(HEIGHTS)]
        regional_time[f'{lo}_to_{hi}']=[dict(Y=float(y),**{
            name:relative(data[hi]['outer'][i,mask],data[lo]['outer'][i,mask])
            for name,mask in regions.items()}) for i,y in enumerate(HEIGHTS)]
    report=dict(
        scope='Saved native trajectory data only. No PDE or symbolic closure solve.',
        coordinates='Each frame uses its own X=a+ell*x,Y=ell*y,rho_X=P*q; heights in this report are fixed native Y.',
        weighting='Uniform xi quadrature and equal sampled height rows; these are shape errors, not physical-area errors.',
        source='runs/geometry_loop_20260928/frame_tau10.mat, frame_tau15.mat, frame_tau20.mat',
        xiSupport=[-1.5,1.5],heightRange=HEIGHTS.tolist(),
        forcedAnchors=['q/A=1 at the peak','sinh and piecewise maps place half heights at xi=+/-1'],
        alignmentTau20Middle=alignment,shapeMetrics=scores,temporalShapeDrift=time,
        regionalTemporalDrift=regional_time,
        lowRank=dict(trainingTaus=[10,15],trainingHeights=HEIGHTS[selected].tolist(),
            centeredSingularValues=s.tolist(),
            centeredEnergyFractions=(s*s/sum(s*s)).tolist(),holdout=hold,
            restriction='Transverse modes are frozen from training; holdout modal amplitudes are measured by projection. This is compression, not forecasting.'),
        geometry={str(t):data[t]['rows'] for t in data},
        simpleShapes=dict(integer='32/(xi^6-3xi^5+3xi^4+3xi^3+28xi^2+32)',
            evenInteger='32/(xi^6+3xi^4+28xi^2+32)',lorentz='1/(1+xi^2)',
            oddReciprocal='(1/f(xi)-1/f(-xi))/2 = (3/32)*xi^3*(1-xi^2) for the proposed integer profile',
            anchoredThreeCoefficientForm='1/f = 1+xi^2+xi^2*(xi^2-1)*(a+b*xi+c*xi^2)',
            integerABC=[1/8,-3/32,1/32],
            interpretation='Zero odd difference at peak/half-height anchors is imposed by normalization. Its amplitude and between-anchor pattern are measured independently.'))
    arrays=dict(xi=XI,outerXi=OUTER_XI,heights=HEIGHTS,trainingMean=mean,residualModes=vt[:3])
    for t in data:
        for name,F in data[t]['shapes'].items(): arrays[f'tau{t}_{name}']=F
        arrays[f'tau{t}_outer']=data[t]['outer']
    return report,arrays


def plots(data,report,arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    colors=plt.cm.viridis(np.linspace(.1,.9,7))
    mid=np.where((HEIGHTS>=4)&(HEIGHTS<=256))[0]
    fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for idx,color in zip(mid,colors):
        F=data[20]['shapes']['sinh'][idx]
        ax[0,0].plot(XI,F,color=color,label=f'Y={HEIGHTS[idx]:g}')
        ax[0,1].plot(XI,100*(F-f0(XI)),color=color)
        pos=XI>=0
        odd=(1/F-1/F[::-1])/2
        ax[1,0].plot(XI[pos],odd[pos],color=color)
    ax[0,0].plot(XI,f0(XI),'k--',lw=1.4,label='integer profile')
    ax[0,0].set(title='One shape across native heights (tau=20)',xlabel='xi',ylabel='gradient / peak')
    ax[0,0].legend(fontsize=8,ncol=2)
    ax[0,1].set(title='Remaining shape difference',xlabel='xi',ylabel='100 x (measured - integer)')
    xp=XI[XI>=0]
    ax[1,0].plot(xp,3/32*xp**3*(1-xp*xp),'k--',lw=1.5)
    ax[1,0].set(title='Residual asymmetry after geometric alignment',xlabel='positive xi',ylabel='odd part of reciprocal profile')
    for t,color in zip([10,15,20],['#8a6fba','#d38b33','#237fa2']):
        rows=[r for r in report['shapeMetrics'] if r['tau']==t and r['Y']>0]
        ax[1,1].semilogx([r['Y'] for r in rows],[100*r['integerF0RelativeL2'] for r in rows],'o-',color=color,label=f'tau={t}')
    ax[1,1].axvspan(4,256,color='gray',alpha=.1)
    ax[1,1].set(title='Fixed integer shape: time and height dependence',xlabel='native Y',ylabel='relative shape error (%)')
    ax[1,1].legend()
    for aa in ax.ravel(): aa.grid(alpha=.18)
    fig.savefig(OUT/'profile_shape_atlas.png',dpi=190); plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    for name,style in zip(['affine_halfwidth','piecewise_halfwidth','sinh'],['-',':','--']):
        mat=data[20]['shapes'][name][mid]
        ax[0].plot(XI,np.std(mat,axis=0)*100,style,label=name.replace('_',' '))
    ax[0].set(title='Cross-height spread after alignment',xlabel='canonical coordinate',ylabel='pointwise standard deviation (%)')
    ax[0].legend(fontsize=8)
    for j in range(2):ax[1].plot(XI,arrays['residualModes'][j],label=f'centered mode {j+1}')
    ax[1].set(title='Frozen modes from tau=10,15 only',xlabel='xi',ylabel='unit-norm shape mode')
    ax[1].legend(fontsize=8)
    h=report['lowRank']['holdout']
    groups=list(h); x=np.arange(len(groups)); width=.24
    for j,(key,label) in enumerate([('meanOnlyRelativeL2','mean'),('meanPlus1ModesRelativeL2','+ 1 mode'),('meanPlus2ModesRelativeL2','+ 2 modes')]):
        ax[2].bar(x+(j-1)*width,[100*h[g][key] for g in groups],width,label=label)
    ax[2].set_xticks(x,['height\nholdout','tau20\nmiddle','tau20\nnear wall','tau20\nouter'])
    ax[2].set(title='Compression on held-out profiles',ylabel='relative shape error (%)')
    ax[2].legend(fontsize=8)
    for aa in ax:aa.grid(alpha=.18)
    fig.savefig(OUT/'profile_compression.png',dpi=190);plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(13,3.8),constrained_layout=True)
    for t,color in zip([10,15,20],['#8a6fba','#d38b33','#237fa2']):
        rows=[r for r in report['shapeMetrics'] if r['tau']==t and r['Y']>=4]
        for j in range(3):
            ax[j].semilogx([r['Y'] for r in rows],
                [r['reciprocalThreeCoefficients']['1.5']['abc'][j] for r in rows],
                'o-',color=color,label=f'tau={t}')
    for j,target in enumerate([1/8,-3/32,1/32]):
        ax[j].axhline(target,color='black',ls='--',label=['1/8','-3/32','1/32'][j])
        ax[j].axvspan(4,256,color='gray',alpha=.08)
        ax[j].set(title=['even coefficient a','odd coefficient b','even coefficient c'][j],
                  xlabel='native Y')
        ax[j].grid(alpha=.18);ax[j].legend(fontsize=8)
    fig.suptitle('Reciprocal shape: 1/f = 1 + xi^2 + xi^2 (xi^2-1) (a+b xi+c xi^2)')
    fig.savefig(OUT/'profile_reciprocal_coefficients.png',dpi=190);plt.close(fig)
    # Physical near-field fan: real X,Y, with each row normalized by its peak.
    p=trajectory_frame(20);field,density=interp_fields(p)
    meta=loadmat(OLD/'frame_tau20.mat',simplify_cells=True)
    a,ell=map(float,[meta['a'],meta['ell']])
    y=np.linspace(0,48,161); x=np.linspace(-30,100,421)
    peak=[peak_measure(p,float(v/ell),field,density) for v in y]
    A=np.array([r['A'] for r in peak]);m=a+ell*np.array([r['m'] for r in peak])
    F=field(y/ell,(x-a)/ell)/A[:,None]
    fig,ax=plt.subplots(figsize=(8,5.5),constrained_layout=True)
    im=ax.pcolormesh(x,y,F,shading='auto',cmap='viridis',vmin=0,vmax=1)
    co=ax.contour(x,y,F,levels=[.25,.5,.75],colors='white',linewidths=.9)
    ax.clabel(co,fmt='%.2f',fontsize=9)
    ax.plot(m,y,'r--',lw=1.4,label='measured peak line')
    ax.plot(y,y,'r:',alpha=.7,label='X = Y (reference)')
    ax.set(xlabel='native X',ylabel='native Y',title='Profile geometry: widening and right-skew (tau=20)')
    ax.legend(loc='upper left',fontsize=8)
    fig.colorbar(im,ax=ax,label='rho_X(X,Y) / max_X rho_X(X,Y)')
    fig.savefig(OUT/'profile_geometry_map.png',dpi=190);plt.close(fig)


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    data=profile_data()
    report,arrays=analyze(data)
    (OUT/'shape_atlas.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(OUT/'shape_atlas.npz',**arrays)
    plots(data,report,arrays)
    print(json.dumps({'alignment':report['alignmentTau20Middle'],
        'centeredEnergyFractions':report['lowRank']['centeredEnergyFractions'],
        'holdout':report['lowRank']['holdout'],
        'tau20Representative':[r for r in report['shapeMetrics'] if r['tau']==20 and r['Y'] in [0,4,16,64,256]]},indent=2))
