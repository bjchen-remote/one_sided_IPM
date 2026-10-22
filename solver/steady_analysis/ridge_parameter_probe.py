"""Small, predeclared amplitude-coordinate skeleton families; no PDE fitting."""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator, RectBivariateSpline
from profile_geometry_measure import load_frame, measure_profile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/ridge_parameter_20260928'
EXPONENTS = [1.5, 5/3, 2.0]


def basis(v, p, density=False):
    if density:
        return (v**(p-1)-1)/(2**(p-1)-1), 2*(1-1/v)
    return (v**p-1)/(2**p-1), v-1


def errors(pred, observed, mask):
    d = pred[mask]-observed[mask]
    return dict(relativeL2=float(np.linalg.norm(d)/np.linalg.norm(observed[mask])),
                relativeRMS=float(np.sqrt(np.mean((d/observed[mask])**2))))


def extract():
    saved = json.loads((ROOT/'runs/ridge_amplitude_20260928/amplitude_detail.json').read_text())
    requested = np.unique(np.r_[1., 2., 4., 8., 16., np.geomspace(1.01,16.,121)])
    frames = []
    for old in saved['timeFrames']:
        tau = int(round(old['tau']))
        f = load_frame(tau)
        aq = np.array([r['AOverA0Native'] for r in old['profileRows']])
        ay = np.array([r['Ynative'] for r in old['profileRows']])
        assert np.all(np.diff(aq)<0)
        yofv = PchipInterpolator(1/aq, ay, extrapolate=False)
        YY = yofv(requested)
        if np.any(~np.isfinite(YY)):
            raise ValueError('Amplitude level would extrapolate beyond saved curves')
        q = RectBivariateSpline(f['y'],f['x'],f['q'])
        Q = RectBivariateSpline(f['y'],f['x'],f['Q'])
        rows = [measure_profile(f,float(y/f['ell']),q,Q) for y in YY]
        A = np.array([r['AInner']*f['P'] for r in rows])
        m = np.array([r['mInner']*f['ell'] for r in rows])
        C = np.array([r['Qpeak']*f['P']*f['ell'] for r in rows])
        W = np.array([.5*(r['levels']['0.5']['widthLeftNative']+r['levels']['0.5']['widthRightNative']) for r in rows])
        Wc = np.array([r['curvatureWidthInner']*f['ell'] for r in rows])
        asym = np.array([r['halfAsymmetry'] for r in rows])
        v = A[0]/A
        assert np.all(np.diff(v)>0)
        # Anchors are measured (interpolated) at B=1/2, not fitted model parameters.
        raw = {'Y':YY,'m':m-m[0],'W':W-W[0],'Wc':Wc-Wc[0],'C':C-C[0]}
        anchors = {key:float(PchipInterpolator(v,val)(2.)) for key,val in raw.items()}
        normalized = {key:val/anchors[key] for key,val in raw.items()}
        fits = []
        for p in EXPONENTS:
            parts = {}
            for key,val in normalized.items():
                b1,b2 = basis(v,p,key=='C')
                train = (v>1)&(v<=4)
                weight = 1/(1+abs(val[train]))
                aa = (b1-b2)[train]*weight
                bb = (val-b2)[train]*weight
                theta = float(aa@bb/(aa@aa))
                pred = theta*b1+(1-theta)*b2
                parts[key] = dict(theta=theta,
                    training=errors(pred,val,train),
                    heldout=errors(pred,val,(v>4)&(v<=8)),
                    outerDiagnostic=errors(pred,val,(v>8)&(v<=16)),
                    sampledPositive=bool(np.all(pred[v>1]>0)),
                    leadingCoefficientPositive=theta>0)
            fits.append(dict(p=p,formalHomogeneousAlpha=1-1/p,parts=parts))
        xi = np.linspace(-1.5,1.5,241)
        shape = np.array([q.ev(np.full_like(xi,YY[i]/f['ell']),
                              (m[i]+W[i]*xi)/f['ell'])/(A[i]/f['P']) for i in range(len(v))])
        jhalf = int(np.argmin(abs(v-2)))
        shape_change = np.linalg.norm(shape-shape[jhalf],axis=1)/np.linalg.norm(shape[jhalf])
        frame = dict(tau=tau,ell=f['ell'],A0=float(A[0]),m0=float(m[0]),W0=float(W[0]),Wc0=float(Wc[0]),C0=float(C[0]),
            anchors=anchors,rows=[dict(v=float(v[i]),Y=float(YY[i]),A=float(A[i]),
                m=float(m[i]),W=float(W[i]),Wc=float(Wc[i]),C=float(C[i]),skew=float(asym[i]),
                affineShapeRelativeChangeFromHalfAmplitude=float(shape_change[i])) for i in range(len(v))],
            models=fits,
            affineShapeChange=dict(trainingMaximum=float(shape_change[v<=4].max()),
                heldoutMaximum=float(shape_change[(v>4)&(v<=8)].max())))
        integer_ridge=(7*v**1.5-8*v+1)/(14*np.sqrt(2)-15)
        frame['integerRidgeCandidate']=dict(
            formula='m=m0+Mstar*(7*t^3-8*t^2+1), t=sqrt(A0/A), Mstar=(mHalf-m0)/(14*sqrt(2)-15)',
            Mstar=float(anchors['m']/(14*np.sqrt(2)-15)),
            training=errors(integer_ridge,normalized['m'],(v>1)&(v<=4)),
            heldout=errors(integer_ridge,normalized['m'],(v>4)&(v<=8)),
            outerDiagnostic=errors(integer_ridge,normalized['m'],(v>8)&(v<=16)),
            status='Post-hoc small-integer compression of the fitted training coefficients; not an exact-constant identification or blind validation.')
        frames.append(frame)
        np.savez_compressed(OUT/f'tau{tau}.npz',v=v,Y=YY,A=A,m=m,W=W,Wc=Wc,C=C,
                            skew=asym,xi=xi,shape=shape)
        print('measured amplitude skeleton tau',tau,flush=True)
    # Test common normalized laws: tau10 coefficients, only measured scalar anchors
    # are updated in tau20. No height-dependent correction is introduced.
    transfer = []
    target = frames[-1]
    v = np.array([r['v'] for r in target['rows']])
    for model in frames[0]['models']:
        parts = {}
        for key in ['Y','m','W','Wc','C']:
            vals = np.array([r[key] for r in target['rows']])
            if key!='Y': vals = vals-target[key+'0']
            vals = vals/target['anchors'][key]
            b1,b2 = basis(v,model['p'],key=='C')
            a = model['parts'][key]['theta']
            pred = a*b1+(1-a)*b2
            parts[key] = dict(trainingBand=errors(pred,vals,(v>1)&(v<=4)),
                             heldoutBand=errors(pred,vals,(v>4)&(v<=8)))
        transfer.append(dict(p=model['p'],parts=parts))
    result = dict(coordinate='v=A(0)/A(Y); native X,Y; C=rho at ridge relative to saved wall density origin',
        halfAmplitudeAnchors='Y,m-m0,W-W0,C-C0 evaluated at v=2',
        objective='One linear coefficient per curve, anchored at v=1 and v=2; weighted least squares on 1<v<=4',
        geometryFamily='a*(v^p-1)/(2^p-1)+(1-a)*(v-1)',
        densityFamily='a*(v^(p-1)-1)/(2^(p-1)-1)+(1-a)*2*(1-1/v)',
        frames=frames,transferTau10To20=transfer,
        caveats=['Finite-time trajectory data, not converged global steady states.',
            'Half-amplitude anchors are identities of normalization, not independent fit evidence.',
            'Candidate exponents are predeclared diagnostics, not selected asymptotic exponents.',
            'A finite skeleton of ridge curves does not imply closure of the transverse mother shape.'])
    (OUT/'amplitude_parameter.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.5,4),constrained_layout=True)
    for frame in frames:
        v=np.array([r['v'] for r in frame['rows']])
        val=(np.array([r['m'] for r in frame['rows']])-frame['m0'])/frame['anchors']['m']
        mask=(v>1)&(v<=8)
        model=(7*v**1.5-8*v+1)/(14*np.sqrt(2)-15)
        line,=axes[0].plot(np.sqrt(v[mask]),val[mask],label=f"tau {frame['tau']}")
        axes[1].plot(np.sqrt(v[mask]),100*(model[mask]/val[mask]-1),color=line.get_color())
    vv=np.linspace(1,8,401)
    axes[0].plot(np.sqrt(vv),(7*vv**1.5-8*vv+1)/(14*np.sqrt(2)-15),'k--',label='integer cubic candidate')
    for ax in axes:
        ax.axvline(2,color='gray',ls=':',lw=1)
        ax.set_xlabel('t = sqrt(A(0) / A)')
        ax.grid(alpha=.2)
    axes[0].set_ylabel('ridge displacement / half-amplitude displacement')
    axes[1].set_ylabel('candidate relative displacement error (%)')
    axes[0].legend(fontsize=8)
    axes[0].set_title('One amplitude parameter for the ridge')
    axes[1].set_title('Dotted line: start of amplitude holdout')
    fig.savefig(OUT/'ridge_cubic_candidate.png',dpi=170)
    plt.close(fig)


if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    extract()
