"""Independent saved-field audit of ridge amplitude, coordinates, and wall jets."""
from pathlib import Path
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline,PchipInterpolator,RectBivariateSpline

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_amplitude_20260928'
YTEST=np.array([0.,.05,.1,.25,.5,1.,2.,4.,8.,16.])


def peak(x,q,mode='cubic',stride=1,phase=0):
    ids=np.unique(np.r_[0,np.arange(phase,len(x),stride),len(x)-1]).astype(int)
    xx,qq=x[ids],q[ids]
    if mode=='pchip':
        i=np.argmax(qq);return float(xx[i]),float(qq[i])
    cs=CubicSpline(xx,qq,extrapolate=False)
    roots=cs.derivative().roots(extrapolate=False)
    points=np.r_[xx[0],roots[(roots>=xx[0])&(roots<=xx[-1])],xx[-1]]
    k=np.argmax(cs(points));return float(points[k]),float(cs(points[k]))


def relative(a,b):
    return dict(relative_l2=float(np.linalg.norm(a-b)/np.linalg.norm(b)),
                normalized_amplitude_linf=float(np.max(abs(a-b))),
                max_relative=float(np.max(abs(a/b-1))))


def main():
    output=[];native_curves={};inner_curves={}
    grid=np.linspace(0,16,401)
    for tau in [10,15,20]:
        f=loadmat(ROOT/f'runs/geometry_loop_20260928/frame_tau{tau}.mat',simplify_cells=True)
        x=f['a']+f['ell']*f['x'];y=f['ell']*f['y'];q=f['P']*f['q']
        field=RectBivariateSpline(y,x,q)
        def get(Y,mode='cubic',stride=1,phase=0):
            return peak(x,field(Y,x)[0],mode,stride,phase)
        m0,A0=get(0)
        original_A=np.array([peak(x,row)[1] for row in q])
        def kinematic(ff,Y):
            mm,aa=peak(x,ff(Y,x)[0])
            qyy=float(ff.ev(Y,mm,dx=2,dy=0))
            correction=float(-ff.ev(Y,mm,dx=1,dy=1)**2/ff.ev(Y,mm,dx=0,dy=2))
            return dict(A1=float(ff.ev(Y,mm,dx=1,dy=0)),A2=qyy+correction,
                        qyy=qyy,ridge_motion_correction=correction)
        geometry=[get(float(Y)) for Y in YTEST]
        grid_A=np.array([get(float(Y))[1]/A0 for Y in grid])
        inner_A=np.array([get(float(f['ell']*Y))[1]/A0 for Y in grid])
        native_curves[tau]=grid_A;inner_curves[tau]=inner_A
        native_rows=[]
        for Y,(m,A) in zip(YTEST,geometry):
            alt={}
            for mode,stride,phase in [('cubic',2,0),('cubic',2,1),('pchip',1,0)]:
                am,aA=get(float(Y),mode,stride,phase)
                alt[f'{mode}_stride{stride}_phase{phase}']=dict(m=am,A=aA,A_relative_change=aA/A-1)
            local_fits=[]
            for nodes in [7,11,15]:
                ids=np.sort(np.argsort(abs(y-Y))[:nodes]);yy=y[ids]-Y
                scale=max(abs(yy))
                for degree in [3,4]:
                    coef=np.polynomial.polynomial.polyfit(yy/scale,original_A[ids],degree)
                    local_fits.append(dict(nodes=nodes,degree=degree,
                                           interval=[float(y[ids[0]]),float(y[ids[-1]])],
                                           A1=float(coef[1]/scale),A2=float(2*coef[2]/scale**2)))
            native_rows.append(dict(Y=float(Y),m=m,m_displacement=m-m0,A=A,A_over_A0=A/A0,
                                    fixed_X0_value=float(field(Y,0)[0,0]),
                                    inner_same_number_A_over_A0=get(float(f['ell']*Y))[1]/A0,
                                    alternatives=alt,kinematic=kinematic(field,float(Y)),
                                    local_envelope_polynomial_derivatives=local_fits))
        # Exact envelope derivatives of the saved bicubic interpolant.
        qy=float(field.ev(0,m0,dx=1,dy=0));qyy=float(field.ev(0,m0,dx=2,dy=0))
        qxy=float(field.ev(0,m0,dx=1,dy=1));qxx=float(field.ev(0,m0,dx=0,dy=2))
        wall=dict(A0=A0,P_algorithm=float(f['P']),relative_A0_minus_P=A0/float(f['P'])-1,
                  A1=qy,A2=qyy-qxy*qxy/qxx,qyy=qyy,ridge_motion_correction=-qxy*qxy/qxx,
                  A1_over_A0=qy/A0,A2_over_A0=(qyy-qxy*qxy/qxx)/A0)
        # Fit only original vertical nodes, with the wall value fixed.
        sel=y<=1.6;yf=y[sel];af=np.array([peak(x,row)[1] for row in q[sel]])
        fits=[]
        for high in [.05,.1,.2,.4,.8,1.6]:
            mask=yf<=high;yy=yf[mask];av=af[mask];scale=float(yy[-1])
            for deg in [2,3,4,5]:
                if len(yy)<deg+1:continue
                powers=np.arange(1,deg+1)
                mat=(yy[:,None]/scale)**powers
                coeff=np.linalg.lstsq(mat,av-A0,rcond=None)[0]
                a1=float(coeff[0]/scale);a2=float(2*coeff[1]/scale**2)
                fits.append(dict(window_requested=[0,high],actual_max=float(yy[-1]),nodes=len(yy),degree=deg,
                                 A1=a1,A2=a2,A1_relative_to_bicubic=a1/qy-1,
                                 A2_relative_to_bicubic=a2/wall['A2']-1,
                                 fit_linf=float(np.max(abs(mat@coeff-(av-A0))))))
        envelope_methods=[]
        for stride in [1,2,3]:
            for phase in range(stride):
                ids=np.unique(np.r_[0,np.arange(phase,len(yf),stride),len(yf)-1]).astype(int)
                for kind in ['cubic','pchip']:
                    cs=(CubicSpline(yf[ids],af[ids]) if kind=='cubic' else PchipInterpolator(yf[ids],af[ids]))
                    envelope_methods.append(dict(method=kind,stride=stride,phase=phase,
                                                 A1=float(cs(0,1)),A2_right=float(cs(0,2))))
        # Alternate source-field y interpolation (not just the envelope spline).
        source_methods=[]
        for stride in [2,3]:
            for phase in range(stride):
                ids=np.unique(np.r_[0,np.arange(phase,len(y),stride),len(y)-1]).astype(int)
                ff=RectBivariateSpline(y[ids],x,q[ids])
                sm,sa=peak(x,ff(0,x)[0]);d1=float(ff.ev(0,sm,dx=1,dy=0))
                d2=float(ff.ev(0,sm,dx=2,dy=0)-ff.ev(0,sm,dx=1,dy=1)**2/ff.ev(0,sm,dx=0,dy=2))
                source_methods.append(dict(stride=stride,phase=phase,A1=d1,A2=d2))
        vertical=[]
        yy=np.unique(np.r_[y[y<=16],16.])
        for X in [0.,1.,4.,16.]:
            col=field(yy,X)[:,0];cs=CubicSpline(yy,col,extrapolate=False)
            rr=cs.derivative().roots(extrapolate=False)
            points=np.r_[0.,rr[(rr>=0)&(rr<=16)],16.]
            yymax=float(points[np.argmax(cs(points))]);maximum=float(cs(yymax))
            vertical.append(dict(X=X,Y_of_max_on_0_16=yymax,vertical_max=maximum,
                                 horizontal_max_at_that_Y=get(yymax)[1]))
        output.append(dict(tau=float(f['tau']),ell=float(f['ell']),native_rows=native_rows,
                           coordinate_same_number_difference=relative(inner_A,grid_A),
                           wall_derivatives=wall,wall_polynomial_window_fits=fits,
                           envelope_spline_derivatives=envelope_methods,
                           source_y_subsampling_derivatives=source_methods,
                           fixed_X_vertical_maxima=vertical,
                           native_original_y_first_nodes=y[:14].tolist()))
    comparisons=[]
    for a,b in [(10,15),(15,20),(10,20)]:
        comparisons.append(dict(first=a,second=b,
                                same_native_Y=relative(native_curves[b],native_curves[a]),
                                same_inner_y=relative(inner_curves[b],inner_curves[a])))
    report=dict(definition='A(Y)=max_X rho_X(X,Y); m(Y)=argmax_X rho_X(X,Y), optionally subtract m(0)',
                scope='Saved native fields; no PDE evolution, no closure search.',
                comparison_grid=dict(coordinate_interval=[0,16],uniform_points=401,
                                     normalization='divide each frame by its own CubicSpline A(0)'),
                frames=output,time_comparisons=comparisons)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'amplitude_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'time_comparisons':comparisons,
                      'frames':[{'tau':r['tau'],'wall':r['wall_derivatives'],
                                 'coordinates':r['coordinate_same_number_difference'],
                                 'rows':[{k:v for k,v in n.items() if k!='alternatives'} for n in r['native_rows']]} for r in output]},indent=2))


if __name__=='__main__':main()
