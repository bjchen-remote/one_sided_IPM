"""One bounded shared-parameter fit of the user-specified level-set family.

f is the observed wall inverse, never extrapolated or modified. This is a
geometry-only finite-window experiment, not a steady closure or rho fit.
"""
import os,json
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline
from scipy.optimize import minimize
from streamfunction_sources import load_native,psi_over_y

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_levelset_20260929'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    p=load_native(20);x=p['x'];yn=p['y'];H=psi_over_y(p)
    sigma=float(p['rates'][0]-p['rates'][1]);v=float(p['rates'][2])
    keep=abs(x)<=8192;xx=x[keep];Hw=H[0,keep]
    if not np.all(np.diff(Hw)>0):raise ValueError('The requested wall inverse is not monotone on the declared window')
    f=PchipInterpolator(Hw,xx,extrapolate=False)
    h=np.unique(np.r_[np.linspace(-24,-2,45),np.linspace(-2,2,81),np.linspace(2,8,25)])
    yy=np.array([.0625,.125,.25,.5,1,2,4,8,16,32,64,128,256,512,1024.])
    train=np.arange(len(yy))%2==0;Y=1024.
    hrows=PchipInterpolator(yn,H,axis=0,extrapolate=False)(yy)[:,keep]
    observations=[];row_support=[]
    for y,row in zip(yy,hrows):
        if not np.all(np.diff(row)>0):raise ValueError(f'Nonmonotone source H at y={y}')
        if h[0]<row[0] or h[-1]>row[-1]:raise ValueError('Selected h levels are not observed at every height')
        inverse=PchipInterpolator(row,xx,extrapolate=False)
        observations.append(inverse(h));row_support.append([float(row[0]),float(row[-1])])
    observations=np.array(observations)
    HH,YY=np.meshgrid(h,yy);scale=1+np.abs(f(HH))+YY
    # Bounds guarantee f support for every h in [-24,8], y in [0,1024].
    shiftlo=-.48*(Hw[-1]-h[-1]);shifthi=.48*(h[0]-Hw[0])
    hscale=8.;wscale=v*hscale/sigma**3
    bounds=[(shiftlo,shifthi),(shiftlo,shifthi),(-10*hscale**2/wscale,10*hscale**2/wscale),(-5*hscale,5*hscale),(0,100*hscale**2)]
    def unpack(z):
        alpha,beta=z[0]/Y,z[1]/Y**2;a2,b1,b2=z[2]*wscale/hscale**2,z[3]/hscale,z[4]/hscale**2
        a0=(alpha+.5)/sigma
        a1=a0*b1+(2*alpha+1)*v/sigma**3
        return alpha,beta,a0,a1,a2,b1,b2
    def denominator_minimum(z):
        *_,b1,b2=unpack(z);locations=[h[0],h[-1]]
        if b2>0 and h[0]<=-b1/(2*b2)<=h[-1]:locations.append(-b1/(2*b2))
        return min(1+b1*t+b2*t*t for t in locations)
    def model(z,hg,yg,derivative=False):
        alpha,beta,a0,a1,a2,b1,b2=unpack(z)
        shifted=hg-alpha*yg-beta*yg**2
        if shifted.min()<Hw[0] or shifted.max()>Hw[-1]:raise ValueError('No extrapolation of f is allowed')
        den=1+b1*hg+b2*hg**2;num=a0+a1*hg+a2*hg**2
        w=num/den
        wp=((a1+2*a2*hg)*den-num*(b1+2*b2*hg))/den**2
        return (f(shifted,1)+yg*wp) if derivative else (f(shifted)+yg*w)
    def raw_objective(z):
        residual=(model(z,HH,YY)-observations)/scale
        return float(np.mean(residual[train]**2))
    ch=np.unique(np.r_[np.linspace(-24,8,97),np.linspace(-2,2,65)])
    cy=np.unique(np.r_[0,yy,np.linspace(0,Y,17)])
    CH,CY=np.meshgrid(ch,cy)
    def monotone(z):
        if denominator_minimum(z)<=0:return np.full(CH.size,-100.-abs(denominator_minimum(z)))
        return (model(z,CH,CY,True)/f(CH,1)-1e-6).ravel()
    ah=np.unique(np.r_[np.linspace(-24,8,513),np.linspace(-2,2,257)])
    ay=np.unique(np.r_[0,yy,np.geomspace(.001,Y,81),np.linspace(0,Y,65)])
    AH,AY=np.meshgrid(ah,ay)
    initial=np.zeros(5);initial_objective=raw_objective(initial)
    initial_diagnostics=dict(parameters=dict(zip(('alpha','beta','a0','a1','a2','b1','b2'),map(float,unpack(initial)))),
        weightedTrainRms=float(np.sqrt(initial_objective)),minimumMonotonicityConstraint=float(monotone(initial).min()),
        note='This is a feasible initial model, not a fitted solution; it shows that solver failure does not establish emptiness of the constrained family.')
    best=dict(z=initial.copy(),objective=initial_objective,updates=0)
    def objective(z):
        dmin=denominator_minimum(z)
        if dmin<=1e-5:return 1e3+1e3*(dmin-1e-5)**2
        loss=raw_objective(z)
        if loss<best['objective'] and monotone(z).min()>0 and model(z,AH,AY,True).min()>0:
            best.update(z=z.copy(),objective=loss,updates=best['updates']+1)
        return loss/initial_objective
    fit=minimize(objective,initial,method='SLSQP',bounds=bounds,
        constraints=[dict(type='ineq',fun=monotone),dict(type='ineq',fun=lambda z:denominator_minimum(z)-1e-5)],
        options=dict(maxiter=150,ftol=1e-10,disp=False))
    selected=best['z'];params=unpack(selected);prediction=model(selected,HH,YY)
    xh=model(selected,AH,AY,True)
    alpha,beta,a0,a1,a2,b1,b2=params
    argument=AH-alpha*AY-beta*AY**2;den=1+b1*ah+b2*ah**2
    # Evaluate the quadratic denominator's actual minimum on the h interval.
    candidates=[h[0],h[-1]]
    if b2>0 and h[0]<=-b1/(2*b2)<=h[-1]:candidates.append(-b1/(2*b2))
    dmin=min(1+b1*t+b2*t*t for t in candidates)
    fp,fpp=float(f(0,1)),float(f(0,2));wp0=a1-a0*b1
    q0=2*(-alpha*fp+a0)/fp
    qx0=2*(wp0*fp-a0*fpp)/fp**3
    wallcs=CubicSpline(x,H[0]);hc1=float(wallcs(0,1));hc2=float(wallcs(0,2))
    jetchecks=[]
    for radius,degree in [(.1,3),(.2,4),(.4,4)]:
        m=abs(x)<radius;co=np.polynomial.polynomial.polyfit(x[m],H[0,m],degree)
        jetchecks.append(dict(xRadius=radius,degree=degree,points=int(m.sum()),
            fprime=1/float(co[1]),fsecond=-2*float(co[2])/float(co[1])**3))
    heights=[]
    for i,y in enumerate(yy):
        err=prediction[i]-observations[i];normalized=err/scale[i]
        heights.append(dict(y=float(y),role='train' if train[i] else 'held_out',
            sourceHSupport=row_support[i],xRmsError=float(np.sqrt(np.mean(err**2))),
            maximumAbsoluteXError=float(abs(err).max()),weightedRms=float(np.sqrt(np.mean(normalized**2)))))
    def summary(mask):
        err=prediction[mask]-observations[mask]
        return dict(weightedRms=float(np.sqrt(np.mean((err/scale[mask])**2))),xRmsError=float(np.sqrt(np.mean(err**2))),maximumAbsoluteXError=float(abs(err).max()))
    valid=bool(xh.min()>0 and dmin>0 and argument.min()>=Hw[0] and argument.max()<=Hw[-1])
    result=dict(family='X(h,y)=f(h-alpha*y-beta*y^2)+y*(a0+a1*h+a2*h^2)/(1+b1*h+b2*h^2)',
        source=p['source'],coordinates='inner x,y,H; H=psi/y, wall H=-saved centered u1',
        domain=dict(h=[float(h[0]),float(h[-1])],y=[0,Y],sourceXWindow=[float(xx[0]),float(xx[-1])],
            wallInverseHSupport=[float(Hw[0]),float(Hw[-1])],wallXAtHBounds=[float(f(h[0])),float(f(h[-1]))],
            scope='A declared finite observation window, not the full saved domain or the whole half-plane. Shared parameters only within this window.',
            noExtrapolation=True,sourceAllSelectedRowsStrictlyMonotone=True),
        parameters=dict(zip(('alpha','beta','a0','a1','a2','b1','b2'),map(float,params))),
        constraints=dict(core='a0=(alpha+1/2)/sigma; a1=a0*b1+(2alpha+1)*v/sigma^3',
            denominator='b2>=0; exact quadratic minimum on the declared h interval constrained positive. No claim outside that interval.',
            variableScaling=dict(Ymax=Y,Hscale=hscale,Wscale=wscale,objectiveScale=initial_objective),
            optimizerBounds=dict(zip(('alpha*Ymax','beta*Ymax^2','a2*Hscale^2/Wscale','b1*Hscale','b2*Hscale^2'),bounds)),
            monotonicity='X_h/fprime(h)>=1e-6 on a predeclared collocation mesh; checked again on a finer independent mesh. Mesh audit is not a proof between nodes.',
            support='Bounds on alpha*Ymax and beta*Ymax^2 guarantee shifted h remains in the observed f support over the entire declared rectangle.'),
        optimizer=dict(method='one scaled SLSQP run from the exact linear feasible baseline; earlier unscaled failure retained separately',
            success=bool(fit.success),message=str(fit.message),iterations=int(fit.nit),selectedObjective=float(best['objective']),
            terminalScaledObjective=float(fit.fun),validFiniteWindowFit=valid,acceptedFeasibleImprovements=best['updates'],
            selection='Best objective among candidates passing the independent dense support/monotonicity/denominator gate, including the initial baseline. Not a claimed optimum.'),
        feasibleInitialPoint=initial_diagnostics,
        audit=dict(minimumXh=float(xh.min()),minimumDenominator=float(dmin),shiftedHRange=[float(argument.min()),float(argument.max())],
            minimumConstraintAtSelectedCandidate=float(monotone(selected).min()),
            predictedXRange=[float(prediction.min()),float(prediction.max())],
            fractionPredictionsOutsideSourceXWindow=float(np.mean((prediction<xx[0])|(prediction>xx[-1])))),
        core=dict(sigma=sigma,v=v,targetF0=0,targetFprime=1/sigma,targetFsecond=2*v/sigma**3,
            pchipF0=float(f(0)),pchipFprime=fp,pchipFsecond=fpp,
            sourceCubicInverseFprime=1/hc1,sourceCubicInverseFsecond=-hc2/hc1**3,localWallPolynomialChecks=jetchecks,
            formalCoefficientConstraintResidual0=float(a0-(alpha+.5)/sigma),
            formalCoefficientConstraintResidual1=float(wp0-(2*alpha+1)*v/sigma**3),
            actualPchipModelXy0=float(-alpha*fp+a0),actualPchipModelQ0=q0,actualPchipModelQx0=qx0,
            note='The saved wall is not replaced by exact steady jets. Coefficient constraints use the target equations, while actual interpolated-f core checks retain the observed discrepancy. PCHIP second derivative is interpolation-sensitive.'),
        train=summary(train),heldOut=summary(~train),heights=heights,
        limitations=['This is a single bounded local optimization, not a proof of global optimality.',
            'Finite-window support was chosen before fitting; no adaptive shrinking is used to rescue an invalid fit.',
            'The PCHIP inverse is C1 and has piecewise second derivatives; no full Poisson/transport residual or steady closure is claimed.',
            'rho was neither fitted nor altered.'])
    (OUT/'global_fit.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(OUT/'global_fit.npz',h=h,y=yy,trainingMask=train,sourceX=observations,predictedX=prediction,
        parameters=np.array(params),parameterNames=np.array(['alpha','beta','a0','a1','a2','b1','b2']),
        wallH=Hw,wallX=xx,wallInverseKnots=f.x,wallInverseCoefficients=f.c,
        auditH=ah,auditY=ay,auditXh=xh)
    report(result)
    print(json.dumps({k:result[k] for k in ('parameters','optimizer','audit','core','train','heldOut')},indent=2))

def report(j):
    lines=['# 指定 H 等值线族的共享参数辅助检验','',
        '**只报告通过稠密单调性、插值支持和分母门控的最佳已观测候选，不把优化终止点自动当作结果。** 这是有限预算内的辅助检验，不是最优性或函数族不可能性证明。', '',
        '**范围：有限观测窗，不是全半平面 global。** 固定 f 是保存墙面 H=−u1 的单调反函数，未修改其原点 jets，未外推。早先未缩放尝试的失败单独保留；本结果来自用户指定的精确可行基线及一次尺度归一化有界 SLSQP。a0、a1 受精确核心关系约束，实际剩五个自由参数。','',
        '`X=f(h−αy−βy²)+y(a0+a1h+a2h²)/(1+b1h+b2h²)`。','',
        f"定义域：h∈{j['domain']['h']}、y∈{j['domain']['y']}；源 x 窗={j['domain']['sourceXWindow']}；固定 f 的 h 支持={j['domain']['wallInverseHSupport']}。墙面 h 两端对应 x={j['domain']['wallXAtHBounds']}。每个选定高度的源 H 都在该 x 窗严格递增。",'',
        f"优化返回：{j['optimizer']['message']}；已选候选的几何有效性标志={j['optimizer']['validFiniteWindowFit']}。此标志要求细网格 Xh>0、分母无零点及 f 支持有效；不表示优化收敛。保留的可行改善共 {j['optimizer']['acceptedFeasibleImprovements']} 次。",'',
        f"初始点满足单调性约束（最小归一化余量 {j['feasibleInitialPoint']['minimumMonotonicityConstraint']:.7g}），但训练加权 RMS 为 {j['feasibleInitialPoint']['weightedTrainRms']:.7g}，并非合理拟合。它说明约束族不为空；后续局部优化的失败不能解释为全族被排除。",'',
        '| 参数 | 值 |','|---|---:|']
    for k,v in j['parameters'].items():lines.append(f'| {k} | {v:.12g} |')
    lines+=['',f"核查：min Xh={j['audit']['minimumXh']:.7g}；分母解析最小值={j['audit']['minimumDenominator']:.7g}；移位 h 范围={j['audit']['shiftedHRange']}。单调性细网格核查不能代替连续域证明。",'',
        f"**可定义且单调不等于形状合适。** 本候选训练加权 RMS={j['train']['weightedRms']:.7g}、留出={j['heldOut']['weightedRms']:.7g}，仍有很大偏离；预测 X 范围={j['audit']['predictedXRange']}，部分预测甚至超出本次可靠源 x 窗。没有将此可行候选称为成功的整体几何骨架。",'',
        '目标函数为各训练点的 `(Xmodel−Xsource)/(1+|f(h)|+y)` 的均方；多高度共享完全相同的系数。留出高度从未进入目标函数。','',
        '| 高度 y | 用途 | 加权 RMS | X RMS | 最大绝对 X 误差 |','|---:|---|---:|---:|---:|']
    for r in j['heights']:lines.append(f"| {r['y']:g} | {r['role']} | {r['weightedRms']:.6g} | {r['xRmsError']:.6g} | {r['maximumAbsoluteXError']:.6g} |")
    c=j['core']
    lines+=['','核心关系没有通过修改 f 强制成立：',
        f"目标 f′(0)={c['targetFprime']:.9g}、f″(0)={c['targetFsecond']:.9g}。固定 PCHIP 给 {c['pchipFprime']:.9g}、{c['pchipFsecond']:.9g}；直接墙面 H 的三次插值反函数 jets 给 {c['sourceCubicInverseFprime']:.9g}、{c['sourceCubicInverseFsecond']:.9g}。局部多项式诊断另存 JSON，表明二阶 jet 对插值处理更敏感。",'',
        f"a0、a1 按目标 σ、v 关系解析固定，但用实际 PCHIP f 反算得到 q(0)={c['actualPchipModelQ0']:.9g}、qx(0)={c['actualPchipModelQx0']:.9g}，应与规范 1、0 区分。",'',
        '这是几何辅助检验。固定 PCHIP 仅 C1，高阶插值导数不足以据此宣称完整 Poisson／输运闭合；没有拟合或改动 ρ。若有效性标志失败，本次有限窗拟合不成立，未缩域重试，也未增加参数或更换函数族。', '',
        '数据：global_fit.json、global_fit.npz。复现脚本：solver/steady_analysis/streamfunction_levelset_global_fit.py。']
    (OUT/'GLOBAL_FIT.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
