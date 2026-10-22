function result=test_ipm_c_gaussian_background()
%TEST_IPM_C_GAUSSIAN_BACKGROUND Derivatives, moving gauge, exact t0 observables.
p=struct('initialCondition','localized_gaussian_primitive','symmetryMode','half_plane', ...
    'gaussianAmplitude',.7,'gaussianCenterX',0,'gaussianWidthX',.75, ...
    'gaussianWidthY',1,'gaussianSkewWeight',.2);
X=[-1.3,-.2,.1,.8,1.7];Y=[0,.2,.7,1.1,.3];
scale=struct('logC_l',.37,'logC_omega',-.21,'X_shift',.43, ...
    'canonicalTime',2.1,'physicalTime',1.7);
before=scale;b=ipm_c_legacy_fixture('background',p,X,Y,scale);
xp=exp(-scale.logC_l)*(X-scale.X_shift);yp=exp(-scale.logC_l)*Y;
original=ipm.field.initialDensity(struct('X',xp,'Y',yp,'ny',1,'nx',numel(X)),p);
assert(max(abs(b.rho-exp(scale.logC_omega)*original))<1e-15 && isequaln(before,scale), ...
    'polar:testBackgroundDatum','Background must be the unchanged original primitive in current coordinates.');
h=1e-4;
xp1=ipm_c_legacy_fixture('background',p,X+h,Y,scale);xm1=ipm_c_legacy_fixture('background',p,X-h,Y,scale);
yp1=ipm_c_legacy_fixture('background',p,X,Y+h,scale);ym1=ipm_c_legacy_fixture('background',p,X,Y-h,scale);
pp=ipm_c_legacy_fixture('background',p,X+h,Y+h,scale);pm=ipm_c_legacy_fixture('background',p,X+h,Y-h,scale);
mp=ipm_c_legacy_fixture('background',p,X-h,Y+h,scale);mm=ipm_c_legacy_fixture('background',p,X-h,Y-h,scale);
errors=struct('rx',max(abs((xp1.rho-xm1.rho)/(2*h)-b.rx)), ...
    'ry',max(abs((yp1.rho-ym1.rho)/(2*h)-b.ry)), ...
    'rxx',max(abs((xp1.rho-2*b.rho+xm1.rho)/h^2-b.rxx)), ...
    'rxy',max(abs((pp.rho-pm.rho-mp.rho+mm.rho)/(4*h^2)-b.rxy)), ...
    'ryy',max(abs((yp1.rho-2*b.rho+ym1.rho)/h^2-b.ryy)));
assert(max(structfun(@(v)v,errors))<2e-7, ...
    'polar:testBackgroundDerivatives','Analytic background gradient/Hessian failed independent finite differences.');
cl=.28;cr=.17;cw=-.14;Sdot=cl*scale.X_shift+cr;
sp=scale;sm=scale;
sp.logC_l=scale.logC_l+h*cl;sm.logC_l=scale.logC_l-h*cl;
sp.logC_omega=scale.logC_omega+h*cw;sm.logC_omega=scale.logC_omega-h*cw;
sp.X_shift=scale.X_shift+h*Sdot;sm.X_shift=scale.X_shift-h*Sdot;
bp=ipm_c_legacy_fixture('background',p,X,Y,sp);bm=ipm_c_legacy_fixture('background',p,X,Y,sm);
tau=-(cl*X+cr).*b.rx-cl*Y.*b.ry+cw*b.rho;
tauError=max(abs((bp.rho-bm.rho)/(2*h)-tau));
assert(tauError<2e-8,'polar:testBackgroundTau', ...
    'Background time derivative must include Sdot=cl*S+cr at nonzero scale/shift.');

% A local deterministic analytic t0 fixture replaces experimental MAT files.
% Reference values come from the unchanged Cartesian wall functionals before
% any polar transfer. This checks representation and offsets without reading
% or modifying a stored trajectory.
[g,cart,initial,shared]=ipm_c_test_fixture(65);
bg=ipm_c_legacy_fixture('background',p,g.X,g.Y,cart.scale);
bw=ipm_c_legacy_fixture('background',p,cart.originalWallOps.x,zeros(size(cart.originalWallOps.x)),cart.scale);
assert(isequaln(bg.rho,initial) && isequaln(bw.rho,shared.wallRho), ...
    'polar:testBackgroundBitwiseDatum','t0 background must reproduce direct original samples exactly.');
opts=struct('interpolation','spline','fieldOffset',bg.rho,'wallOffset',bw.rho);
m=ipm_c.gauge.measure(initial,g.r,cart.originalWallOps,opts);
assert(isequaln([m.I,m.P,m.C],shared.IPC) && ...
    isequaln(m.q,shared.wallQ),'polar:testBackgroundBitwiseReferences', ...
    'At v=0 offset measurement must reproduce original I/P/C and wall q bit for bit.');
opts.rateField=.3*bg.rho;opts.rateFieldOffset=.3*bg.rho;opts.wallRateOffset=.3*bw.rho;
mr=ipm_c.gauge.measure(initial,g.r,cart.originalWallOps,opts);
expectedRate=(.3*bw.rho)*cart.originalWallOps.Dx';
assert(isequaln(mr.rate.qDot,expectedRate), ...
    'polar:testBackgroundRateOffset','Complete RHS offsets must use the same original wall derivative.');
alias=opts;alias=rmfield(alias,'wallRateOffset');alias.rateWallOffset=.3*bw.rho;
ma=ipm_c.gauge.measure(initial,g.r,cart.originalWallOps,alias);
assert(isequaln(ma.rate,mr.rate),'polar:testBackgroundRateAlias','Rate-offset spellings must agree.');
result=struct('passed',true,'derivativeErrors',errors,'tauDerivativeError',tauError, ...
    'fixture','deterministic_analytic_t0_without_external_frame','t0ReferencesBitwiseEqual',true,'originalIPC',[m.I,m.P,m.C], ...
    'referenceReset',false,'nativePdeCalls',0,'poissonCalls',0);
end
