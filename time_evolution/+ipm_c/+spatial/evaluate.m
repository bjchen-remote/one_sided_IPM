function [F,flow,info]=evaluate(rho,ops,scale)
%EVALUATE Centered density W=R-a and its origin amplitude a.
% Poisson, transport and the original I/P/C gauges see W; constants have
% zero spatial derivatives. The scalar a uses the same explicit RK stages.
g=ops.polarGrid;cfg=ops.polarCfg;cart=ops.cart;options=cfg.rhsOptions;
assert(isfield(scale,'densityOrigin')&&isfinite(scale.densityOrigin)&& ...
    all(rho(1,:)==0),'ipm_c:CenteredDensityState', ...
    'The centered density has a unique zero origin and a finite scalar amplitude.');
d=ipm_c.field.derivatives(rho,g,options);poissonOptions=options;
if isfield(cfg,'elliptic')
    names=fieldnames(cfg.elliptic);
    for k=1:numel(names),poissonOptions.(names{k})=cfg.elliptic.(names{k});end
end
[physical,pinfo]=ipm_c.spatial.poisson(d.rx,g,poissonOptions);
method='linear';if isfield(options,'measurementInterpolation'),method=options.measurementInterpolation;end
measurementOptions=struct('interpolation',method);
measurement=ipm_c.gauge.measure(rho,g.r,cart.originalWallOps,measurementOptions);
signedX=[-flipud(g.r(2:end));g.r];
wallU=[flipud(physical.u1(2:end,end));physical.u1(:,1)];
Uref=interp1(signedX,wallU,measurement.C,'linear');
guess=[cart.flow.c_l,cart.flow.c_r];
if isfield(ops.rescaling,'gaugeRateGuess'),guess=ops.rescaling.gaugeRateGuess;end
if isfield(options,'gaugeGuess'),guess=options.gaugeGuess;end
gaugeOptions=options;gaugeOptions.measurementOptions=measurementOptions;
transportOptions=options;if isfield(cfg,'transport'),transportOptions.transport=cfg.transport;end
% Prepare fixed-density WENO derivatives once for this RHS. The direct path
% is retained for independent operator equivalence and timing checks.
useTransportCache=~isfield(options,'disableTransportCache') || ~options.disableTransportCache;
if useTransportCache
    transportEvaluator=ipm_c.field.transport(rho,[],[],g,transportOptions);
end
try
    [F,rates,ginfo]=ipm_c.gauge.solve(rho,g.r,cart.originalWallOps,@base,Uref,guess,gaugeOptions);
catch cause
    if any(strcmp(cause.identifier,{'polar:GaugeRootUnresolved','polar:GaugeFullRhsAudit'}))
        failure=MException('ipm:UnresolvedWallIntegralGauge','Polar I/P/C solve failed: %s',cause.message);
        failure=addCause(failure,cause);throw(failure);
    end
    rethrow(cause);
end
transportOrigin=F(1,1);
F=F-transportOrigin;F(1,:)=0;
densityOriginRate=transportOrigin+rates.c_omega*scale.densityOrigin;
vx=(physical.u1-Uref)+rates.c_l*g.X+rates.delta;
vy=physical.u2+rates.c_l*g.Y;vy(:,[1,end])=0;
[cfl,position]=ipm_c.spatial.transportRate(vx,vy,g);
q=ipm_c.field.derivatives(F,g,options);
p=pinfo.linearSolveBackwardError;
assert(isfinite(p)&&p<=1e-10&&max(abs(ginfo.auditResidual))<=ginfo.auditTolerance, ...
    'ipm_c:SpatialAudit','The complete original gauge and B-spline linear solve must pass.');
flow=physical;flow.fullQ=q.rx;
flow.densityOriginRate=densityOriginRate;flow.transportOrigin=transportOrigin;
flow.c_l=rates.c_l;flow.c_x=rates.c_l;flow.c_y=rates.c_l;
flow.c_omega=rates.c_omega;flow.c_r=rates.c_r;
flow.canonicalCL=rates.c_l;flow.canonicalCX=rates.c_l;flow.canonicalCY=rates.c_l;
flow.canonicalCOmega=rates.c_omega;flow.canonicalCR=rates.c_r;flow.timeSpeed=1;
flow.conservativeSource=2*rates.c_l+rates.c_omega;
flow.canonicalConservativeSource=flow.conservativeSource;
flow.transportU1=vx;flow.transportU2=vy;flow.generatorCflRate=cfl;
flow.generatorReferenceVelocity=Uref;flow.generatorTranslationRate=rates.delta;
flow.generatorRootCondition=ginfo.rootCondition;flow.generatorRootTolerance=ginfo.rootTolerance;
flow.generatorRootResidual=max(abs(ginfo.rootResidual));
flow.gaugeResidual=ginfo.auditResidual;
flow.generatorIntegralResidual=ginfo.auditResidual(1);flow.generatorPeakResidual=ginfo.auditResidual(2);
flow.generatorPhaseResidual=ginfo.auditResidual(3);
flow.poissonResidual=p;flow.poissonSolveInfo=pinfo;
flow.poissonSolveInfo.relativeResidual=p;flow.poissonSolveInfo.algebraicResidual=p;
flow.rhoXInf=max(abs(d.rx),[],'all');flow.rhoYInf=max(abs(d.ry),[],'all');
flow.velocityInf=max(hypot(physical.u1,physical.u2),[],'all');flow.wallIntegral=measurement.I;
flow.globalWallPeak=measurement.P;flow.wallCentroid=measurement.C;
flow.wallIntegralRelativeDrift=measurement.relativeIntegralDrift;
flow.wallCentroidReference=cart.originalWallOps.rescaling.referenceWallGradientCentroid;
flow.wallCentroidRateResidual=ginfo.auditResidual(3);
info=struct('kind','centered_density_polar_spatial_backend_v1','poisson',pinfo,'gauge',ginfo, ...
    'maximumMappedPolarCflRate',cfl,'diagnosticTransportRateIncludingOrigin',cfl, ...
    'cflPosition',position,'originalWallFunctionalGridPreserved',true,'referencesReset',false, ...
    'physicalTimeReparameterized',false,'sameContinuousEquation',true,'sameSemiDiscreteEquation',false, ...
    'stateRepresentation','density_minus_origin','physicalRhsIndependentOfClocks',true);
flow.ipmCInfo=info;
    function B=base(cl,delta)
        bx=(physical.u1-Uref)+cl*g.X+delta;by=physical.u2+cl*g.Y;by(:,[1,end])=0;
        if useTransportCache
            B=transportEvaluator(bx,by);
        else
            B=ipm_c.field.transport(rho,bx,by,g,transportOptions);
        end
    end
end
