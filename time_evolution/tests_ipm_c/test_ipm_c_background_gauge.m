function result=test_ipm_c_background_gauge()
%TEST_IPM_C_BACKGROUND_GAUGE Analytic rate offsets include actual c_omega.
% Synthetic complete-RHS test only: no Darcy/Poisson/native RHS/time step.
r=linspace(0,2,21)';theta=linspace(0,pi,11);
X=r*cos(theta);Y=r*sin(theta);a=.02;
R=X-a*X.^3+.04*Y.^2;Q=1-3*a*X.^2;
x=linspace(-1.9,1.9,81);n=numel(x);h=x(2)-x(1);
Dx=spdiags([-ones(n,1),zeros(n,1),ones(n,1)],[-1,0,1],n,n)/(2*h);
Dx(1,:)=0;Dx(1,1:3)=[-3,4,-1]/(2*h);
Dx(end,:)=0;Dx(end,end-2:end)=[1,-4,3]/(2*h);
settings=struct('wallIntegralRadius',1,'cOmegaGaugeFloor',1e-14, ...
    'pointStrainConditionFloor',1e-12,'referenceWallIntegral',17.25, ...
    'referenceWallGradientCentroid',.17,'peakPositionRestoringRate',.5,'pinX',0);
ops=struct('x',x,'Dx',Dx,'rescaling',settings);savedOps=ops;
wallR=x-a*x.^3;wallQ=1-3*a*x.^2;
stateOptions=struct('interpolation','linear','fieldOffset',R,'wallOffset',wallR);
seenAmplitude=[];
options=struct('rootTolerance',1e-10,'auditTolerance',1e-9, ...
    'measurementOptions',stateOptions,'rateMeasurementOffsets',@offsets);
callback=@(cl,delta)-((cl+.23)*X+(delta-.16)).*Q+.07*R;
[F,rates,info]=ipm_c.gauge.solve(R,r,ops,callback,.31,[.1,-.2],options);
assert(info.analyticRateOffsetsApplied&&info.fullAuditIncludesActualAmplitudeBackgroundRate&& ...
    info.fullMeasurement.rate.fieldOffsetUsed&&info.initialMeasurement.fieldOffsetUsed, ...
    'polar:testBackgroundOffsetUse','Both state and complete rate offsets must be used.');
assert(abs(rates.c_omega)>1e-5&&any(seenAmplitude==rates.c_omega)&&any(seenAmplitude==0), ...
    'polar:testBackgroundActualAmplitude','Base probes use zero amplitude and the full audit uses actual amplitude.');
assert(max(abs(info.auditResidual))<1e-9&&isequaln(ops,savedOps), ...
    'polar:testBackgroundFullAudit','The total RHS must satisfy the unchanged full I/P/C audit.');
expected=callback(rates.c_l,rates.delta)+rates.c_omega*R;
assert(isequal(F,expected),'polar:testBackgroundTotalF','Offsets must not alter the returned total F.');

% A deliberately incomplete background rate still leaves the total F
% unchanged, but maps c_omega*R through the coarse interpolation. The full
% audit must detect that error instead of silently omitting b_tau.
badOptions=options;badOptions.rateMeasurementOffsets=@missing_amplitude;
rejected=false;
try
    ipm_c.gauge.solve(R,r,ops,callback,.31,[.1,-.2],badOptions);
catch problem
    rejected=strcmp(problem.identifier,'polar:GaugeFullRhsAudit');
end
assert(rejected,'polar:testBackgroundMissingTimeRate', ...
    'An omitted amplitude background rate must fail the independent complete audit.');
result=struct('passed',true,'kind','synthetic_complete_background_rate_gauge_test', ...
    'nativeRhsCalls',0,'poissonCalls',0,'timeAdvanced',0,'rates',rates, ...
    'auditResidual',info.auditResidual,'missingAmplitudeRateRejected',rejected, ...
    'originalReferencesPreserved',isequaln(ops,savedOps));

    function [polarOffset,wallOffset]=offsets(cl,delta,cw)
        seenAmplitude(end+1)=cw;
        polarOffset=callback(cl,delta)+cw*R;
        wallOffset=-((cl+.23)*x+(delta-.16)).*wallQ+(.07+cw)*wallR;
    end
    function [polarOffset,wallOffset]=missing_amplitude(cl,delta,~)
        [polarOffset,wallOffset]=offsets(cl,delta,0);
    end
end
