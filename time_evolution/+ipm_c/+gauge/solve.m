function [fullF,rates,info] = solve(R,r,wallOps,baseEvaluator,Uref,originalGuess,options)
%POLAR_SOLVE_GAUGE Solve the original I/P/C semantic gauge on a polar RHS.
% baseEvaluator(c_l,delta) returns the COMPLETE polar transport base RHS B,
% with c_omega=0 and centered translation delta=c_r+Uref. Uref is frozen
% once per state from the complete Darcy field at the original wall C.
% originalGuess is [c_l,c_r] in physical gauge coordinates, or a scalar
% nominal c_l. No persisted reference is regenerated or projected here.
%
% [fullF,rates,info] returns B+c_omega*R, physical c_r=delta-Uref, and an
% independent complete-RHS I/P/C audit. options.measurementOptions controls
% ipm_c.gauge.measure (default data-linear interpolation). Explicit pchip
% is allowed diagnostically, but its nonlinear mapping can fail the audit.
% The continuum gauge formulas are retained; the spatial WENO operator is
% that supplied by baseEvaluator, not asserted identical to Cartesian WENO.
% Optional rateMeasurementOffsets(cl,delta,c_omega) returns the analytic
% polar/wall offsets of the COMPLETE rate field. The background time rate
% must therefore be included with c_omega=0 for base probes and with the
% actual c_omega for the independent complete-RHS audit.

if nargin<7,options=struct();end
assert(isstruct(options) && isscalar(options), ...
    'polar:GaugeOptions','Gauge options must be a scalar structure.');
assert(isa(baseEvaluator,'function_handle'), ...
    'polar:GaugeCallback','baseEvaluator must be a function handle.');
validateattributes(Uref,{'numeric'},{'real','finite','scalar'});
validateattributes(originalGuess,{'numeric'},{'real','finite','vector','nonempty'});
assert(ismember(numel(originalGuess),[1,2]), ...
    'polar:GaugeGuess','originalGuess must be nominal c_l or [c_l,c_r].');
measurementOptions=struct();
if isfield(options,'measurementOptions'),measurementOptions=options.measurementOptions;end
if isfield(measurementOptions,'rateField')
    measurementOptions=rmfield(measurementOptions,'rateField');
end
rateOffsets=[];
if isfield(options,'rateMeasurementOffsets')
    rateOffsets=options.rateMeasurementOffsets;
    assert(isa(rateOffsets,'function_handle'),'polar:GaugeRateOffsets', ...
        'rateMeasurementOffsets must return analytic polar and original-wall rate offsets.');
end
rootTolerance=1e-6;auditTolerance=1e-5;wallTolerance=1e-10;
maximumEvaluations=160;maximumIterations=14;
names={'rootTolerance','auditTolerance','wallTolerance'};
limits=[rootTolerance,auditTolerance,wallTolerance];
for k=1:numel(names)
    if isfield(options,names{k})
        value=options.(names{k});
        assert(isnumeric(value) && isscalar(value) && isfinite(value) && ...
            value>0 && value<=limits(k),'polar:GaugeTolerance', ...
            'Gauge gates may be tightened but not relaxed.');
        limits(k)=value;
    end
end
rootTolerance=limits(1);auditTolerance=limits(2);wallTolerance=limits(3);
[m,preparedRate]=ipm_c.gauge.measure(R,r,wallOps,measurementOptions);
kappa=wallOps.rescaling.peakPositionRestoringRate;
assert(isnumeric(kappa) && isscalar(kappa) && isfinite(kappa) && kappa>=0, ...
    'polar:GaugeRestoringRate','The original centroid restoring rate must be finite and nonnegative.');
phase=kappa*(m.C-wallOps.rescaling.referenceWallGradientCentroid);
target=[0;0;phase];evaluations=0;iterations=0;polished=false;
originalGuess=originalGuess(:)';nominalL=originalGuess(1);
% Use exactly the original wall measurement operators for the cold seed.
[~,physicalRate]=evaluate([0;0]);
continuous=[observe_wall(-m.x.*m.q),observe_wall(-m.q),observe_wall(m.wallValues)];
cold=lsqminnorm(continuous,-physicalRate-target);
seeds=[cold(1:2)';nominalL,-nominalL*m.C-phase;1,0;0,0;-1,0];
warmGuess=[];
if numel(originalGuess)==2
    warmGuess=originalGuess;
elseif isfield(wallOps.rescaling,'gaugeRateGuess') && ...
        isnumeric(wallOps.rescaling.gaugeRateGuess) && ...
        numel(wallOps.rescaling.gaugeRateGuess)==2 && ...
        all(isfinite(wallOps.rescaling.gaugeRateGuess))
    warmGuess=reshape(wallOps.rescaling.gaugeRateGuess,1,2);
end
if ~isempty(warmGuess)
    centered=[warmGuess(1),warmGuess(2)+Uref];
    if max(abs(centered))<=8*max(1,max(abs(cold(1:2))))
        seeds=[centered;seeds];
    else
        seeds=[seeds;centered];
    end
end
found=false;condition=NaN;last=[NaN;NaN];lastX=[NaN;NaN];
for seedIndex=1:size(seeds,1)
    z=seeds(seedIndex,:)';
    for iteration=1:maximumIterations
        iterations=iterations+1;
        if evaluations>=maximumEvaluations,break;end
        [g,v,base,baseMeasurement]=evaluate(z);
        merit=max(abs(g));last=g;lastX=z;
        if merit<=rootTolerance
            if merit>1e-10 && ~polished
                J=jacobian(z);trial=z-lsqminnorm(J,g);
                [trialG,trialV,trialBase,trialMeasurement]=evaluate(trial);
                polished=true;
                if max(abs(trialG))<merit
                    z=trial;g=trialG;v=trialV;base=trialBase;
                    baseMeasurement=trialMeasurement;
                end
            end
            found=true;break;
        end
        J=jacobian(z);condition=rcond(J);step=-lsqminnorm(J,g);
        stepLimit=4*max(1,norm(z,inf));
        if norm(step,inf)>stepLimit,step=step*(stepLimit/norm(step,inf));end
        accepted=false;
        for damping=2.^(-(0:7))
            trial=z+damping*step;[trialG,~]=evaluate(trial);
            if max(abs(trialG))<merit*(1-1e-3*damping) || ...
                    max(abs(trialG))<=rootTolerance
                z=trial;accepted=true;break;
            end
        end
        if ~accepted
            for lambda=max(1,norm(J,'fro')^2)*10.^(-4:1)
                trial=z-(J'*J+lambda*eye(2))\(J'*g);
                [trialG,~]=evaluate(trial);
                if max(abs(trialG))<merit,z=trial;accepted=true;break;end
            end
        end
        if ~accepted,break;end
    end
    if found || evaluations>=maximumEvaluations,break;end
end
if ~found
    error('polar:GaugeRootUnresolved', ...
        'Polar I/P/C root unresolved after %d evaluations: cl=%g delta=%g G=[%g %g].', ...
        evaluations,lastX(1),lastX(2),last(1),last(2));
end
if ~isfinite(condition),condition=rcond(jacobian(z));end
cOmega=-v(2);cL=z(1);delta=z(2);cR=delta-Uref;
% Reinvoke the full callback, then sample the actual full polar RHS. This
% catches callback inconsistencies and the nonlinearity of explicit pchip.
auditBase=checked_base(cL,delta);evaluations=evaluations+1;
auditBaseMeasurement=rate_measurement(auditBase,cL,delta,0,true);
wallDifference=max(abs(auditBaseMeasurement.rate.wallValues- ...
    baseMeasurement.rate.wallValues))/max(1,max(abs(auditBaseMeasurement.rate.wallValues)));
fullF=auditBase+cOmega*R;
audit=rate_measurement(fullF,cL,delta,cOmega,true);residual=audit.rate.normalized+target;
if wallDifference>wallTolerance || any(~isfinite(residual)) || ...
        max(abs(residual))>auditTolerance
    error('polar:GaugeFullRhsAudit', ...
        ['Polar complete RHS audit failed: interpolation=%s wall=%g, ' ...
        'I/P/C=[%g %g %g]. PCHIP is not data-linear.'], ...
        m.interpolation,wallDifference,residual(1),residual(2),residual(3));
end
rates=struct('c_l',cL,'delta',delta,'c_r',cR,'c_omega',cOmega,'Uref',Uref);
info=struct('definition','original_IPC_formulas_on_polar_complete_RHS', ...
    'initialMeasurement',m,'fullMeasurement',audit,'rootResidual',g, ...
    'auditResidual',residual,'rootTolerance',rootTolerance, ...
    'auditTolerance',auditTolerance,'wallRelativeDifference',wallDifference, ...
    'rootCondition',condition,'iterations',iterations,'evaluations',evaluations, ...
    'seedIndex',seedIndex,'warmStartUsed',~isempty(warmGuess), ...
    'originalGuess',originalGuess,'Uref',Uref,'target',target, ...
    'referenceWallIntegral',wallOps.rescaling.referenceWallIntegral, ...
    'referenceWallGradientCentroid',wallOps.rescaling.referenceWallGradientCentroid, ...
    'interpolation',m.interpolation,'dataLinearInterpolation',m.dataLinearInterpolation, ...
    'netRootPolished',polished,'baseRhs',auditBase, ...
    'analyticRateOffsetsApplied',~isempty(rateOffsets), ...
    'fullAuditIncludesActualAmplitudeBackgroundRate',~isempty(rateOffsets));

    function values=observe_wall(wallValues)
        dq=ipm_c.gauge.wallDerivative(wallValues,wallOps.Dx);
        values=[m.integral.weights*dq'/m.I; ...
            m.peak.weights*dq(m.peak.indices)'/m.P; ...
            m.centroid.weights*dq(m.centroid.indices)'];
    end
    function value=rate_measurement(rhs,cl,translation,cw,independentAudit)
        if nargin<5,independentAudit=false;end
        rateOptions=measurementOptions;rateOptions.rateField=rhs;
        if ~isempty(rateOffsets)
            [polarOffset,wallOffset]=rateOffsets(cl,translation,cw);
            assert(isnumeric(polarOffset)&&isreal(polarOffset)&& ...
                isequal(size(polarOffset),size(R))&&all(isfinite(polarOffset),'all')&& ...
                isnumeric(wallOffset)&&isreal(wallOffset)&&isvector(wallOffset)&& ...
                numel(wallOffset)==numel(wallOps.x)&&all(isfinite(wallOffset),'all'), ...
                'polar:GaugeRateOffsets','Analytic rate offsets must match full polar/original-wall fields.');
            rateOptions.rateFieldOffset=polarOffset;
            rateOptions.wallRateOffset=reshape(wallOffset,1,[]);
        end
        if independentAudit
            % Both callback-repeat and full-RHS audits rebuild observations
            % from the actual density, independently of the probe closure.
            value=ipm_c.gauge.measure(R,r,wallOps,rateOptions);
        else
            value=preparedRate(rhs,rateOptions);
        end
    end
    function rhs=checked_base(cl,translation)
        rhs=baseEvaluator(cl,translation);
        assert(isnumeric(rhs) && isreal(rhs) && isequal(size(rhs),size(R)) && ...
            all(isfinite(rhs),'all'),'polar:GaugeCallbackSize', ...
            'baseEvaluator must return one finite complete polar base RHS.');
    end
    function [gauge,values,rhs,measured]=evaluate(z)
        evaluations=evaluations+1;rhs=checked_base(z(1),z(2));
        measured=rate_measurement(rhs,z(1),z(2),0);values=measured.rate.normalized;
        gauge=[values(1)-values(2);values(3)+phase];
        assert(all(isfinite(gauge)),'polar:GaugeFinite','Gauge rates must be finite.');
    end
    function J=jacobian(z)
        J=zeros(2,2);
        for column=1:2
            h=1e-4*max(1,abs(z(column)));e=zeros(2,1);e(column)=h;
            plus=evaluate(z+e);minus=evaluate(z-e);J(:,column)=(plus-minus)/(2*h);
        end
    end
end
