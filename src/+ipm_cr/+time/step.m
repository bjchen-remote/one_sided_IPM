function [out,info]=step(state,ctx,finalTau,controls)
%STEP One explicit step; numerical verification belongs to the test suite.
% CFL selects the initial interval. RK4 partitions ordinary continuous peak
% crossings of original spline knots by four secant iterates of S_x=0.
% Numerical diagnostics never accept/reject a state or select a fallback.
if nargin<4,controls=struct();end
policy=struct('cfl',ctx.opts.cfl,'maxDt',ctx.opts.maxDt, ...
    'timeIntegrator',ctx.opts.timeIntegrator);
for name=fieldnames(controls)',policy.(name{1})=controls.(name{1});end
method=policy.timeIntegrator;
h=min([policy.maxDt,.95*policy.cfl/state.rhs.rate, ...
    finalTau-state.scale.canonicalTime]);
initialH=h;probeH=h;hasResidual=isfield(state,'eta');calls=0;peaks=struct([]);
callbacks=struct('stage',@stage,'finish',@finish,'notePeak',@note_peak);
started=tic;
[out,records]=trial_at(h);
partition=struct('kind','direct_tableau_interval','fitted',false, ...
    'initialDt',h,'segmentDt',h,'knotX',NaN,'knotIndex',NaN,'direction',0, ...
    'secantIterations',0,'extraRhsEvaluations',0,'knotSlope',NaN,'rawPeakKnotDefect',NaN, ...
    'coordinateIdentityUlps',256,'startClassification',struct(),'endClassification',struct());
if strcmp(method,'rk4')
    [out,records,h,partition]=fit_knot(state,out,records,h,ctx.wall.x,@trial_at,partition);
end
partition.extraRhsEvaluations=calls-4;
if strcmp(method,'ssprk3'),partition.extraRhsEvaluations=calls-3;end
[maximumRate,index]=max(records.rates);
[~,clock]=ipm_cr.time.stageScale(state.scale,pack(out.scale),records.physicalClockIncrement);
order=4;stages=4;scheme='rk4_fixed4_spline_knot_partition_v4';
if strcmp(method,'ssprk3'),order=3;stages=3;scheme='ssprk3_direct_hybrid_v3';end
info=struct('method',scheme,'timeIntegrator',method, ...
    'temporalOrder',order,'spatialOrderCertified',false,'dt',h,'seconds',toc(started), ...
    'initialProposedDt',initialH,'timePartition',partition, ...
    'rate',state.rhs.rate,'maximumStageRate',maximumRate, ...
    'realizedMaximumStageCfl',h*maximumRate,'configuredCfl',policy.cfl, ...
    'cflStepFormulaApplied',true,'cflStepFormula','min(maxDt,0.95*cfl/currentRate,remainingTau)', ...
    'stageRates',records.rates,'stageCfl',h*records.rates, ...
    'stageAbscissae',records.abscissae,'stageCount',stages, ...
    'maximumCflStageIndex',index,'cflPosition',out.rhs.cflPosition, ...
    'maximumStageCflPosition',records.positions{index}, ...
    'stagePeakDiagnostics',records.peaks,'allAttemptPeakDiagnostics',peaks, ...
    'stageGaugeResidualInfinity',records.gaugeResidualInfinity, ...
    'maximumStageGaugeResidual',max(records.gaugeResidualInfinity), ...
    'stagePhysicalClockRates',records.clockRates,'rhsEvaluations',calls, ...
    'initialRhsReused',true,'terminalRhsEvaluated',true,'terminalRhsReusedNextStep',true, ...
    'cflScope','initial interval from current RHS; ordinary spline-knot partition; stage CFL recorded without rejection', ...
    'physicalClockIncrement',records.physicalClockIncrement,'physicalClockUpdate',clock, ...
    'compensatedPhysicalClock',true,'originAmplitudeAdvancedWithSameTableau',true, ...
    'sspMethod',strcmp(method,'ssprk3'),'analyticBackgroundVariableChange',hasResidual);

    function [value,record]=trial_at(interval)
        probeH=interval;
        [value,record]=ipm_cr.time.trial(state,interval,method,ctx.mesh,policy,callbacks);
    end
    function [f,W]=stage(field,z,previous,wantDiagnostics,increment)
        calls=calls+1;
        staged=state;staged.scale=ipm_cr.time.stageScale(state.scale,z,increment);staged.rhs=previous;
        if hasResidual,staged.eta=field;else,staged.W=field;end
        [f,W]=ipm_cr.rhs(staged,ctx,wantDiagnostics);
        note_peak(f,z(5),(z(5)-state.scale.canonicalTime)/probeH,false,true);
    end
    function s=finish(old,field,W,z,f,step,increment)
        s=old;s.W=W;s.rhs=f;
        s.scale=ipm_cr.time.stageScale(old.scale,z,increment);
        s.scale.canonicalTime=old.scale.canonicalTime+step;
        if hasResidual,s.eta=field;end
        s.step=old.step+1;
    end
    function note_peak(f,tau,abscissa,reused,returned)
        p=ipm_cr.time.peak(f);p.canonicalTime=tau;p.abscissa=abscissa;
        p.initialRhsReused=reused;p.rhsReturned=returned;
        if isempty(peaks),peaks=p;else,peaks(end+1)=p;end
    end
end

function [out,records,h,p]=fit_knot(initial,out,records,h,knots,evaluate,p)
% Only a unique, locally concave peak crossing adjacent spline cells defines
% this partition. A disconnected maximizing-branch change evolves directly.
p.kind='RK4_original_spline_knot_fixed4_secant';
a=ipm_cr.time.peak(initial.rhs);b=ipm_cr.time.peak(out.rhs);
direction=sign(b.x-a.x);p.direction=direction;
if ~a.available||~b.available||direction==0,return;end
[ca,p.startClassification]=peak_cell(a.x,knots,direction,true);
[cb,p.endClassification]=peak_cell(b.x,knots,direction,false);
if cb-ca~=direction||any([records.peaks.tieNonunique])|| ...
        a.curvature>=0||b.curvature>=0,return;end
node=ca+(direction>0);[v0,c0]=knot_slope(initial.rhs,node,knots);
[v1,c1]=knot_slope(out.rhs,node,knots);
if ~(c0<0&&c1<0&&v0*v1<=0),return;end
p.fitted=true;p.knotIndex=node;p.knotX=knots(node);
previousH=0;previousValue=v0;value=v1;
for iteration=1:4
    if value==0,break;end
    nextH=h-value*(h-previousH)/(value-previousValue);
    if nextH==h,break;end
    previousH=h;previousValue=value;h=nextH;
    [out,records]=evaluate(h);value=knot_slope(out.rhs,node,knots);
    p.secantIterations=iteration;
end
p.segmentDt=h;p.knotSlope=value;p.rawPeakKnotDefect=out.rhs.measurement.peak.x-knots(node);
end
function [cell,info]=peak_cell(value,knots,direction,isStart)
% Coordinate identity only: outgoing ownership at a start knot, incoming
% ownership at an end knot. Neither the measured peak nor state is changed.
knots=knots(:)';[distance,node]=min(abs(knots-value));
tolerance=256*eps(max([1,abs(value),abs(knots(node))]));onKnot=distance<=tolerance;
classified=value;
if onKnot
    classified=knots(node);cell=node-(direction<0);
    if ~isStart,cell=node-(direction>0);end
else,cell=find(knots<value,1,'last');end
cell=max(1,min(numel(knots)-1,cell));
info=struct('observedX',value,'classifiedX',classified,'coordinateDifference',classified-value, ...
    'onKnot',onKnot,'tolerance',tolerance,'cellIndex',cell);
end
function [slope,curvature]=knot_slope(f,node,knots)
q=f.measurement.q;peak=f.measurement.peak;
if isfield(peak,'evaluationAnchor'),q=q-peak.evaluationAnchor;end
M=peak.splineSecondDerivatives;width=knots(node+1)-knots(node);
slope=(q(node+1)-q(node))/width-width*(2*M(node)+M(node+1))/6;
curvature=M(node);
end
function z=pack(s)
z=[s.logC_l;s.logC_omega;s.physicalTime;s.X_shift;s.canonicalTime;s.densityOrigin];
end
