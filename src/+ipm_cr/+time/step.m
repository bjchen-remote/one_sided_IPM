function [out,info]=step(state,ctx,finalTau,controls)
%STEP Transactional complete-RHS explicit RK4, with legacy SSPRK3 also available.
% CONTROLS contains only cfl, maxDt, minDt and timeIntegrator overrides.
% Missing controls inherit the prepared context's time policy. Every stage and
% fresh endpoint use the complete W -> q -> Poisson -> I/P/C -> WENO RHS.
% The identical accepted initial RHS is reused, then the fresh endpoint is
% retained for the next step. All stages AND the endpoint obey the CFL.
% The stated order is temporal only; no spatial-order certificate follows.
% With state.eta, integrate that exact background-subtracted variable
% using rhs.residualRhs. rhs reconstructs full centered W at each scalar
% stage; state.W and rhs.F retain their full-W meanings.
% State/RHS context identity is checked before the cached initial RHS is used.
% RK4 steps are split at observed original-spline peak-cell crossings.
% Event probes use the same complete explicit tableau from the accepted
% initial state; the original peak functional is never altered or projected.
timer=tic;
assert(nargin>=3&&isa(ctx,'ipm_cr.context'),'ipm_cr:StateContext', ...
    'Use the prepared context that owns the accepted state.');
g=ctx.mesh;wallOps=ctx.wall;
key=ipm_cr.internal.key(g,wallOps,ctx.opts,ctx.identity);
assert(isfield(state,'ctxKey')&&isequaln(state.ctxKey,key), ...
    'ipm_cr:StateContext','The accepted RHS belongs to another geometry or method. Reinitialize or refresh explicitly.');
assert(isfield(state,'ctxId')&&strcmp(state.ctxId,ctx.id)&& ...
    isequaln(key,ctx.key),'ipm_cr:StateContext', ...
    'The accepted state and step must use the same run context.');
if nargin<4,controls=struct();end
assert(isstruct(controls)&&isscalar(controls),'ipm_cr:StepOptions','Use scalar time controls.');
names=fieldnames(controls);allowed={'cfl','maxDt','minDt','timeIntegrator'};
assert(all(ismember(names,allowed)),'ipm_cr:StepOptions', ...
    'Only time controls can change on an existing context.');
policy=struct();
for k=1:numel(allowed)
    name=allowed{k};if isfield(ctx.opts,name),policy.(name)=ctx.opts.(name);end
end
for k=1:numel(names),policy.(names{k})=controls.(names{k});end
controls=policy;
method='rk4';if isfield(controls,'timeIntegrator'),method=char(controls.timeIntegrator);end
assert(ismember(method,{'rk4','ssprk3'}),'ipm_cr:ExplicitTimeIntegrator', ...
    'Use rk4 or the explicitly selected legacy ssprk3 method.');
for name={'cfl','maxDt','minDt'}
    validateattributes(controls.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(finalTau,{'numeric'},{'scalar','real','finite'});
validate_rate(state.rhs.rate);
initial=state;initialZ=pack(state.scale);
hasResidual=isfield(state,'eta');
assert(all(isfinite(initialZ))&&isfinite(state.scale.physicalTimeLow), ...
    'ipm_cr:ExplicitScale','The accepted initial scale and clock must be finite.');
remaining=finalTau-state.scale.canonicalTime;
assert(remaining>0,'ipm_cr:ExplicitFinalTime','The endpoint must follow the accepted state.');
% The predictor margin never replaces the binding complete-stage audit.
dt=min([controls.maxDt,.95*controls.cfl/max(state.rhs.rate,realmin),remaining]);
retries=0;rhsEvaluations=0;retryReasons={};
retryDiagnostics=struct([]);retryCauses={};maximumRetries=12;
peakDiagnostics=struct([]);
currentPeakProbe=0;
wallKnots=[];
if isstruct(wallOps)&&isfield(wallOps,'x'),wallKnots=wallOps.x(:);end
currentStageTau=state.scale.canonicalTime;
callbacks=struct('stage',@stage,'finish',@finish,'notePeak',@note_peak);
while true
    if dt<controls.minDt||initial.scale.canonicalTime+dt<=initial.scale.canonicalTime
        retry_error(MException('ipm_cr:ExplicitMinimumStep', ...
            'No admissible explicit step remains: dt %.17g, minDt %.17g.',dt,controls.minDt));
    end
    before=rhsEvaluations;
    currentStageTau=initial.scale.canonicalTime;
    currentPeakProbe=0;peakEvent=empty_peak_event();
    try
        [trial,records]=ipm_cr.time.trial(initial,dt,method,g,controls,callbacks);
        if strcmp(method,'rk4')&&~isempty(trial)&&inside_cfl(max(records.rates),dt)
            eventCallsBefore=rhsEvaluations;
            [trial,records,dt,peakEvent]=ipm_cr.time.event( ...
                initial,trial,records,dt,wallKnots,controls,@event_probe,peakEvent);
            peakEvent.extraRhsEvaluations=rhsEvaluations-eventCallsBefore;
        end
    catch problem
        if ~recoverable(problem)
            if isempty(retryCauses),rethrow(problem);else,retry_error(problem);end
        end
        nextDt=dt/2;
        note_retry(problem,dt,nextDt,rhsEvaluations-before);
        if retries>=maximumRetries,retry_error(problem);end
        retryReasons{end+1}=problem.identifier; %#ok<AGROW>
        dt=nextDt;retries=retries+1;continue;
    end
    maximumRate=max(records.rates);
    if ~isempty(trial)&&inside_cfl(maximumRate,dt)
        out=trial;acceptedAttemptRhsEvaluations=rhsEvaluations-before;break;
    end
    nextDt=min(.95*dt,.95*controls.cfl/maximumRate);
    problem=MException('ipm_cr:ExplicitStageCfl', ...
        'Complete stage CFL %.17g exceeds configured bound %.17g at dt %.17g (maximum rate %.17g).', ...
        dt*maximumRate,controls.cfl,dt,maximumRate);
    note_retry(problem,dt,nextDt,rhsEvaluations-before);
    if retries>=maximumRetries,retry_error(problem);end
    retryReasons{end+1}=problem.identifier; %#ok<AGROW>
    dt=nextDt;retries=retries+1;
end
[~,maximumStageIndex]=max(records.rates);
if strcmp(method,'rk4'),order=4;stages=4;label='rk4_full_hybrid_cfl_peak_events_v2';
else,order=3;stages=3;label='ssprk3_full_hybrid_cfl_v1';end
info=struct('method',label,'timeIntegrator',method,'temporalOrder',order, ...
    'spatialOrderCertified',false,'dt',dt,'seconds',toc(timer), ...
    'rate',initial.rhs.rate,'maximumStageRate',maximumRate, ...
    'realizedMaximumStageCfl',dt*maximumRate,'configuredCfl',controls.cfl, ...
    'retries',retries,'retryReasons',{retryReasons},'retryDiagnostics',retryDiagnostics, ...
    'maximumStepRetries',maximumRetries,'explicitCflLimitApplied',true, ...
    'cflPosition',out.rhs.cflPosition, ...
    'maximumStageCflPosition',records.positions{maximumStageIndex}, ...
    'maximumCflStageIndex',maximumStageIndex,'stageCount',stages, ...
    'stageRates',records.rates,'stageCfl',dt*records.rates, ...
    'stageAbscissae',records.abscissae, ...
    'stagePeakDiagnostics',records.peaks,'allAttemptPeakDiagnostics',peakDiagnostics, ...
    'peakDiagnosticsAreReadOnly',true, ...
    'peakEvent',peakEvent,'originalPeakCellEventsEnabled',strcmp(method,'rk4'), ...
    'stageGaugeResidualInfinity',records.gaugeResidualInfinity, ...
    'stagePhysicalClockRates',records.clockRates, ...
    'maximumStageGaugeResidual',max(records.gaugeResidualInfinity), ...
    'rhsEvaluations',rhsEvaluations,'acceptedAttemptRhsEvaluations',acceptedAttemptRhsEvaluations, ...
    'initialRhsReused',true,'terminalRhsEvaluated',true,'terminalRhsReusedNextStep',true, ...
    'cflScope','cached initial RHS, complete intermediate stages and fresh endpoint', ...
    'physicalClockIncrement',records.physicalClockIncrement,'compensatedPhysicalClock',true, ...
    'originAmplitudeAdvancedWithSameTableau',true,'sspMethod',strcmp(method,'ssprk3'), ...
    'analyticBackgroundVariableChange',hasResidual);

    function [f,rhoActual]=stage(field,z,previous,wantDiagnostics,physicalIncrement)
        if nargin<4,wantDiagnostics=false;end
        currentStageTau=z(5);
        speed=exp(z(2)-z(1));
        assert(all(isfinite(z))&&isfinite(speed)&&speed>0, ...
            'ipm_cr:ExplicitScale','Every explicit scalar stage and clock speed must be finite.');
        % evaluate synchronizes its density; RK4 avoids a redundant exchange.
        rhsEvaluations=rhsEvaluations+1;
        peakIndex=note_peak([],currentStageTau, ...
            (currentStageTau-initial.scale.canonicalTime)/dt,false,false);
        staged=initial;staged.scale=ipm_cr.time.stageScale(initial.scale,z,physicalIncrement);staged.rhs=previous;
        if hasResidual,staged.eta=field;else,staged.W=field;end
        [f,rhoActual]=ipm_cr.rhs(staged,ctx,wantDiagnostics);
        snapshot=ipm_cr.time.peak(f);names=fieldnames(snapshot);
        for k=1:numel(names),peakDiagnostics(peakIndex).(names{k})=snapshot.(names{k});end
        peakDiagnostics(peakIndex).rhsReturned=true;
        validate_rate(f.rate);
    end
    function s=finish(old,q,rhoActual,z,f,h,increment)
        assert(all(isfinite(z))&&isfinite(increment)&&increment>0, ...
            'ipm_cr:ExplicitScale','The scalar state and weighted clock increment must be finite.');
        s=old;s.W=rhoActual;s.rhs=f;s.scale=ipm_cr.time.stageScale(old.scale,z,increment);
        if hasResidual,s.eta=q;end
        s.scale.canonicalTime=old.scale.canonicalTime+h;
        assert(isfinite(s.scale.physicalTime)&&isfinite(s.scale.physicalTimeLow)&& ...
            (s.scale.physicalTime>old.scale.physicalTime|| ...
            (s.scale.physicalTime==old.scale.physicalTime&&s.scale.physicalTimeLow>old.scale.physicalTimeLow)), ...
            'ipm_cr:ExplicitClock','The compensated physical clock must remain finite and increase.');
        s.step=old.step+1;
        assert(all(isfinite(q.core),'all')&&all(isfinite(q.polar),'all')&& ...
            all(isfinite(rhoActual.core),'all')&&all(isfinite(rhoActual.polar),'all')&& ...
            rhoActual.core(1,g.core.originIndex)==0,'ipm_cr:ExplicitState', ...
            'A complete accepted state must remain finite and centered.');
    end
    function yes=inside_cfl(rate,h)
        yes=h*rate<=controls.cfl*(1+64*eps);
    end
    function [probedState,probedRecords]=event_probe(h,probeIndex)
        savedDt=dt;savedProbe=currentPeakProbe;
        dt=h;currentPeakProbe=probeIndex;
        try
            [probedState,probedRecords]=ipm_cr.time.trial(initial,h,'rk4',g,controls,callbacks);
        catch problem
            dt=savedDt;currentPeakProbe=savedProbe;rethrow(problem);
        end
        dt=savedDt;currentPeakProbe=savedProbe;
        if isempty(probedState)||~inside_cfl(max(probedRecords.rates),h)
            error('ipm_cr:ExplicitStageCfl', ...
                'A complete peak-event RK4 probe exceeded the binding stage CFL.');
        end
    end
    function index=note_peak(f,tau,abscissa,reused,returned)
        entry=ipm_cr.time.peak(f);
        entry.attempt=retries+1;entry.canonicalTime=tau;entry.abscissa=abscissa;
        entry.initialRhsReused=reused;entry.rhsEvaluation=rhsEvaluations;
        entry.peakEventProbe=currentPeakProbe;
        entry.rhsReturned=returned;
        if isempty(peakDiagnostics),peakDiagnostics=entry;
        else,peakDiagnostics(end+1)=entry;end %#ok<AGROW>
        index=numel(peakDiagnostics);
    end
    function note_retry(problem,h,nextDt,calls)
        entry=struct('attempt',retries+1,'dt',h,'nextDt',nextDt, ...
            'stageCanonicalTime',currentStageTau,'rhsCallsInAttempt',calls, ...
            'identifier',problem.identifier,'message',problem.message,'stack',problem.stack);
        if isempty(retryDiagnostics),retryDiagnostics=entry;
        else,retryDiagnostics(end+1)=entry;end %#ok<AGROW>
        retryCauses{end+1}=problem;
    end
    function retry_error(problem)
        text=retry_text(retryDiagnostics);
        failure=MException(problem.identifier, ...
            'Explicit step stopped after %d retries; last dt %.17g. %s Retry diagnostics: %s', ...
            retries,dt,problem.message,text);
        for causeIndex=1:numel(retryCauses)
            failure=addCause(failure,retryCauses{causeIndex});
        end
        if isempty(retryCauses)||~isequal(problem,retryCauses{end})
            failure=addCause(failure,problem);
        end
        throw(failure);
    end
end

function text=retry_text(records)
if isempty(records),text='none';return;end
items=cell(1,numel(records));
for j=1:numel(records)
    r=records(j);
    items{j}=sprintf('[attempt %d dt %.17g tau %.17g %s: %s]', ...
        r.attempt,r.dt,r.stageCanonicalTime,r.identifier,r.message);
end
text=strjoin(items,' ');
end

function e=empty_peak_event()
e=struct('limited',false,'located',false,'knotX',NaN,'direction',0, ...
    'rootResidual',NaN,'probeCount',0,'extraRhsEvaluations',0, ...
    'initialProposedDt',NaN,'acceptedDt',NaN,'rootTimeBracket',[NaN,NaN], ...
    'stageGuardCuts',0,'status','no_observed_peak_cell_crossing');
end

function yes=recoverable(problem)
yes=any(strcmp(problem.identifier,{'polar:GaugeRootUnresolved','polar:GaugeFullRhsAudit', ...
    'polar:WallDegenerate','ipm_cr:ExplicitGaugeAudit','ipm_cr:ExplicitGaugeRoot','ipm_cr:ExplicitGaugePeak','ipm_cr:AccuracyGate', ...
    'ipm_cr:ExplicitScale','ipm_cr:ExplicitClock','ipm_cr:ExplicitState', ...
    'ipm_cr:ExplicitTransportNonfinite','ipm_cr:ExplicitSynchronization', ...
    'ipm_cr:ExplicitStageCfl'}));
end
function validate_rate(rate)
assert(isnumeric(rate)&&isreal(rate)&&isscalar(rate)&&isfinite(rate)&&rate>=0, ...
    'ipm_cr:ExplicitRate','The explicit transport plus reaction rate must be finite and nonnegative.');
end
function z=pack(s)
z=[s.logC_l;s.logC_omega;s.physicalTime;s.X_shift;s.canonicalTime;s.densityOrigin];
end
