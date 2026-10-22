function [out,info]=step(state,g,wallOps,options,finalTau,ctx)
%STEP Transactional complete-RHS explicit RK4, with legacy SSPRK3 also available.
% OPTIONS.timeIntegrator is 'rk4' (default) or 'ssprk3'. Every stage and
% fresh endpoint use the complete W -> q -> Poisson -> I/P/C -> WENO RHS.
% The identical accepted initial RHS is reused, then the fresh endpoint is
% retained for the next step. All stages AND the endpoint obey the CFL.
% The stated order is temporal only; no spatial-order certificate follows.
% With state.residual, integrate that exact background-subtracted variable
% using flow.residualRhs. evolve.residual reconstructs full centered W at
% each scalar stage; state.rho and flow.F retain their full-W meanings.
% State/RHS context identity is checked before the cached initial RHS is used.
% RK4 steps are split at observed original-spline peak-cell crossings.
% Event probes use the same complete explicit tableau from the accepted
% initial state; the original peak functional is never altered or projected.
timer=tic;
key=ipm_cr.internal.key(g,wallOps,options);
assert(isfield(state,'ctxKey')&&isequaln(state.ctxKey,key), ...
    'ipm_cr:StateContext','The accepted RHS belongs to another geometry or method. Reinitialize or refresh explicitly.');
if nargin<6||isempty(ctx)
    ctx=ipm_cr.context(g,wallOps,options);
    state.ctxId=ctx.id;
else
    assert(isfield(state,'ctxId')&&strcmp(state.ctxId,ctx.id)&& ...
        isequaln(key,ctx.key),'ipm_cr:StateContext', ...
        'The accepted state and step must use the same run context.');
end
method='rk4';if isfield(options,'timeIntegrator'),method=char(options.timeIntegrator);end
assert(ismember(method,{'rk4','ssprk3'}),'ipm_cr:ExplicitTimeIntegrator', ...
    'Use rk4 or the explicitly selected legacy ssprk3 method.');
for name={'cfl','maxDt','minDt'}
    validateattributes(options.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(finalTau,{'numeric'},{'scalar','real','finite'});
validate_rate(state.flow.rate);
initial=state;initialZ=pack(state.scale);
hasResidual=isfield(state,'residual');
assert(all(isfinite(initialZ))&&isfinite(state.scale.physicalTimeLow), ...
    'ipm_cr:ExplicitScale','The accepted initial scale and clock must be finite.');
remaining=finalTau-state.scale.canonicalTime;
assert(remaining>0,'ipm_cr:ExplicitFinalTime','The endpoint must follow the accepted state.');
% The predictor margin never replaces the binding complete-stage audit.
dt=min([options.maxDt,.95*options.cfl/max(state.flow.rate,realmin),remaining]);
retries=0;rhsEvaluations=0;retryReasons={};
retryDiagnostics=struct([]);retryCauses={};maximumRetries=12;
peakDiagnostics=struct([]);
currentPeakProbe=0;
wallKnots=[];
if isstruct(wallOps)&&isfield(wallOps,'x'),wallKnots=wallOps.x(:);end
currentStageTau=state.scale.canonicalTime;
callbacks=struct('stage',@stage,'finish',@finish,'notePeak',@note_peak);
while true
    if dt<options.minDt||initial.scale.canonicalTime+dt<=initial.scale.canonicalTime
        retry_error(MException('ipm_cr:ExplicitMinimumStep', ...
            'No admissible explicit step remains: dt %.17g, minDt %.17g.',dt,options.minDt));
    end
    before=rhsEvaluations;
    currentStageTau=initial.scale.canonicalTime;
    currentPeakProbe=0;peakEvent=empty_peak_event();
    try
        [trial,records]=ipm_cr.time.trial(initial,dt,method,g,options,callbacks);
        if strcmp(method,'rk4')&&~isempty(trial)&&inside_cfl(max(records.rates),dt)
            eventCallsBefore=rhsEvaluations;
            [trial,records,dt,peakEvent]=ipm_cr.time.event( ...
                initial,trial,records,dt,wallKnots,options,@event_probe,peakEvent);
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
    nextDt=min(.95*dt,.95*options.cfl/maximumRate);
    problem=MException('ipm_cr:ExplicitStageCfl', ...
        'Complete stage CFL %.17g exceeds configured bound %.17g at dt %.17g (maximum rate %.17g).', ...
        dt*maximumRate,options.cfl,dt,maximumRate);
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
    'rate',initial.flow.rate,'maximumStageRate',maximumRate, ...
    'realizedMaximumStageCfl',dt*maximumRate,'configuredCfl',options.cfl, ...
    'retries',retries,'retryReasons',{retryReasons},'retryDiagnostics',retryDiagnostics, ...
    'maximumStepRetries',maximumRetries,'explicitCflLimitApplied',true, ...
    'cflPosition',out.flow.cflPosition, ...
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

    function [f,rhoActual]=stage(field,z,previous,wantDiagnostics)
        if nargin<4,wantDiagnostics=false;end
        currentStageTau=z(5);
        speed=exp(z(2)-z(1));
        assert(all(isfinite(z))&&isfinite(speed)&&speed>0, ...
            'ipm_cr:ExplicitScale','Every explicit scalar stage and clock speed must be finite.');
        % evaluate synchronizes its density; RK4 avoids a redundant exchange.
        rhsEvaluations=rhsEvaluations+1;
        peakIndex=note_peak([],currentStageTau, ...
            (currentStageTau-initial.scale.canonicalTime)/dt,false,false);
        if hasResidual
            stageScale=unpack(z,initial.scale);
            [f,rhoActual]=ipm_cr.evolve.residual( ...
                field,stageScale,g,wallOps,options,[previous.c_l,previous.c_r],wantDiagnostics,ctx);
        else
            [f,rhoActual]=ipm_cr.evolve.rhs( ...
                field,g,wallOps,options,[previous.c_l,previous.c_r],[],wantDiagnostics,ctx);
        end
        snapshot=ipm_cr.time.peak(f);names=fieldnames(snapshot);
        for k=1:numel(names),peakDiagnostics(peakIndex).(names{k})=snapshot.(names{k});end
        peakDiagnostics(peakIndex).rhsReturned=true;
        validate_rate(f.rate);
    end
    function s=finish(old,q,rhoActual,z,f,h,increment)
        assert(all(isfinite(z))&&isfinite(increment)&&increment>0, ...
            'ipm_cr:ExplicitScale','The scalar state and weighted clock increment must be finite.');
        s=old;s.rho=rhoActual;s.flow=f;s.scale=unpack(z,old.scale);
        if hasResidual,s.residual=q;end
        s.scale.canonicalTime=old.scale.canonicalTime+h;
        [s.scale.physicalTime,s.scale.physicalTimeLow]=clock_add( ...
            old.scale.physicalTime,old.scale.physicalTimeLow,increment);
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
        yes=h*rate<=options.cfl*(1+64*eps);
    end
    function [probedState,probedRecords]=event_probe(h,probeIndex)
        savedDt=dt;savedProbe=currentPeakProbe;
        dt=h;currentPeakProbe=probeIndex;
        try
            [probedState,probedRecords]=ipm_cr.time.trial(initial,h,'rk4',g,options,callbacks);
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
    'polar:WallDegenerate','ipm_cr:ExplicitGaugeAudit','ipm_cr:ExplicitGaugeRoot','ipm_cr:AccuracyGate', ...
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
function s=unpack(z,s)
s.logC_l=z(1);s.logC_omega=z(2);s.physicalTime=z(3);s.X_shift=z(4);s.canonicalTime=z(5);s.densityOrigin=z(6);
end
function [high,low]=clock_add(high,low,increment)
% Preserve the historical compensated positive-increment clock arithmetic.
s=high+increment;b=s-high;e=(high-(s-b))+(increment-b);
t=low+e;high=s+t;low=(s-high)+t;
end
