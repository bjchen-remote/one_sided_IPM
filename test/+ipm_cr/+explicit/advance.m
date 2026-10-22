function [out,info]=advance(state,g,wallOps,options,finalTau)
%ADVANCE Complete-RHS explicit RK4, with legacy SSPRK3 also available.
% OPTIONS.timeIntegrator is 'rk4' (default) or 'ssprk3'. Every stage and
% fresh endpoint use the complete W -> q -> Poisson -> I/P/C -> WENO RHS.
% The identical accepted initial RHS is reused, then the fresh endpoint is
% retained for the next step. All stages AND the endpoint obey the CFL.
% The stated order is temporal only; no spatial-order certificate follows.
% With state.residual, integrate that exact background-subtracted variable
% using flow.residualRhs. evaluateResidual reconstructs full centered W at
% each scalar stage; state.rho and flow.F retain their full-W meanings.
% RK4 steps are split at observed original-spline peak-cell crossings.
% Event probes use the same complete explicit tableau from the accepted
% initial state; the original peak functional is never altered or projected.
timer=tic;
method='rk4';if isfield(options,'timeIntegrator'),method=char(options.timeIntegrator);end
assert(ismember(method,{'rk4','ssprk3'}),'ipm_cr:ExplicitTimeIntegrator', ...
    'Use rk4 or the explicitly selected legacy ssprk3 method.');
for name={'cfl','maxDt','minDt'}
    validateattributes(options.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(finalTau,{'numeric'},{'scalar','real','finite'});
validate_rate(state.flow.rate);
initial=state;initialZ=pack(state.scale);
hasResidual=isfield(state,'residual');initialGrid=state.rho;
if hasResidual,initialGrid=state.residual;end
assert(all(isfinite(initialZ))&&isfinite(state.scale.physicalTimeLow), ...
    'ipm_cr:ExplicitScale','The accepted initial scale and clock must be finite.');
remaining=finalTau-state.scale.canonicalTime;
assert(remaining>0,'ipm_cr:ExplicitFinalTime','The endpoint must follow the accepted state.');
% The predictor margin never replaces the binding complete-stage audit.
dt=min([options.maxDt,.95*options.cfl/max(state.flow.rate,realmin),remaining]);
retries=0;rhsEvaluations=0;retryReasons={};tableau=[];
retryDiagnostics=struct([]);retryCauses={};maximumRetries=12;
peakDiagnostics=struct([]);
currentPeakProbe=0;
wallKnots=[];
if isstruct(wallOps)&&isfield(wallOps,'x'),wallKnots=wallOps.x(:);end
currentStageTau=state.scale.canonicalTime;
if strcmp(method,'rk4'),tableau=ipm_cr.explicit.rk4Tableau();end
while true
    if dt<options.minDt||initial.scale.canonicalTime+dt<=initial.scale.canonicalTime
        retry_error(MException('ipm_cr:ExplicitMinimumStep', ...
            'No admissible explicit step remains: dt %.17g, minDt %.17g.',dt,options.minDt));
    end
    before=rhsEvaluations;
    currentStageTau=initial.scale.canonicalTime;
    currentPeakProbe=0;peakEvent=empty_peak_event();
    try
        if strcmp(method,'rk4'),[trial,records]=attempt_rk4(initial,dt);
        else,[trial,records]=attempt_ssprk3(initial,dt);end
        if strcmp(method,'rk4')&&~isempty(trial)&&inside_cfl(max(records.rates),dt)
            [trial,records,dt,peakEvent]=limit_peak_event(trial,records,dt);
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

    function [s,records]=attempt_rk4(old,h)
        s=[];z=initialZ;f0=old.flow;k0=scalar_rate(z,f0);
        note_peak(f0,z(5),0,true,true);
        records=record_start(f0,k0(3));r0=grid_rhs(f0);
        q1=add(initialGrid,1,r0,h*tableau.A(2,1));z1=z+h*tableau.A(2,1)*k0;
        f1=stage(q1,z1,f0);k1=scalar_rate(z1,f1);
        records=record_add(records,f1,k1(3),tableau.c(2));
        if ~inside_cfl(f1.rate,h),return;end
        r1=grid_rhs(f1);
        q2=add(initialGrid,1,r1,h*tableau.A(3,2));z2=z+h*tableau.A(3,2)*k1;
        f2=stage(q2,z2,f1);k2=scalar_rate(z2,f2);
        records=record_add(records,f2,k2(3),tableau.c(3));
        if ~inside_cfl(f2.rate,h),return;end
        r2=grid_rhs(f2);
        q3=add(initialGrid,1,r2,h*tableau.A(4,3));z3=z+h*tableau.A(4,3)*k2;
        f3=stage(q3,z3,f2);k3=scalar_rate(z3,f3);
        records=record_add(records,f3,k3(3),tableau.c(4));
        if ~inside_cfl(f3.rate,h),return;end
        b=tableau.b;r3=grid_rhs(f3);
        qNew=struct('core',initialGrid.core+h*(b(1)*r0.core+b(2)*r1.core+ ...
            b(3)*r2.core+b(4)*r3.core), ...
            'polar',initialGrid.polar+h*(b(1)*r0.polar+b(2)*r1.polar+ ...
            b(3)*r2.polar+b(4)*r3.polar));
        zNew=z+h*(b(1)*k0+b(2)*k1+b(3)*k2+b(4)*k3);
        qNew=ipm_cr.explicit.synchronize(qNew,g);
        [fEnd,rhoActual]=stage(qNew,zNew,f3,true);
        records=record_add(records,fEnd,exp(zNew(2)-zNew(1)),1);
        if ~inside_cfl(fEnd.rate,h),return;end
        increment=h*(b(1)*k0(3)+b(2)*k1(3)+b(3)*k2(3)+b(4)*k3(3));
        s=finish(old,qNew,rhoActual,zNew,fEnd,h,increment);records.physicalClockIncrement=increment;
    end
    function [s,records]=attempt_ssprk3(old,h)
        % Preserve the historical SSPRK3 density/scalar arithmetic.
        s=[];z=initialZ;f0=old.flow;k0=scalar_rate(z,f0);
        note_peak(f0,z(5),0,true,true);
        records=record_start(f0,k0(3));
        q1=add(initialGrid,1,grid_rhs(f0),h);z1=z+h*k0;
        q1=ipm_cr.explicit.synchronize(q1,g);
        f1=stage(q1,z1,f0);k1=scalar_rate(z1,f1);
        records=record_add(records,f1,k1(3),1);
        if ~inside_cfl(f1.rate,h),return;end
        q2=add(initialGrid,.75,add(q1,1,grid_rhs(f1),h),.25);z2=.75*z+.25*(z1+h*k1);
        q2=ipm_cr.explicit.synchronize(q2,g);
        f2=stage(q2,z2,f1);k2=scalar_rate(z2,f2);
        records=record_add(records,f2,k2(3),.5);
        if ~inside_cfl(f2.rate,h),return;end
        qNew=add(initialGrid,1/3,add(q2,1,grid_rhs(f2),h),2/3);zNew=(z+2*(z2+h*k2))/3;
        qNew=ipm_cr.explicit.synchronize(qNew,g);
        [fEnd,rhoActual]=stage(qNew,zNew,f2,true);
        records=record_add(records,fEnd,exp(zNew(2)-zNew(1)),1);
        if ~inside_cfl(fEnd.rate,h),return;end
        increment=h*(k0(3)/6+k1(3)/6+2*k2(3)/3);
        s=finish(old,qNew,rhoActual,zNew,fEnd,h,increment);records.physicalClockIncrement=increment;
    end
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
            [f,rhoActual]=ipm_cr.explicit.evaluateResidual( ...
                field,stageScale,g,wallOps,options,[previous.c_l,previous.c_r],wantDiagnostics);
        else
            [f,rhoActual]=ipm_cr.explicit.evaluate( ...
                field,g,wallOps,options,[previous.c_l,previous.c_r],[],wantDiagnostics);
        end
        snapshot=peak_snapshot(f);names=fieldnames(snapshot);
        for k=1:numel(names),peakDiagnostics(peakIndex).(names{k})=snapshot.(names{k});end
        peakDiagnostics(peakIndex).rhsReturned=true;
        validate_rate(f.rate);
    end
    function r=grid_rhs(f)
        if hasResidual,r=f.residualRhs;else,r=f.F;end
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
    function [s,r,h,event]=limit_peak_event(s,r,h)
        % Only event searching changes H; all trial states still originate
        % at INITIAL and are accepted by the ordinary outer-step gate.
        event=empty_peak_event();event.initialProposedDt=h;event.acceptedDt=h;
        callsBefore=rhsEvaluations;probes=0;guardCuts=0;
        start=peak_snapshot(initial.flow);
        if isempty(wallKnots)||~start.available||~isfinite(start.x),return;end
        while true
            crossing=find_crossing(r,h,start.x);
            if ~crossing.found,break;end
            event.knotX=crossing.knot;event.direction=crossing.direction;
            xtol=64*eps(max(1,abs(crossing.knot)));
            endpoint=peak_snapshot(s.flow);
            fStart=crossing.direction*(start.x-crossing.knot);
            fEnd=crossing.direction*(endpoint.x-crossing.knot);
            if abs(fEnd)<=xtol
                event.located=true;event.rootResidual=endpoint.x-crossing.knot;
                event.status='located_peak_cell_boundary';break;
            end
            if fEnd<0
                % An internal RK prediction crossed, but its true substep
                % endpoint did not bracket that crossing. A half step is a
                % guard, not a falsely reported located physical event.
                h=h/2;guardCuts=guardCuts+1;
                if h<options.minDt||initial.scale.canonicalTime+h<=initial.scale.canonicalTime
                    error('ipm_cr:PeakEventMinimumStep', ...
                        'A peak-cell stage guard requires dt below minDt.');
                end
                if guardCuts>24
                    error('ipm_cr:PeakEventUnresolved','Peak-cell stage guards did not resolve a crossing.');
                end
                [s,r]=probe(h);event.status='stage_crossing_guard';
                event.rootResidual=s.flow.measurement.peak.x-crossing.knot;
                continue;
            end
            % The endpoint brackets the first observed knot in this
            % direction. Locate endpoint.x(h)=knot using complete RK4
            % substeps, with safeguarded secant/bisection interpolation.
            lo=0;hi=h;flo=fStart;fhi=fEnd;
            bestState=s;bestRecords=r;bestH=h;bestResidual=fEnd;
            located=false;slowProgress=0;
            for iteration=1:64
                if abs(bestResidual)<=xtol,located=true;break;end
                timeTolerance=16*eps(max(1,abs(initial.scale.canonicalTime)+hi));
                if hi-lo<=timeTolerance
                    if abs(bestResidual)<=xtol,located=true;end
                    break;
                end
                width=hi-lo;probeH=lo-flo*width/(fhi-flo);
                % A nearly linear peak trajectory legitimately puts its
                % next secant very near a bracket endpoint. Preserve it
                % unless rounding leaves the bracket, or two consecutive
                % probes failed to reduce the best residual by 10 percent.
                if ~isfinite(probeH)||probeH<=lo||probeH>=hi||slowProgress>=2
                    probeH=lo+width/2;slowProgress=0;
                end
                if probeH<options.minDt||initial.scale.canonicalTime+probeH<=initial.scale.canonicalTime
                    error('ipm_cr:PeakEventMinimumStep', ...
                        'A located peak-cell event requires dt below minDt.');
                end
                [candidateState,candidateRecords]=probe(probeH);
                snapshot=peak_snapshot(candidateState.flow);
                value=crossing.direction*(snapshot.x-crossing.knot);
                if ~isfinite(value)
                    error('ipm_cr:PeakEventUnresolved','The endpoint peak is nonfinite during event localization.');
                end
                previousBest=abs(bestResidual);
                if abs(value)<previousBest
                    bestState=candidateState;bestRecords=candidateRecords;
                    bestH=probeH;bestResidual=value;
                end
                if abs(bestResidual)>.9*previousBest,slowProgress=slowProgress+1;
                else,slowProgress=0;end
                if value<0,lo=probeH;flo=value;
                else,hi=probeH;fhi=value;end
            end
            located=located||abs(bestResidual)<=xtol;
            if ~located
                error('ipm_cr:PeakEventUnresolved', ...
                    ['The original maximum did not continuously reach knot %.17g; ' ...
                    'best x residual %.3e. A nonunique peak or unresolved event cannot be projected.'], ...
                    crossing.knot,bestResidual);
            end
            s=bestState;r=bestRecords;h=bestH;
            event.located=true;
            event.rootResidual=crossing.direction*bestResidual;
            event.rootTimeBracket=[lo,hi];
            event.status='located_peak_cell_boundary';break;
        end
        event.limited=h<event.initialProposedDt;
        event.acceptedDt=h;event.probeCount=probes;
        event.extraRhsEvaluations=rhsEvaluations-callsBefore;
        event.stageGuardCuts=guardCuts;

        function [probedState,probedRecords]=probe(step)
            probes=probes+1;savedDt=dt;savedProbe=currentPeakProbe;
            dt=step;currentPeakProbe=probes;
            try
                [probedState,probedRecords]=attempt_rk4(initial,step);
            catch problem
                dt=savedDt;currentPeakProbe=savedProbe;rethrow(problem);
            end
            dt=savedDt;currentPeakProbe=savedProbe;
            if isempty(probedState)||~inside_cfl(max(probedRecords.rates),step)
                error('ipm_cr:ExplicitStageCfl', ...
                    'A complete peak-event RK4 probe exceeded the binding stage CFL.');
            end
        end
    end
    function crossing=find_crossing(r,h,x0)
        crossing=struct('found',false,'knot',NaN,'direction',0);
        firstTime=Inf;
        % Exclude only a local ulp band around the accepted starting knot.
        % A step leaving a located event then seeks the NEXT knot, avoiding
        % repeated zero-length steps without changing the measured peak.
        if numel(wallKnots)<3,return;end
        knots=wallKnots(2:end-1);
        tolerance=64*eps(max(1,abs(knots)));
        for j=2:numel(r.peaks)
            p=r.peaks(j);
            if ~p.available||~isfinite(p.x)||r.abscissae(j)<=0,continue;end
            if p.x>x0
                index=find(knots>x0+tolerance&knots<=p.x+tolerance,1,'first');direction=1;
            elseif p.x<x0
                index=find(knots<x0-tolerance&knots>=p.x-tolerance,1,'last');direction=-1;
            else,continue;end
            if isempty(index),continue;end
            estimate=h*r.abscissae(j)*abs((knots(index)-x0)/(p.x-x0));
            if estimate<firstTime
                crossing=struct('found',true,'knot',knots(index),'direction',direction);
                firstTime=estimate;
            end
        end
    end
    function index=note_peak(f,tau,abscissa,reused,returned)
        entry=peak_snapshot(f);
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

function r=record_start(f,clockRate)
r=struct('rates',f.rate,'positions',{{f.cflPosition}},'abscissae',0, ...
    'peaks',peak_snapshot(f), ...
    'gaugeResidualInfinity',max(abs(f.gaugeResidual)),'clockRates',clockRate,'physicalClockIncrement',NaN);
end
function r=record_add(r,f,clockRate,abscissa)
r.rates(end+1)=f.rate;r.positions{end+1}=f.cflPosition;
r.abscissae(end+1)=abscissa;r.gaugeResidualInfinity(end+1)=max(abs(f.gaugeResidual));
r.clockRates(end+1)=clockRate;
r.peaks(end+1)=peak_snapshot(f);
end
function p=peak_snapshot(f)
% Observe existing measurements only; missing diagnostics never reject a step.
p=struct('x',NaN,'cellIndex',NaN,'cellFraction',NaN, ...
    'tieNonunique',false,'curvature',NaN,'available',false);
if ~isstruct(f)||~isscalar(f)||~isfield(f,'measurement')|| ...
        ~isstruct(f.measurement)||~isfield(f.measurement,'peak')|| ...
        ~isstruct(f.measurement.peak)||~isscalar(f.measurement.peak)
    return;
end
source=f.measurement.peak;names={'x','cellIndex','cellFraction','tieNonunique','curvature'};
for k=1:numel(names)
    if isfield(source,names{k}),p.(names{k})=source.(names{k});end
end
p.available=true;
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
function r=add(a,wa,b,wb)
r=struct('core',wa*a.core+wb*b.core,'polar',wa*a.polar+wb*b.polar);
end
function z=pack(s)
z=[s.logC_l;s.logC_omega;s.physicalTime;s.X_shift;s.canonicalTime;s.densityOrigin];
end
function s=unpack(z,s)
s.logC_l=z(1);s.logC_omega=z(2);s.physicalTime=z(3);s.X_shift=z(4);s.canonicalTime=z(5);s.densityOrigin=z(6);
end
function k=scalar_rate(z,f)
k=[f.c_l;f.c_omega;exp(z(2)-z(1));f.c_l*z(4)+f.c_r;1;f.transportOrigin+f.c_omega*z(6)];
assert(all(isfinite(k))&&k(3)>0,'ipm_cr:ExplicitScale', ...
    'Every scalar rate must be finite, with a positive physical clock speed.');
end
function [high,low]=clock_add(high,low,increment)
% Preserve the historical compensated positive-increment clock arithmetic.
s=high+increment;b=s-high;e=(high-(s-b))+(increment-b);
t=low+e;high=s+t;low=(s-high)+t;
end
