function [summary,state,ctx]=run(overrides)
%RUN Accepted-only, resumable explicit hybrid trajectory using the public API.
% Use config('case') or a JSON case filename. Legacy flat inputs are resolved
% only by config. Stops occur at accepted endpoints; tauEnd is an exact cap,
% while physical-time thresholds report any accepted-step overshoot.
if nargin<1,overrides=struct();end
if ischar(overrides)||isstring(overrides),overrides=jsondecode(fileread(overrides));end
root=ipm_cr.io.identity('root');
[opts,resolved,requested]=ipm_cr.config(overrides,'run');folder=opts.outputDir;
if ~isfolder(folder),mkdir(folder);end
file=java.io.RandomAccessFile(fullfile(folder,'owner.lock'),'rw');channel=file.getChannel();lock=channel.tryLock();
assert(~isempty(lock),'ipm_cr:ExplicitRunnerOwner','This output folder already has a live owner.');
ownerCleanup=onCleanup(@()release_owner(lock,channel,file)); %#ok<NASGU>
logfid=fopen(fullfile(folder,'run.log'),'a');assert(logfid>=0,'ipm_cr:ExplicitRunnerLog');
logCleanup=onCleanup(@()fclose(logfid)); %#ok<NASGU>
threads=maxNumCompThreads(1);threadCleanup=onCleanup(@()maxNumCompThreads(threads)); %#ok<NASGU>
timer=tic;state=[];ctx=[];g=[];wallOps=[];method=[];pointer=[];priorWall=0;
auditHistory=struct([]);stopState=[];stopRecords=struct([]);stopDirty=false;
historyColumns=ipm_cr.io.history('columns');history=zeros(0,numel(historyColumns));
manifest=struct([]);authority=struct();
summary=struct('kind','ipm_cr_run_v2','status','initializing', ...
    'outputDir',folder,'targetTau',opts.tauEnd,'startedUtc',ipm_cr.io.identity('utc'), ...
    'timeIntegrator','rk4','explicitCflLimitApplied',true,'threads',1, ...
    'runControls',opts,'steadyStateCertified',false,'steadyCandidate',false, ...
    'matlabVersion',version,'computer',computer,'sourceInventory',manifest, ...
    'authority',authority,'trajectoryAccuracyCertified',false, ...
    'gaugeResidualHistoryScope','Maximum over all complete RK stages and accepted endpoint, with initial endpoint included.', ...
    'timingScope','MATLAB invocation initialization, accepted and retried complete explicit steps, checkpoint readback/SHA, history and output; external MATLAB startup is excluded.', ...
    'strongResidualScope','Native-array source reconstruction has no independent off-grid strong PDE audit unless the returned solver explicitly supplies one.');
logline('START targetTau=%.17g wallSeconds=%.3f folder=%s',opts.tauEnd,opts.wallSeconds,folder);
try
    manifest=ipm_cr.io.identity('source');summary.sourceInventory=manifest;
    cursor=fullfile(folder,'latest.json');
    if isfile(cursor)
        assert(opts.resume,'ipm_cr:ExplicitRunnerResume','Existing accepted data require resume=true or a new output directory.');
        [checkpoint,pointer]=ipm_cr.io.checkpoint('load',folder);
        assert(isequaln(checkpoint.sourceInventory,manifest), ...
            'ipm_cr:ExplicitRunnerSourceChanged','A recorded source changed; use a new run for a changed implementation.');
        method=checkpoint.explicitOptions;
        authority=ipm_cr.io.identity('authority',root,method.initial);
        assert(isequaln(checkpoint.authority,authority), ...
            'ipm_cr:ExplicitRunnerAuthorityChanged','The original immutable datum authority changed.');
        state=checkpoint.state;g=checkpoint.geometry;wallOps=checkpoint.wallOps;
        history=checkpoint.history;
        if isfield(checkpoint,'stopState')&&isfield(checkpoint,'runnerOptions')&& ...
                isfield(checkpoint.runnerOptions,'stop')&&isequaln(checkpoint.runnerOptions.stop,opts.stop)
            stopState=checkpoint.stopState;
            if isfield(checkpoint,'stopRecords'),stopRecords=checkpoint.stopRecords;end
        elseif isfield(checkpoint,'stopState')&&~isempty(checkpoint.stopState)
            summary.stopPolicyReset=true;
        end
        ctx=ipm_cr.context(g,wallOps,method);
        auditHistory=checkpoint.auditHistory;
        assert(isequal(checkpoint.historyColumns,historyColumns), ...
            'ipm_cr:ExplicitRunnerHistory','Checkpoint history schema differs.');
        names=fieldnames(requested);
        for k=1:numel(names)
            assert(isequaln(requested.(names{k}),method.(names{k})), ...
                'ipm_cr:ExplicitRunnerOptionsChanged','Resume cannot change numerical option %s.',names{k});
        end
        priorWall=checkpoint.accumulatedWallSeconds;
        if isfield(pointer,'accumulatedWallSeconds'),priorWall=max(priorWall,pointer.accumulatedWallSeconds);end
        previousSummary=fullfile(folder,'summary.json');
        if isfile(previousSummary)
            previous=jsondecode(fileread(previousSummary));
            if isfield(previous,'acceptedCheckpoint')&&isstruct(previous.acceptedCheckpoint)&& ...
                isfield(previous.acceptedCheckpoint,'sha256')&&strcmp(previous.acceptedCheckpoint.sha256,pointer.sha256)&& ...
                isfield(previous,'accumulatedWallSeconds')
                priorWall=max(priorWall,previous.accumulatedWallSeconds);
            end
        end
        oldF=state.rhs.F;
        state=ipm_cr.refresh(state,ctx);
        summary.resumeRefreshCenteredFInfinity=max([max(abs(oldF.core-state.rhs.F.core),[],'all'), ...
            max(abs(oldF.polar-state.rhs.F.polar),[],'all')]);
        summary.resumed=true;
        logline('RESUME step=%d tau=%.17g checkpoint=%s',state.step,state.scale.canonicalTime,pointer.sha256);
    else
        authority=ipm_cr.io.identity('authority',root,resolved.initial);
        [state,ctx]=ipm_cr.init(resolved);g=ctx.mesh;wallOps=ctx.wall;method=ctx.opts;
        history=ipm_cr.io.history('row',state,struct(),toc(timer));summary.resumed=false;
        [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
            manifest,authority,priorWall+toc(timer),pointer,auditHistory,stopState,stopRecords);
        logline('INITIAL step=0 tau=0 core=%dx%d annulus=%dx%d cfl=%.6g', ...
            g.core.ny,g.core.nx,g.polar.nr,g.polar.ntheta,method.cfl);
    end
    summary.authority=authority;
    if opts.stop.steady.enabled,observe_steady();end
    assert(opts.tauEnd>=state.scale.canonicalTime,'ipm_cr:ExplicitRunnerTarget', ...
        'The requested target precedes the accepted checkpoint.');
    startStep=state.step;startTau=state.scale.canonicalTime;
    summary.initializationSeconds=toc(timer);summary.startStep=startStep;summary.startTau=startTau;
    nextCheckpoint=(floor(state.scale.canonicalTime/opts.checkpointEveryTau)+1)*opts.checkpointEveryTau;
    nextLog=(floor(state.scale.canonicalTime/opts.logEveryTau)+1)*opts.logEveryTau;
    lastCheckpointWall=toc(timer);status='running';
    nextObservation=(floor(state.scale.canonicalTime/opts.stop.steady.sampleEveryTau)+1)*opts.stop.steady.sampleEveryTau;
    while state.scale.canonicalTime<opts.tauEnd&&state.step<opts.maxSteps
        if toc(timer)>=opts.wallSeconds,status='wall_budget';break;end
        if isfile(fullfile(folder,'STOP')),status='stopped';break;end
        if ~isempty(opts.stop.physicalTime)&& ...
                (state.scale.physicalTime-opts.stop.physicalTime)+state.scale.physicalTimeLow>=0
            status='physical_time';break;
        end
        if summary.steadyCandidate,status='steady_candidate';break;end
        [trial,info]=ipm_cr.step(state,ctx,opts.tauEnd);
        assert(info.explicitCflLimitApplied&&info.realizedMaximumStageCfl<=method.cfl*(1+64*eps), ...
            'ipm_cr:ExplicitRunnerCfl','A trial violated the configured complete all-stage explicit CFL.');
        state=trial;history(end+1,:)=ipm_cr.io.history('row',state,info,priorWall+toc(timer)); %#ok<AGROW>
        if opts.stop.steady.enabled&&state.scale.canonicalTime>=nextObservation
            observe_steady();
            nextObservation=(floor(state.scale.canonicalTime/opts.stop.steady.sampleEveryTau)+1)*opts.stop.steady.sampleEveryTau;
        end
        if state.scale.canonicalTime>=nextCheckpoint||toc(timer)-lastCheckpointWall>=opts.checkpointWallSeconds
            [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
                manifest,authority,priorWall+toc(timer),pointer,auditHistory,stopState,stopRecords);
            nextCheckpoint=(floor(state.scale.canonicalTime/opts.checkpointEveryTau)+1)*opts.checkpointEveryTau;
            lastCheckpointWall=toc(timer);stopDirty=false;
        end
        if state.scale.canonicalTime>=nextLog
            log_progress(state,info);
            nextLog=(floor(state.scale.canonicalTime/opts.logEveryTau)+1)*opts.logEveryTau;
        end
    end
    if opts.stop.steady.enabled,observe_steady();end
    if strcmp(status,'running')
        if summary.steadyCandidate,status='steady_candidate';
        elseif state.scale.canonicalTime>=opts.tauEnd,status='target_reached';else,status='maximum_steps';end
    end
    if isempty(pointer)||pointer.sequence~=state.step||stopDirty
        [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
            manifest,authority,priorWall+toc(timer),pointer,auditHistory,stopState,stopRecords);
    end
    summary=ipm_cr.io.history('summary',summary,status,state,g,method,history,historyColumns,pointer,priorWall,toc(timer));
    summary.checkpointAudits=auditHistory;summary.steadyObservations=stopRecords;
    if ~isempty(opts.stop.physicalTime)
        summary.physicalTimeStopOvershoot=(state.scale.physicalTime-opts.stop.physicalTime)+state.scale.physicalTimeLow;
    end
    ipm_cr.io.json(fullfile(folder,'case.json'),struct('format','ipm_cr_resolved_case_v1', ...
        'method',method,'run',opts,'authority',authority,'source',manifest));
    if opts.save.history
        writetable(array2table(history,'VariableNames',historyColumns),fullfile(folder,'history.csv'));
    end
    if opts.save.data
        if opts.save.history,packageHistory=history;else,packageHistory=[];end
        meta=struct('status',status,'label',opts.save.label,'acceptedCheckpoint',pointer, ...
            'stop',opts.stop,'steadyCandidate',summary.steadyCandidate);
        summary.dataPackage=ipm_cr.pack(state,ctx,folder,meta,packageHistory,historyColumns);
    end
    summary.invocationWallSeconds=toc(timer);summary.accumulatedWallSeconds=priorWall+toc(timer);
    summary.endedUtc=ipm_cr.io.identity('utc');
    ipm_cr.io.history('write',folder,summary,history,historyColumns);
    logline('END status=%s tau=%.17g steps=%d wall=%.6f maxStageCfl=%.6g gauge=%.3e', ...
        status,state.scale.canonicalTime,state.step,summary.invocationWallSeconds, ...
        summary.maximumAcceptedStageCfl,summary.maximumRecordedGaugeResidual);
catch problem
    % Preserve the durable accepted pointer/slots exactly as published.
    % The recovery payload is evidence only; it is never made the cursor.
    errorRecord=struct('identifier',problem.identifier,'message',problem.message, ...
        'stack',problem.stack,'utc',ipm_cr.io.identity('utc'),'invocationWallSeconds',toc(timer));
    summary.error=errorRecord;
    summary.checkpointAudits=auditHistory;
    if ~isempty(state)&&~isempty(method)
        summary=ipm_cr.io.history('summary',summary,'failed',state,g,method,history,historyColumns,pointer,priorWall,toc(timer));
    else
        summary.status='failed';summary.invocationWallSeconds=toc(timer);summary.endedUtc=ipm_cr.io.identity('utc');
    end
    lastAcceptedInMemory=state; %#ok<NASGU>
    errorPath=fullfile(folder,['error_',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')),'.mat']);
    save(errorPath,'errorRecord','lastAcceptedInMemory','-v7.3');
    summary.errorEvidence=errorPath;ipm_cr.io.history('write',folder,summary,history,historyColumns);
    logline('FAILED %s: %s',problem.identifier,problem.message);
    if opts.rethrowOnFailure,rethrow(problem);end
end

    function observe_steady()
        if ~isempty(stopState)&&isfield(stopState,'lastStep')&&stopState.lastStep==state.step
            summary.steadyCandidate=stopState.lastRecord.candidate;return;
        end
        policy=rmfield(opts.stop.steady,{'enabled','sampleEveryTau'});
        [stopState,record]=ipm_cr.diag.steady(state,ctx,policy,stopState);
        if isempty(stopRecords),stopRecords=record;else,stopRecords(end+1)=record;end
        stopDirty=true;summary.steadyCandidate=record.candidate;
        ipm_cr.io.json(fullfile(folder,'steady.json'),stopRecords);
    end
    function logline(format,varargin)
        line=sprintf(format,varargin{:});fprintf(logfid,'[%s] %s\n',ipm_cr.io.identity('utc'),line);
        if opts.verbose,fprintf('%s\n',line);end
    end
    function log_progress(s,info)
        p=s.rhs.poisson;
        logline('ACCEPT step=%d tau=%.12g dt=%.6g maxStageCFL=%.6g retry=%d gauge=%.3e backward=%.3e strong=%g Q=%.3e wall=%.3fs', ...
            s.step,s.scale.canonicalTime,info.dt,info.realizedMaximumStageCfl,info.retries, ...
            max(abs(s.rhs.gaugeResidual)),get(p,'linearSolveBackwardError',NaN), ...
            get(p,'strongResidual',NaN),max([max(abs(s.rhs.Q.core),[],'all'),max(abs(s.rhs.Q.polar),[],'all')]),toc(timer));
    end
end

function [pointer,audits]=commit(folder,state,g,wall,method,opts,history,names,source,authority,seconds,pointer,audits,stopState,stopRecords)
cp=struct('state',state,'geometry',g,'wallOps',wall,'explicitOptions',method, ...
    'runnerOptions',opts,'stopState',stopState,'stopRecords',stopRecords,'history',history,'historyColumns',{names}, ...
    'sourceInventory',source,'authority',authority,'auditHistory',audits, ...
    'accumulatedWallSeconds',seconds);
[pointer,cp]=ipm_cr.io.checkpoint('save',folder,cp,pointer);audits=cp.auditHistory;
end

function value=get(s,name,default)
if isfield(s,name),value=s.(name);else,value=default;end
end

function release_owner(lock,channel,file)
lock.release();channel.close();file.close();
end
