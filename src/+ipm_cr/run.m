function [summary,state,ctx]=run(overrides)
%RUN Accepted-only, resumable explicit hybrid trajectory using the public API.
% Options: outputDir, tauEnd=200, wallSeconds=600, checkpointEveryTau=1,
% checkpointWallSeconds=60, logEveryTau=1, maxSteps=2e6, resume=true,
% verbose=true, rethrowOnFailure=true, and explicitOptions. The initialize
% settings geometry/solverOptions/transport/cfl/maxDt/minDt/timeIntegrator/background are also accepted
% directly. Only complete accepted states may replace the durable pointer.
if nargin<1,overrides=struct();end
root=ipm_cr.io.identity('root');
[opts,requested]=settings(overrides,root);folder=opts.outputDir;
if ~isfolder(folder),mkdir(folder);end
file=java.io.RandomAccessFile(fullfile(folder,'owner.lock'),'rw');channel=file.getChannel();lock=channel.tryLock();
assert(~isempty(lock),'ipm_cr:ExplicitRunnerOwner','This output folder already has a live owner.');
ownerCleanup=onCleanup(@()release_owner(lock,channel,file)); %#ok<NASGU>
logfid=fopen(fullfile(folder,'run.log'),'a');assert(logfid>=0,'ipm_cr:ExplicitRunnerLog');
logCleanup=onCleanup(@()fclose(logfid)); %#ok<NASGU>
threads=maxNumCompThreads(1);threadCleanup=onCleanup(@()maxNumCompThreads(threads)); %#ok<NASGU>
timer=tic;state=[];ctx=[];g=[];wallOps=[];method=[];pointer=[];priorWall=0;
auditHistory=struct([]);
historyColumns=ipm_cr.io.history('columns');history=zeros(0,numel(historyColumns));
manifest=struct([]);authority=struct();
summary=struct('kind','ipm_cr_run_v2','status','initializing', ...
    'outputDir',folder,'targetTau',opts.tauEnd,'startedUtc',ipm_cr.io.identity('utc'), ...
    'timeIntegrator','rk4','explicitCflLimitApplied',true,'threads',1, ...
    'matlabVersion',version,'computer',computer,'sourceInventory',manifest, ...
    'authority',authority,'trajectoryAccuracyCertified',false, ...
    'gaugeResidualHistoryScope','Maximum over all complete RK stages and accepted endpoint, with initial endpoint included.', ...
    'timingScope','MATLAB invocation initialization, accepted and retried complete explicit steps, checkpoint readback/SHA, history and output; external MATLAB startup is excluded.', ...
    'strongResidualScope','Native-array source reconstruction has no independent off-grid strong PDE audit unless the returned solver explicitly supplies one.');
logline('START targetTau=%.17g wallSeconds=%.3f folder=%s',opts.tauEnd,opts.wallSeconds,folder);
try
    manifest=ipm_cr.io.identity('source');authority=ipm_cr.io.identity('authority');
    summary.sourceInventory=manifest;summary.authority=authority;
    cursor=fullfile(folder,'latest.json');
    if isfile(cursor)
        assert(opts.resume,'ipm_cr:ExplicitRunnerResume','Existing accepted data require resume=true or a new output directory.');
        [checkpoint,pointer]=ipm_cr.io.checkpoint('load',folder);
        assert(isequaln(checkpoint.sourceInventory,manifest), ...
            'ipm_cr:ExplicitRunnerSourceChanged','A recorded source changed; use a new run for a changed implementation.');
        assert(isequaln(checkpoint.authority,authority), ...
            'ipm_cr:ExplicitRunnerAuthorityChanged','The original immutable datum authority changed.');
        state=checkpoint.state;g=checkpoint.geometry;wallOps=checkpoint.wallOps;
        method=checkpoint.explicitOptions;history=checkpoint.history;
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
        [state,ctx]=ipm_cr.init(requested);g=ctx.mesh;wallOps=ctx.wall;method=ctx.opts;
        history=ipm_cr.io.history('row',state,struct(),toc(timer));summary.resumed=false;
        [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
            manifest,authority,priorWall+toc(timer),pointer,auditHistory);
        logline('INITIAL step=0 tau=0 core=%dx%d annulus=%dx%d cfl=%.6g', ...
            g.core.ny,g.core.nx,g.polar.nr,g.polar.ntheta,method.cfl);
    end
    assert(opts.tauEnd>=state.scale.canonicalTime,'ipm_cr:ExplicitRunnerTarget', ...
        'The requested target precedes the accepted checkpoint.');
    startStep=state.step;startTau=state.scale.canonicalTime;
    summary.initializationSeconds=toc(timer);summary.startStep=startStep;summary.startTau=startTau;
    nextCheckpoint=(floor(state.scale.canonicalTime/opts.checkpointEveryTau)+1)*opts.checkpointEveryTau;
    nextLog=(floor(state.scale.canonicalTime/opts.logEveryTau)+1)*opts.logEveryTau;
    lastCheckpointWall=toc(timer);status='running';
    while state.scale.canonicalTime<opts.tauEnd&&state.step<opts.maxSteps
        if toc(timer)>=opts.wallSeconds,status='wall_budget';break;end
        if isfile(fullfile(folder,'STOP')),status='stopped';break;end
        [trial,info]=ipm_cr.step(state,ctx,opts.tauEnd);
        assert(info.explicitCflLimitApplied&&info.realizedMaximumStageCfl<=method.cfl*(1+64*eps), ...
            'ipm_cr:ExplicitRunnerCfl','A trial violated the configured complete all-stage explicit CFL.');
        state=trial;history(end+1,:)=ipm_cr.io.history('row',state,info,priorWall+toc(timer)); %#ok<AGROW>
        if state.scale.canonicalTime>=nextCheckpoint||toc(timer)-lastCheckpointWall>=opts.checkpointWallSeconds
            [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
                manifest,authority,priorWall+toc(timer),pointer,auditHistory);
            nextCheckpoint=(floor(state.scale.canonicalTime/opts.checkpointEveryTau)+1)*opts.checkpointEveryTau;
            lastCheckpointWall=toc(timer);
        end
        if state.scale.canonicalTime>=nextLog
            log_progress(state,info);
            nextLog=(floor(state.scale.canonicalTime/opts.logEveryTau)+1)*opts.logEveryTau;
        end
    end
    if strcmp(status,'running')
        if state.scale.canonicalTime>=opts.tauEnd,status='target_reached';else,status='maximum_steps';end
    end
    if isempty(pointer)||pointer.sequence~=state.step
        [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
            manifest,authority,priorWall+toc(timer),pointer,auditHistory);
    end
    summary=ipm_cr.io.history('summary',summary,status,state,g,method,history,historyColumns,pointer,priorWall,toc(timer));
    summary.checkpointAudits=auditHistory;
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

function [opts,requested]=settings(overrides,root)
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:ExplicitRunnerOptions');
opts=struct('outputDir',fullfile(root,'results','explicit_tau200'), ...
    'tauEnd',200,'wallSeconds',600,'checkpointEveryTau',1,'checkpointWallSeconds',60, ...
    'logEveryTau',1,'maxSteps',2000000,'resume',true,'verbose',true,'rethrowOnFailure',true);
methodFields={'geometry','solverOptions','transport','cfl','maxDt','minDt','timeIntegrator','background'};requested=struct();
if isfield(overrides,'explicitOptions'),requested=overrides.explicitOptions;end
assert(isstruct(requested)&&isscalar(requested),'ipm_cr:ExplicitRunnerOptions');
for alias={'mesh','poisson'}
    name=alias{1};target='geometry';if strcmp(name,'poisson'),target='solverOptions';end
    if isfield(requested,name)
        assert(~isfield(requested,target),'ipm_cr:ExplicitRunnerOptions','Duplicate alias %s.',target);
        requested.(target)=requested.(name);requested=rmfield(requested,name);
    end
end
assert(all(ismember(fieldnames(requested),methodFields)),'ipm_cr:ExplicitRunnerOptions', ...
    'explicitOptions contains an unknown numerical option.');
names=fieldnames(overrides);
for k=1:numel(names)
    name=names{k};
    if strcmp(name,'explicitOptions'),continue;
    elseif strcmp(name,'mesh'),requested.geometry=overrides.(name);
    elseif strcmp(name,'poisson'),requested.solverOptions=overrides.(name);
    elseif ismember(name,methodFields),requested.(name)=overrides.(name);
    else
        assert(isfield(opts,name),'ipm_cr:ExplicitRunnerOptions','Unknown runner option %s.',name);
        opts.(name)=overrides.(name);
    end
end
opts.outputDir=char(opts.outputDir);
assert(startsWith(opts.outputDir,filesep),'ipm_cr:ExplicitRunnerOutput','outputDir must be absolute.');
for name={'tauEnd','wallSeconds','checkpointEveryTau','checkpointWallSeconds','logEveryTau'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(opts.maxSteps,{'numeric'},{'scalar','integer','positive'});
for name={'resume','verbose','rethrowOnFailure'}
    validateattributes(opts.(name{1}),{'logical'},{'scalar'});
end
end

function [pointer,audits]=commit(folder,state,g,wall,method,opts,history,names,source,authority,seconds,pointer,audits)
cp=struct('state',state,'geometry',g,'wallOps',wall,'explicitOptions',method, ...
    'runnerOptions',opts,'history',history,'historyColumns',{names}, ...
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
