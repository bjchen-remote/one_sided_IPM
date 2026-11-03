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
auditHistory=struct([]);remeshHistory=struct([]);
historyColumns=ipm_cr.io.history('columns');history=zeros(0,numel(historyColumns));
manifest=struct([]);authority=struct();
summary=struct('kind','ipm_cr_run_v2','status','initializing', ...
    'outputDir',folder,'targetTau',opts.tauEnd,'startedUtc',ipm_cr.io.identity('utc'), ...
    'timeIntegrator','rk4','cflStepFormulaApplied',true,'threads',1, ...
    'runControls',opts,'steadyStateCertified',false, ...
    'matlabVersion',version,'computer',computer,'sourceInventory',manifest, ...
    'authority',authority,'trajectoryAccuracyCertified',false, ...
    'gaugeResidualHistoryScope','Maximum over all complete RK stages and accepted endpoint, with initial endpoint included.', ...
    'timingScope','MATLAB invocation initialization, direct complete explicit steps, checkpoint readback/SHA, history and output; external MATLAB startup is excluded.', ...
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
        if isfield(checkpoint,'remeshHistory'),remeshHistory=checkpoint.remeshHistory;end
        verify_mesh_history(folder,remeshHistory);
        savedSeed=[];
        if isfield(checkpoint.runnerOptions,'seed'),savedSeed=checkpoint.runnerOptions.seed;end
        if isempty(opts.seed),opts.seed=savedSeed;end
        if ~isempty(opts.seed)
            assert(~isempty(savedSeed)&&strcmp(opts.seed.sha256,savedSeed.sha256)&& ...
                opts.seed.tau==savedSeed.tau, ...
                'ipm_cr:SeedIdentity','Resume cannot replace the original seed identity or time.');
            % A frozen run can move between machines or directories. The
            % content and actual seed time are identity; its absolute path
            % is only the current location, checked again on every resume.
            assert(isfile(opts.seed.file)&&strcmp(ipm_cr.io.sha(opts.seed.file),opts.seed.sha256), ...
                'ipm_cr:SeedIdentity','The resumed seed file does not match its declared SHA256.');
        end
        ctx=ipm_cr.context(g,wallOps,method);
        auditHistory=checkpoint.auditHistory;
        assert(isequal(checkpoint.historyColumns,historyColumns), ...
            'ipm_cr:ExplicitRunnerHistory','Checkpoint history schema differs.');
        names=fieldnames(requested);
        for k=1:numel(names)
            assert(isequaln(resolved.(names{k}),method.(names{k})), ...
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
        if isempty(opts.seed),[state,ctx]=ipm_cr.init(resolved);
        else,[state,ctx]=from_seed(resolved,opts.seed);end
        g=ctx.mesh;wallOps=ctx.wall;method=ctx.opts;
        if isfield(state,'import'),summary.seedImport=state.import;end
        history=ipm_cr.io.history('row',state,struct(),toc(timer));summary.resumed=false;
        [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
            manifest,authority,priorWall+toc(timer),pointer,auditHistory,remeshHistory);
        logline('INITIAL step=%d tau=%.17g core=%dx%d annulus=%dx%d cfl=%.6g', ...
            state.step,state.scale.canonicalTime,g.core.ny,g.core.nx,g.polar.nr,g.polar.ntheta,method.cfl);
    end
    summary.authority=authority;
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
        if ~isempty(opts.stop.physicalTime)&& ...
                (state.scale.physicalTime-opts.stop.physicalTime)+state.scale.physicalTimeLow>=0
            status='physical_time';break;
        end
        [trial,info]=ipm_cr.step(state,ctx,opts.tauEnd);
        [state,ctx,meshInfo]=ipm_cr.remesh(trial,ctx,method.remesh);
        if ~strcmp(meshInfo.status,'not_triggered')
            entry=record_remesh(folder,g,ctx.mesh,meshInfo);
            if isempty(remeshHistory),remeshHistory=entry;else,remeshHistory(end+1)=entry;end %#ok<AGROW>
            logline('REMESH step=%d tau=%.12g status=%s reason=%s', ...
                state.step,state.scale.canonicalTime,meshInfo.status,meshInfo.reason);
        end
        g=ctx.mesh;method=ctx.opts;
        history(end+1,:)=ipm_cr.io.history('row',state,info,priorWall+toc(timer)); %#ok<AGROW>
        if state.scale.canonicalTime>=nextCheckpoint||toc(timer)-lastCheckpointWall>=opts.checkpointWallSeconds
            [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns, ...
                manifest,authority,priorWall+toc(timer),pointer,auditHistory,remeshHistory);
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
            manifest,authority,priorWall+toc(timer),pointer,auditHistory,remeshHistory);
    end
    summary=ipm_cr.io.history('summary',summary,status,state,g,method,history,historyColumns,pointer,priorWall,toc(timer));
    summary.checkpointAudits=auditHistory;
    summary.remeshHistory=remeshHistory;
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
            'stop',opts.stop);
        summary.dataPackage=ipm_cr.pack(state,ctx,folder,meta,packageHistory,historyColumns);
    end
    summary.invocationWallSeconds=toc(timer);summary.accumulatedWallSeconds=priorWall+toc(timer);
    summary.endedUtc=ipm_cr.io.identity('utc');
    ipm_cr.io.history('write',folder,summary,history,historyColumns);
    logline('END status=%s tau=%.17g steps=%d wall=%.6f maxStageCfl=%.6g gauge=%.3e', ...
        status,state.scale.canonicalTime,state.step,summary.invocationWallSeconds, ...
        summary.maximumStageCfl,summary.maximumRecordedGaugeResidual);
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
        logline('ACCEPT step=%d tau=%.12g dt=%.6g maxStageCFL=%.6g gauge=%.3e backward=%.3e strong=%g Q=%.3e wall=%.3fs', ...
            s.step,s.scale.canonicalTime,info.dt,info.realizedMaximumStageCfl, ...
            max(abs(s.rhs.gaugeResidual)),get(p,'linearSolveBackwardError',NaN), ...
            get(p,'strongResidual',NaN),max([max(abs(s.rhs.Q.core),[],'all'),max(abs(s.rhs.Q.polar),[],'all')]),toc(timer));
    end
end

function [state,ctx]=from_seed(method,seed)
% An explicit cross-method/grid import is a NEW run, never a fake resume.
assert(isfile(seed.file)&&strcmp(ipm_cr.io.sha(seed.file),seed.sha256), ...
    'ipm_cr:SeedIdentity','Seed file does not match its declared SHA256.');
loaded=load(seed.file,'checkpoint');
assert(strcmp(ipm_cr.io.sha(seed.file),seed.sha256),'ipm_cr:SeedIdentity', ...
    'Seed changed during read.');
assert(isfield(loaded,'checkpoint'),'ipm_cr:Seed','A complete accepted checkpoint is required.');
parent=loaded.checkpoint;
assert(all(isfield(parent,{'state','geometry','wallOps','acceptedOnly'}))&&parent.acceptedOnly&& ...
    isfield(parent.state,'W')&&parent.state.scale.canonicalTime==seed.tau, ...
    'ipm_cr:Seed','Seed is not an accepted state at the declared actual tau.');
[~,cfg]=ipm_cr.config(method);[wall,~,reference]=ipm_cr.model.reference(cfg);
assert(isequaln(wall,parent.wallOps),'ipm_cr:SeedAuthority','Seed original wall differs.');
assert(isfield(parent,'explicitOptions')&& ...
    isequaln(parent.explicitOptions.background,method.background), ...
    'ipm_cr:SeedBackground','Import cannot reinterpret the residual under a different background.');
assert(isfield(parent,'authority')&&isfield(parent.authority,'sha256')&& ...
    strcmp(parent.authority.sha256,cfg.authoritySha), ...
    'ipm_cr:SeedAuthority','Seed initial authority differs.');
g=ipm_cr.mesh.build(method.geometry);ctx=ipm_cr.context(g,wall,method);
[state,transfer]=ipm_cr.remesh('transfer',parent.state,parent.geometry,g);
if isfield(state,'remesh'),state=rmfield(state,'remesh');end
state.ctxId=ctx.id;state.ctxKey=ctx.key;
[state.rhs,state.W]=ipm_cr.rhs(state,ctx);
comparison=seed_comparison(parent.state,parent.geometry,state,g,method);
state.import=struct('kind','explicit_cross_method_or_grid_seed_v1', ...
    'file',seed.file,'sha256',seed.sha256,'tau',seed.tau,'parentStep',state.step, ...
    'transfer',transfer,'comparison',comparison, ...
    'sameMethodContinuationClaimed',false,'timeRelabeled',false, ...
    'originalGeometryKind',parent.geometry.kind,'initialGeometryKind',g.kind, ...
    'gaugeResidual',state.rhs.gaugeResidual,'poisson',state.rhs.poisson);
end

function result=seed_comparison(a,ga,b,gb,method)
[x,y]=meshgrid(linspace(gb.core.activeXBounds(1),gb.core.activeXBounds(2),17), ...
    linspace(0,gb.core.activeYBounds(2),17));
r=exp(linspace(log(gb.polar.innerActiveRadius),log(gb.polar.r(end)),17))';r(end)=gb.polar.r(end);
theta=linspace(0,pi,17);px=r*cos(theta);py=r*sin(theta);py(:,[1,end])=0;
X=[x(:);px(:)];Y=[y(:);py(:)];result=struct('queryCount',numel(X));
% Honor the same analytic-background representation as the public sampler.
% Interpolating complete W would instead measure coarse interpolation of H,
% even when the actual transferred residual is exactly unchanged.
av=ipm_cr.diag.fields(a,ga,method,X,Y);bv=ipm_cr.diag.fields(b,gb,method,X,Y);
for name={'W','q','F','Q'}
    key=name{1};aa=av.(key);bb=bv.(key);
    result.(key)=max(abs(aa-bb)./max(1,max(abs(aa),abs(bb))));
end
result.scope='saved parent and fresh import using analytic background plus sampled residual at common physical points; not continuum truth';
end

function [pointer,audits]=commit(folder,state,g,wall,method,opts,history,names,source,authority,seconds,pointer,audits,remeshHistory)
cp=struct('state',state,'geometry',g,'wallOps',wall,'explicitOptions',method, ...
    'runnerOptions',opts,'history',history,'historyColumns',{names}, ...
    'sourceInventory',source,'authority',authority,'auditHistory',audits,'remeshHistory',remeshHistory, ...
    'accumulatedWallSeconds',seconds);
[pointer,cp]=ipm_cr.io.checkpoint('save',folder,cp,pointer);audits=cp.auditHistory;
end

function verify_mesh_history(folder,history)
for k=1:numel(history)
    entry=history(k);
    assert(~isempty(regexp(entry.file,'^mesh/[A-Za-z0-9_]+\.mat$','once')), ...
        'ipm_cr:RemeshStorage','Mesh lineage must reference a local mesh event basename.');
    file=fullfile(folder,strrep(entry.file,'/',filesep));
    assert(isfile(file)&&strcmp(ipm_cr.io.sha(file),entry.sha256), ...
        'ipm_cr:RemeshStorage','Missing or changed mesh lineage file: %s.',entry.file);
end
end

function entry=record_remesh(folder,before,after,info)
% The accepted cursor owns this lineage. Orphan events after a crash cannot
% overwrite an earlier event, and are not mistaken for committed geometry.
directory=fullfile(folder,'mesh');if ~isfolder(directory),mkdir(directory);end
temporary=[tempname(directory),'.mat'];
event=struct('format','ipm_cr_remesh_event_v1','info',info, ...
    'before',before.geometryKey,'after',after.geometryKey);
save(temporary,'event','-v7.3');readback=load(temporary,'event');
assert(isequaln(event,readback.event),'ipm_cr:RemeshStorage','Mesh event readback differs.');
entry=struct('step',info.step,'tau',info.tau,'status',info.status, ...
    'reason',info.reason,'file','','sha256',ipm_cr.io.sha(temporary));
[~,name,extension]=fileparts(temporary);entry.file=['mesh/',name,extension];
end

function value=get(s,name,default)
if isfield(s,name),value=s.(name);else,value=default;end
end

function release_owner(lock,channel,file)
lock.release();channel.close();file.close();
end
