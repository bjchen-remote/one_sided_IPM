function [summary,state,g,wallOps]=run_ipm_cr_explicit(overrides)
%RUN_IPM_CR_EXPLICIT Resumable native hybrid fifth-order/RK4 experiment.
% Options: outputDir, tauEnd=200, wallSeconds=600, checkpointEveryTau=1,
% checkpointWallSeconds=60, logEveryTau=1, maxSteps=2e6, resume=true,
% verbose=true, rethrowOnFailure=true, and explicitOptions. The initialize
% settings geometry/solverOptions/transport/cfl/maxDt/minDt/timeIntegrator/background are also accepted
% directly. Only complete accepted states may replace the durable pointer.
if nargin<1,overrides=struct();end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'time_evolution'));
[opts,requested]=settings(overrides,root);folder=opts.outputDir;
if ~isfolder(folder),mkdir(folder);end
file=java.io.RandomAccessFile(fullfile(folder,'owner.lock'),'rw');channel=file.getChannel();lock=channel.tryLock();
assert(~isempty(lock),'ipm_cr:ExplicitRunnerOwner','This output folder already has a live owner.');
ownerCleanup=onCleanup(@()release_owner(lock,channel,file)); %#ok<NASGU>
logfid=fopen(fullfile(folder,'run.log'),'a');assert(logfid>=0,'ipm_cr:ExplicitRunnerLog');
logCleanup=onCleanup(@()fclose(logfid)); %#ok<NASGU>
threads=maxNumCompThreads(1);threadCleanup=onCleanup(@()maxNumCompThreads(threads)); %#ok<NASGU>
timer=tic;state=[];g=[];wallOps=[];method=[];pointer=[];priorWall=0;
auditHistory=struct([]);
historyColumns=columns();history=zeros(0,numel(historyColumns));
manifest=struct([]);authority=struct();
summary=struct('kind','ipm_cr_uniform_core_annular_explicit_run_v1','status','initializing', ...
    'outputDir',folder,'targetTau',opts.tauEnd,'startedUtc',utc_now(), ...
    'timeIntegrator','rk4','explicitCflLimitApplied',true,'threads',1, ...
    'matlabVersion',version,'computer',computer,'sourceInventory',manifest, ...
    'authority',authority,'trajectoryAccuracyCertified',false, ...
    'gaugeResidualHistoryScope','Maximum over all complete RK stages and accepted endpoint, with initial endpoint included.', ...
    'timingScope','MATLAB invocation initialization, accepted and retried complete explicit steps, checkpoint readback/SHA, history and output; external MATLAB startup is excluded.', ...
    'strongResidualScope','Native-array source reconstruction has no independent off-grid strong PDE audit unless the returned solver explicitly supplies one.');
logline('START targetTau=%.17g wallSeconds=%.3f folder=%s',opts.tauEnd,opts.wallSeconds,folder);
try
    manifest=source_inventory(root);authority=authority_identity();
    summary.sourceInventory=manifest;summary.authority=authority;
    cursor=fullfile(folder,'latest.json');
    if isfile(cursor)
        assert(opts.resume,'ipm_cr:ExplicitRunnerResume','Existing accepted data require resume=true or a new output directory.');
        [checkpoint,pointer]=load_accepted(folder);
        assert(isequaln(checkpoint.sourceInventory,manifest), ...
            'ipm_cr:ExplicitRunnerSourceChanged','A recorded source changed; use a new run for a changed implementation.');
        assert(isequaln(checkpoint.authority,authority), ...
            'ipm_cr:ExplicitRunnerAuthorityChanged','The original immutable datum authority changed.');
        state=checkpoint.state;g=checkpoint.geometry;wallOps=checkpoint.wallOps;
        method=checkpoint.explicitOptions;history=checkpoint.history;
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
        oldF=state.flow.F;
        if isfield(state,'residual')
            [state.flow,state.rho]=ipm_cr.explicit.evaluateResidual(state.residual,state.scale,g,wallOps,method,[state.flow.c_l,state.flow.c_r]);
        else
            state.flow=ipm_cr.explicit.evaluate(state.rho,g,wallOps,method,[state.flow.c_l,state.flow.c_r]);
        end
        summary.resumeRefreshCenteredFInfinity=max([max(abs(oldF.core-state.flow.F.core),[],'all'), ...
            max(abs(oldF.polar-state.flow.F.polar),[],'all')]);
        summary.resumed=true;
        logline('RESUME step=%d tau=%.17g checkpoint=%s',state.step,state.scale.canonicalTime,pointer.sha256);
    else
        [state,g,wallOps,method]=ipm_cr.explicit.initialize(requested);
        history=history_row(state,struct(),toc(timer));summary.resumed=false;
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
        [trial,info]=ipm_cr.explicit.advance(state,g,wallOps,method,opts.tauEnd);
        assert(info.explicitCflLimitApplied&&info.realizedMaximumStageCfl<=method.cfl*(1+64*eps), ...
            'ipm_cr:ExplicitRunnerCfl','A trial violated the configured complete all-stage explicit CFL.');
        state=trial;history(end+1,:)=history_row(state,info,priorWall+toc(timer)); %#ok<AGROW>
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
    summary=finish_summary(summary,status,state,g,method,history,historyColumns,pointer,priorWall,toc(timer));
    summary.checkpointAudits=auditHistory;
    write_summary(folder,summary,history,historyColumns);
    logline('END status=%s tau=%.17g steps=%d wall=%.6f maxStageCfl=%.6g gauge=%.3e', ...
        status,state.scale.canonicalTime,state.step,summary.invocationWallSeconds, ...
        summary.maximumAcceptedStageCfl,summary.maximumRecordedGaugeResidual);
catch problem
    % Preserve the durable accepted pointer/slots exactly as published.
    % The recovery payload is evidence only; it is never made the cursor.
    errorRecord=struct('identifier',problem.identifier,'message',problem.message, ...
        'stack',problem.stack,'utc',utc_now(),'invocationWallSeconds',toc(timer));
    summary.error=errorRecord;
    summary.checkpointAudits=auditHistory;
    if ~isempty(state)&&~isempty(method)
        summary=finish_summary(summary,'failed',state,g,method,history,historyColumns,pointer,priorWall,toc(timer));
    else
        summary.status='failed';summary.invocationWallSeconds=toc(timer);summary.endedUtc=utc_now();
    end
    lastAcceptedInMemory=state; %#ok<NASGU>
    errorPath=fullfile(folder,['error_',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')),'.mat']);
    save(errorPath,'errorRecord','lastAcceptedInMemory','-v7.3');
    summary.errorEvidence=errorPath;write_summary(folder,summary,history,historyColumns);
    logline('FAILED %s: %s',problem.identifier,problem.message);
    if opts.rethrowOnFailure,rethrow(problem);end
end

    function logline(format,varargin)
        line=sprintf(format,varargin{:});fprintf(logfid,'[%s] %s\n',utc_now(),line);
        if opts.verbose,fprintf('%s\n',line);end
    end
    function log_progress(s,info)
        p=s.flow.poisson;
        logline('ACCEPT step=%d tau=%.12g dt=%.6g maxStageCFL=%.6g retry=%d gauge=%.3e backward=%.3e strong=%g Q=%.3e wall=%.3fs', ...
            s.step,s.scale.canonicalTime,info.dt,info.realizedMaximumStageCfl,info.retries, ...
            max(abs(s.flow.gaugeResidual)),get(p,'linearSolveBackwardError',NaN), ...
            get(p,'strongResidual',NaN),max([max(abs(s.flow.Q.core),[],'all'),max(abs(s.flow.Q.polar),[],'all')]),toc(timer));
    end
end

function [opts,requested]=settings(overrides,root)
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:ExplicitRunnerOptions');
opts=struct('outputDir',fullfile(root,'test','results','explicit_tau200'), ...
    'tauEnd',200,'wallSeconds',600,'checkpointEveryTau',1,'checkpointWallSeconds',60, ...
    'logEveryTau',1,'maxSteps',2000000,'resume',true,'verbose',true,'rethrowOnFailure',true);
methodFields={'geometry','solverOptions','transport','cfl','maxDt','minDt','timeIntegrator','background'};requested=struct();
if isfield(overrides,'explicitOptions'),requested=overrides.explicitOptions;end
assert(isstruct(requested)&&isscalar(requested),'ipm_cr:ExplicitRunnerOptions');
assert(all(ismember(fieldnames(requested),methodFields)),'ipm_cr:ExplicitRunnerOptions', ...
    'explicitOptions contains an unknown numerical option.');
names=fieldnames(overrides);
for k=1:numel(names)
    name=names{k};
    if strcmp(name,'explicitOptions'),continue;
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

function names=columns()
names={'step','tau','physicalTime','physicalTimeLow','dt','maximumStageRate','maximumStageCfl', ...
    'configuredCfl','retries','stepSeconds','accumulatedWallSeconds','gaugeResidualInfinity', ...
    'poissonBackwardError','poissonStrongResidual','poissonStrongAvailable','overlapValueMismatch', ...
    'overlapVelocityMismatch','normalFluxMismatch','traceResidual','cl','cw','cr', ...
    'wallI','wallP','wallC','centeredDensityMin','centeredDensityMax','densityOrigin', ...
    'totalDensityMin','totalDensityMax','QInfinity','centeredFInfinity', ...
    'rhsEvaluations','peakEventLimited','peakEventLocated','peakEventProbeCount','peakEventExtraRhs', ...
    'peakTieStages','minimumPeakCurvatureMagnitude'};
end

function row=history_row(s,info,wall)
p=s.flow.poisson;event=get(info,'peakEvent',struct());
peakTies=0;peakCurvature=NaN;
if isfield(info,'stagePeakDiagnostics')&&~isempty(info.stagePeakDiagnostics)
    peaks=info.stagePeakDiagnostics;peakTies=nnz([peaks.tieNonunique]);
    peakCurvature=min(abs([peaks.curvature]));
end
centeredMin=min([min(s.rho.core,[],'all'),min(s.rho.polar,[],'all')]);
centeredMax=max([max(s.rho.core,[],'all'),max(s.rho.polar,[],'all')]);
row=[s.step,s.scale.canonicalTime,s.scale.physicalTime,s.scale.physicalTimeLow, ...
    get(info,'dt',NaN),get(info,'maximumStageRate',NaN),get(info,'realizedMaximumStageCfl',NaN), ...
    get(info,'configuredCfl',NaN),get(info,'retries',0),get(info,'seconds',0),wall, ...
    max(max(abs(s.flow.gaugeResidual)),get(info,'maximumStageGaugeResidual',0)),get(p,'linearSolveBackwardError',NaN), ...
    get(p,'strongResidual',NaN),get(p,'strongResidualAvailable',false), ...
    get(p,'overlapValueMismatch',NaN),get(p,'overlapVelocityMismatch',NaN), ...
    get(p,'normalFluxMismatch',NaN),get(p,'traceResidual',NaN), ...
    s.flow.c_l,s.flow.c_omega,s.flow.c_r,s.flow.measurement.I,s.flow.measurement.P,s.flow.measurement.C, ...
    centeredMin,centeredMax,s.scale.densityOrigin,centeredMin+s.scale.densityOrigin, ...
    centeredMax+s.scale.densityOrigin, ...
    max([max(abs(s.flow.Q.core),[],'all'),max(abs(s.flow.Q.polar),[],'all')]), ...
    max([max(abs(s.flow.F.core),[],'all'),max(abs(s.flow.F.polar),[],'all')]), ...
    get(info,'rhsEvaluations',0),get(event,'limited',false),get(event,'located',false), ...
    get(event,'probeCount',0),get(event,'extraRhsEvaluations',0),peakTies,peakCurvature];
end

function [pointer,auditHistory]=commit(folder,state,g,wallOps,method,opts,history,historyColumns,manifest,authority,wall,pointer,auditHistory)
commitTimer=tic;
spatialAudit=ipm_cr.explicit.audit(state,g);
entry=struct('step',state.step,'tau',state.scale.canonicalTime,'audit',spatialAudit);
if isempty(auditHistory),auditHistory=entry;else,auditHistory(end+1)=entry;end
slot=0;generation=0;
if ~isempty(pointer)
    slot=double(strcmp(pointer.path,fullfile(folder,'accepted_0.mat')));generation=pointer.generation+1;
end
checkpoint=struct('kind','ipm_cr_explicit_accepted_checkpoint_v1','state',state,'rho',state.rho, ...
    'geometry',g,'wallOps',wallOps,'explicitOptions',method,'runnerOptions',opts, ...
    'history',history,'historyColumns',{historyColumns},'sourceInventory',manifest, ...
    'spatialAudit',spatialAudit,'auditHistory',auditHistory, ...
    'authority',authority,'accumulatedWallSeconds',wall,'sequence',state.step, ...
    'acceptedOnly',true,'solverFactorCacheSaved',false);
path=fullfile(folder,sprintf('accepted_%d.mat',slot));temporary=[path,'.writing.mat'];
save(temporary,'checkpoint','-v7.3');check=load(temporary,'checkpoint');
assert(isequaln(check.checkpoint,checkpoint),'ipm_cr:ExplicitCheckpointReadback', ...
    'Accepted checkpoint readback differs from its complete payload.');
[ok,message]=movefile(temporary,path,'f');assert(ok,'ipm_cr:ExplicitCheckpointRename','%s',message);
pointer=struct('kind','ipm_cr_explicit_accepted_pointer_v1','path',path,'generation',generation, ...
    'sequence',state.step,'tau',state.scale.canonicalTime,'sha256',ipm_c.output.sha256(path));
pointer.accumulatedWallSeconds=wall+toc(commitTimer);
write_json(fullfile(folder,'latest.json'),pointer);
end

function [checkpoint,pointer]=load_accepted(folder)
pointer=jsondecode(fileread(fullfile(folder,'latest.json')));
assert(strcmp(pointer.kind,'ipm_cr_explicit_accepted_pointer_v1')&& ...
    any(strcmp(pointer.path,{fullfile(folder,'accepted_0.mat'),fullfile(folder,'accepted_1.mat')}))&& ...
    strcmp(ipm_c.output.sha256(pointer.path),pointer.sha256),'ipm_cr:ExplicitCheckpointIdentity');
loaded=load(pointer.path,'checkpoint');checkpoint=loaded.checkpoint;
assert(strcmp(checkpoint.kind,'ipm_cr_explicit_accepted_checkpoint_v1')&&checkpoint.acceptedOnly&& ...
    checkpoint.sequence==pointer.sequence&&checkpoint.state.step==pointer.sequence&& ...
    checkpoint.state.scale.canonicalTime==pointer.tau&&isequaln(checkpoint.rho,checkpoint.state.rho), ...
    'ipm_cr:ExplicitCheckpointState');
end

function s=finish_summary(s,status,state,g,method,history,names,pointer,priorWall,elapsed)
s.status=status;s.endedUtc=utc_now();s.tau=state.scale.canonicalTime;s.steps=state.step;
s.invocationWallSeconds=elapsed;s.accumulatedWallSeconds=priorWall+elapsed;
s.physicalTime=state.scale.physicalTime;s.physicalTimeLow=state.scale.physicalTimeLow;
s.rates=[state.flow.c_l,state.flow.c_omega,state.flow.c_r];s.explicitOptions=method;
if isfield(method,'timeIntegrator'),s.timeIntegrator=method.timeIntegrator;else,s.timeIntegrator='ssprk3';end
s.coreSize=size(state.rho.core);s.polarSize=size(state.rho.polar);s.activeDofCount=g.activeDofCount;
s.acceptedCheckpoint=pointer;s.historyRows=size(history,1);s.historyColumns=names;
s.maximumAcceptedStageCfl=finite_max(history(:,strcmp(names,'maximumStageCfl')));
s.maximumRecordedGaugeResidual=finite_max(history(:,strcmp(names,'gaugeResidualInfinity')));
s.maximumRecordedPoissonBackwardError=finite_max(history(:,strcmp(names,'poissonBackwardError')));
s.maximumRecordedStrongResidual=finite_max(history(:,strcmp(names,'poissonStrongResidual')));
s.strongResidualAvailableRows=nnz(history(:,strcmp(names,'poissonStrongAvailable')));
s.totalStageCflRetries=sum(history(:,strcmp(names,'retries')));
s.totalRhsEvaluations=sum(history(:,strcmp(names,'rhsEvaluations')));
s.originalPeakEventLimitedSteps=sum(history(:,strcmp(names,'peakEventLimited')));
s.originalPeakEventsLocated=sum(history(:,strcmp(names,'peakEventLocated')));
s.originalPeakEventProbes=sum(history(:,strcmp(names,'peakEventProbeCount')));
s.originalPeakEventExtraRhsEvaluations=sum(history(:,strcmp(names,'peakEventExtraRhs')));
s.nonuniquePeakStageCount=sum(history(:,strcmp(names,'peakTieStages')));
curvature=history(:,strcmp(names,'minimumPeakCurvatureMagnitude'));
curvature=curvature(isfinite(curvature));s.minimumRecordedPeakCurvatureMagnitude=min(curvature);
s.finalDensityRange=history(end,ismember(names,{'totalDensityMin','totalDensityMax'}));
s.finalQInfinity=max([max(abs(state.flow.Q.core),[],'all'),max(abs(state.flow.Q.polar),[],'all')]);
s.targetReached=state.scale.canonicalTime>=s.targetTau;
s.finalPoissonDiagnostics=state.flow.poisson;
end

function write_summary(folder,summary,history,historyColumns)
temporary=fullfile(folder,'summary.writing.mat');save(temporary,'summary','history','historyColumns','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'summary.mat'),'f');assert(ok,'ipm_cr:ExplicitSummary','%s',message);
write_json(fullfile(folder,'summary.json'),summary);
end

function files=source_inventory(root)
trees={fullfile(root,'test','+ipm_cr','+explicit'),fullfile(root,'test','+ipm_cr','+core'), ...
    fullfile(root,'test','+ipm_cr','+polar'),fullfile(root,'test','+ipm_cr','+internal'), ...
    fullfile(root,'test','+ipm_cr','+carrier'),fullfile(root,'time_evolution','+ipm'), ...
    fullfile(root,'time_evolution','+ipm_c')};
paths={mfilename('fullpath'),fullfile(root,'test','+ipm_cr','Solver.m'),fullfile(root,'test','+ipm_cr','config.m')};
paths{1}=[paths{1},'.m'];
for k=1:numel(trees)
    listing=dir(fullfile(trees{k},'**','*.m'));
    for j=1:numel(listing),paths{end+1}=fullfile(listing(j).folder,listing(j).name);end %#ok<AGROW>
end
mexFile=which('ipm_c.field.weno5UnitPairMex');if ~isempty(mexFile),paths{end+1}=mexFile;end
paths{end+1}=fullfile(root,'test','+ipm_cr','+explicit','wenoPairSmoothMex.cpp');
smoothMex=which('ipm_cr.explicit.wenoPairSmoothMex');if ~isempty(smoothMex),paths{end+1}=smoothMex;end
paths=sort(unique(paths));files=struct([]);
for k=1:numel(paths)
    entry=struct('path',paths{k},'sha256',ipm_c.output.sha256(paths{k}));
    if isempty(files),files=entry;else,files(end+1)=entry;end %#ok<AGROW>
end
end

function value=authority_identity()
cfg=ipm_c.config();value=struct('file',cfg.authorityFile,'sha256',ipm_c.output.sha256(cfg.authorityFile));
end

function write_json(path,value)
temporary=[path,'.writing'];fid=fopen(temporary,'w');assert(fid>=0,'ipm_cr:ExplicitJson');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,path,'f');assert(ok,'ipm_cr:ExplicitJson','%s',message);
end

function value=get(s,name,default)
if isfield(s,name),value=s.(name);else,value=default;end
end

function value=finite_max(values)
values=values(isfinite(values));if isempty(values),value=NaN;else,value=max(values);end
end

function value=utc_now()
value=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

function release_owner(lock,channel,file)
lock.release();channel.close();file.close();
end
