function report=solve(cfg,folder)
%SOLVE Direct CFL evolution, accepted-state checkpoints and trajectory output.
assert(startsWith(folder,filesep),'ipm_c:AbsoluteOutput');
assert(cfg.targetTau>0&&cfg.wallSeconds>0&&cfg.maxSteps>=1,'ipm_c:SharedConfiguration');
if ~isfolder(folder),mkdir(folder);end
file=java.io.RandomAccessFile(fullfile(folder,'owner.lock'),'rw');channel=file.getChannel();lock=channel.tryLock();
assert(~isempty(lock),'ipm_c:LiveOwner');cleanup=onCleanup(@()release(lock,channel,file)); %#ok<NASGU>
threads=maxNumCompThreads(cfg.runtime.threads);threadCleanup=onCleanup(@()maxNumCompThreads(threads)); %#ok<NASGU>
timer=tic;manifest=ipm_c.output.manifest();
if isfile(fullfile(folder,'latest.json'))
    [cp,pointer]=ipm_c.output.loadCheckpoint(folder);
    old=jsondecode(fileread(fullfile(folder,'manifest.json')));
    assert(isequaln(old,manifest),'ipm_c:SharedResumeSourceChanged');
    saved=cp.cfg;for field={'targetTau','wallSeconds','maxSteps'},saved.(field{1})=cfg.(field{1});end
    assert(isequaln(saved,cfg),'ipm_c:SharedResumeConfiguration');
    authority=cp.authority;state=cp.state;history=cp.history;holdStart=cp.holdStart;sequence=pointer.sequence;
    state.config=cfg;state.config.time.finalTime=cfg.targetTau;state.config.time.maxSteps=cfg.maxSteps;
    state.ops=ipm_c.spatial.build(cfg,authority.cart,cp.grid);
    state.ops.spatialOperators.remesh=@ipm_c.remesh.apply;
    state.ops.rescaling=cp.rescaling;state.ops.remeshCount=history(end,13);
    [F,state.flow]=ipm.evolve.flow(state.rho,state.ops,state.scale);
    state.rhsCache=ipm.evolve.makeRhsCache(state.rho,F,state.flow,state.ops,state.scale);
    state.timeStep=ipm.evolve.emptyTimestep(state.config.time,state.scale,'restored');
else
    assert(~isfile(fullfile(folder,'latest_shared.json')),'ipm_c:LegacyCheckpoint', ...
        'An archived solver checkpoint cannot be silently converted to this direct CFL architecture.');
    [state,context]=ipm_c.initialize(cfg);
    authority=struct('cart',context.cart,'shared',context.shared, ...
        'authoritySha256',ipm_c.output.sha256(cfg.authorityFile));
    state.config.time.finalTime=cfg.targetTau;state.config.time.maxSteps=cfg.maxSteps;
    history=[];holdStart=NaN;sequence=0;
    if max(abs(state.flow.fullQ),[],'all')<cfg.stationaryQ,holdStart=state.scale.canonicalTime;end
    ipm_c.output.atomicJson(fullfile(folder,'manifest.json'),manifest);
    save(fullfile(folder,'settings_used.mat'),'cfg','authority','-v7.3');
end
assert(strcmp(ipm_c.output.sha256(cfg.authorityFile),authority.authoritySha256),'ipm_c:AuthorityChanged');
status='running';startStep=state.step;startTau=state.scale.canonicalTime;
historyCount=size(history,1);monitorWritten=0;
append_monitor(NaN);publish(true);lastCheckpointTau=state.scale.canonicalTime;
nextOutput=(floor(state.scale.canonicalTime/cfg.time.outputEvery)+1)*cfg.time.outputEvery;
lastCheckpointWall=toc(timer);
try
    while strcmp(status,'running')&&state.scale.canonicalTime<cfg.targetTau&&state.step<cfg.maxSteps&&toc(timer)<cfg.wallSeconds
        if isfile(fullfile(folder,'STOP')),status='stopped';break;end
        [trial,reason]=ipm_c.evolve.advance(state);
        if ~isempty(reason),status=reason;break;end
        state=trial;sequence=sequence+1;append_monitor(state.timeStep.dt);
        if max(abs(state.flow.fullQ),[],'all')<cfg.stationaryQ
            if isnan(holdStart),holdStart=state.scale.canonicalTime;end
        else,holdStart=NaN;end
        due=state.scale.canonicalTime>=nextOutput;
        if due,nextOutput=(floor(state.scale.canonicalTime/cfg.time.outputEvery)+1)*cfg.time.outputEvery;end
        checkpointDue=state.scale.canonicalTime-lastCheckpointTau>=cfg.output.checkpointEvery|| ...
            toc(timer)-lastCheckpointWall>=cfg.output.checkpointWallSeconds;
        if due,ipm_c.output.frame(folder,state,history(1:historyCount,:));end
        if checkpointDue,publish(false);lastCheckpointTau=state.scale.canonicalTime;lastCheckpointWall=toc(timer);end
        if due||checkpointDue
            report=summary();ipm_c.output.atomicJson(fullfile(folder,'report.json'),report);
        end
    end
    if strcmp(status,'running')
        if state.scale.canonicalTime>=cfg.targetTau,status='target_reached';else,status='segment_budget';end
    end
catch problem
    status='failed';publish(true);report=summary();
    report.error=struct('identifier',problem.identifier,'message',problem.message);
    ipm_c.output.atomicJson(fullfile(folder,'report.json'),report);rethrow(problem);
end
publish(true);assert(isequaln(manifest,ipm_c.output.manifest()),'ipm_c:SourceChanged');
report=summary();ipm_c.output.atomicJson(fullfile(folder,'report.json'),report);
    function append_monitor(h)
        f=state.flow;
        row=[state.step,state.scale.canonicalTime,state.scale.physicalTime,h, ...
            max(abs(f.fullQ),[],'all'),max(abs(state.rhsCache.rhoRate+f.densityOriginRate),[],'all'), ...
            f.c_l,f.c_omega,f.c_r,max(abs(state.rho+state.scale.densityOrigin),[],'all'),f.poissonResidual, ...
            max(abs(f.gaugeResidual)),state.ops.remeshCount, ...
            state.scale.logC_l,state.scale.logC_omega,state.scale.X_shift, ...
            f.rhoXInf,f.rhoYInf,f.velocityInf, ...
            sum((state.rho+state.scale.densityOrigin).*state.ops.integrationWeights,'all')];
        if historyCount==0||history(historyCount,1)~=state.step
            historyCount=historyCount+1;
            if historyCount>size(history,1),history(end+4096,numel(row))=NaN;end
            history(historyCount,:)=row;
        end
        if cfg.output.verbose&&(isnan(h)||state.scale.canonicalTime>=nextOutput)
            fprintf('ipm_c CFL step=%d tau=%.12g t=%.12g H=%.5g Q=%.5g remesh=%d\n', ...
                state.step,state.scale.canonicalTime,state.scale.physicalTime,h,row(5),state.ops.remeshCount);
        end
    end
    function publish(frame)
        acceptedHistory=history(1:historyCount,:);
        ipm_c.output.saveCheckpoint(folder,state,authority,cfg,acceptedHistory,holdStart,sequence);
        csv=fullfile(folder,'monitor.csv');
        if monitorWritten==0
            target=fullfile(folder,'monitor.writing.csv');mode='w';
        else,target=csv;mode='a';end
        fid=fopen(target,mode);assert(fid>=0,'ipm_c:MonitorOpen');
        if monitorWritten==0
            fprintf(fid,'step,tau,t,H,Qinf,Finf,cl,cw,cr,rhoInf,poisson,gauge,remesh,logC_l,logC_omega,X_shift,rhoXinf,rhoYinf,uInf,mass\n');
        end
        format=[repmat('%.17g,',1,size(acceptedHistory,2)-1),'%.17g\n'];
        fprintf(fid,format,acceptedHistory(monitorWritten+1:end,:).');fclose(fid);
        if monitorWritten==0,movefile(target,csv,'f');end
        monitorWritten=historyCount;
        ipm_c.output.frame(folder,state,acceptedHistory,frame);
    end
    function result=summary()
        result=struct('kind',cfg.kind,'architecture',cfg.architecture,'status',status, ...
            'method',cfg.time.timeIntegrator,'sharedStepFunction','ipm.evolve.stepRk', ...
            'advanceFunction','ipm_c.evolve.advance','timeControl','full_polar_cfl', ...
            'primaryState','density_minus_origin_with_scalar_origin','tau',state.scale.canonicalTime,'physicalTime',state.scale.physicalTime, ...
            'physicalTimeLow',state.scale.physicalTimeLow,'densityOrigin',state.scale.densityOrigin, ...
            'targetTau',cfg.targetTau,'steps',state.step,'acceptedThisSegment',state.step-startStep, ...
            'canonicalTimeAdvanced',state.scale.canonicalTime-startTau,'elapsedSeconds',toc(timer), ...
            'fullQ',max(abs(state.flow.fullQ),[],'all'),'rates',[state.flow.c_l,state.flow.c_omega,state.flow.c_r], ...
            'remeshCount',state.ops.remeshCount,'lastTimeStep',state.timeStep, ...
            'poissonStrongResidual',state.flow.poissonSolveInfo.maximumOffGridStrongResidualScaledBySource, ...
            'poissonCoreStrongResidual',state.flow.poissonSolveInfo.coreOffGridStrongResidualScaledBySource, ...
            'velocityInfinity',state.flow.velocityInf,'cflPosition',state.flow.ipmCInfo.cflPosition, ...
            'gaugeRootCondition',state.flow.generatorRootCondition, ...
            'wallIntegral',state.flow.wallIntegral,'wallPeak',state.flow.globalWallPeak, ...
            'wallCentroid',state.flow.wallCentroid,'gaugeResidual',state.flow.gaugeResidual, ...
            'wallIntegralRelativeDrift',state.flow.wallIntegralRelativeDrift, ...
            'wallPeakX',state.flow.ipmCInfo.gauge.initialMeasurement.peak.x, ...
            'exponentialAuxiliarySystemUsed',false,'nonlinearStageSolverUsed',false, ...
            'stepRetryUsed',false,'spatialTrajectoryQualified',false, ...
            'longTimeQualified',false,'stationaryHoldPassed', ...
            isfinite(holdStart)&&state.scale.canonicalTime-holdStart>=cfg.holdTau);
    end
end
function release(lock,channel,file)
lock.release();channel.close();file.close();
end
