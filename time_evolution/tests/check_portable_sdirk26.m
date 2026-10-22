function report=check_portable_sdirk26(outputRoot,allowLocalExecution)
% Two genuine SDIRK2 steps, one CP replay, and a separate 2-D area fixture.
if nargin<2,allowLocalExecution=false;end
maxNumCompThreads(1); % Test mode has one numerical thread even through main.m.
assert(~ismac || allowLocalExecution);
assert(startsWith(outputRoot,filesep) && ~isfolder(outputRoot));
packageRoot=fileparts(fileparts(mfilename('fullpath')));
assert(startsWith(which('ipm.solve'),[packageRoot,filesep]));
assert(isempty(which('ipm.field.weno5LocalMex')),'A MEX must not supply this portable test.');
manifest=jsondecode(fileread(fullfile(packageRoot,'MANIFEST.json')));
verify_runtime(packageRoot,manifest.runtimeRows);
mkdir(outputRoot);timer=tic;
report=struct('kind','portable_source26_sdirk2_actual_advance_checkpoint_area_acceptance', ...
    'status','initializing','passed',false,'executed',true,'pureMatlab',true, ...
    'threads',maxNumCompThreads, ...
    'sourceOrRhsMaskUsed',false,'goalHoldPassed',false, ...
    'runtimeTreeSha256',manifest.runtimeTreeSha256, ...
    'timeIntegrator','sdirk2_alexander_practical_v2', ...
    'initialCondition','localized_gaussian_primitive', ...
    'gaussianParameters',[.7,0,.75,1,.2], ...
    'requestedMaximumCanonicalStep',.2, ...
    'realTrajectorySteps',0,'realAdvanceCalls',0, ...
    'syntheticAreaFixtureIsPhysicalEvolution',false,'elapsedSeconds',0);
persist(outputRoot,report);
try
    [opts,~]=ipm.workflow.build_fresh_options(settings(),outputRoot);
    opts.saveResults=false;opts.makePlots=false;opts.livePlot=false;opts.writeVideo=false;
    opts.verbose=false;
    opts.checkpoint=struct('enabled',false,'file',fullfile(outputRoot,'disabled.mat'), ...
        'every',Inf,'atExit',false);
    fresh=ipm.evolve.initialize(opts);
    assert(strcmp(fresh.config.time.timeIntegrator,report.timeIntegrator));
    assert(strcmp(fresh.config.physics.initialCondition,report.initialCondition));
    assert(isequal([fresh.config.physics.gaussianAmplitude, ...
        fresh.config.physics.gaussianCenterX,fresh.config.physics.gaussianWidthX, ...
        fresh.config.physics.gaussianWidthY,fresh.config.physics.gaussianSkewWeight], ...
        report.gaussianParameters));
    assert(strcmp(fresh.config.scaling.gaugeTransportAssembly, ...
        'net_velocity_density_centered_advective'));
    assert(strcmp(fresh.config.scaling.scalingContract,'wall_integral_centroid_practical_v1'));
    assert(fresh.config.remesh.adaptiveRemesh && fresh.config.scaling.wallIntegralRadius==1);
    I0=fresh.ops.rescaling.referenceWallIntegral;
    C0=fresh.ops.rescaling.referenceWallGradientCentroid;
    report.referenceI0=I0;report.referenceC0=C0;report.initialAudit=audit(fresh);

    % Genuine trajectory: retain the original analytic datum and references.
    first=advance_one(fresh);report.realAdvanceCalls=1;report.realTrajectorySteps=1;
    log=ipm.output.initializeLog();[log,~]=ipm.output.record(log,first);
    cursor=struct('nextOutput',first.scale.canonicalTime+first.config.time.outputEvery, ...
        'nextCheckpoint',Inf,'lastCheckpointStep',-1,'lastCheckpointCanonicalTime',NaN);
    cp=ipm.output.makeCheckpoint(first,log,cursor);
    [~,checkpointFile]=ipm.output.writeCheckpoint(cp,fullfile(outputRoot,'accepted_step1.mat'));
    report.checkpointFile=checkpointFile;report.checkpointSha256=local_sha(checkpointFile);
    direct=advance_one(first);report.realAdvanceCalls=2;report.realTrajectorySteps=2;
    [restored,~,~]=ipm.output.restoreCheckpoint(checkpointFile);
    assert(restored.step==first.step && isequaln(restored.scale,first.scale));
    replay=advance_one(restored);report.realAdvanceCalls=3;
    assert(direct.step==2 && replay.step==2 && direct.scale.canonicalTime>0);
    assert(abs(first.timeStep.dt-.2)<=8*eps(.2) && ...
        abs(direct.timeStep.dt-.2)<=8*eps(.2) && ...
        abs(direct.scale.canonicalTime-.4)<=16*eps(.4), ...
        'Default-step acceptance requires two actual h=.2 steps.');
    assert(direct.ops.rescaling.referenceWallIntegral==I0 && ...
        replay.ops.rescaling.referenceWallIntegral==I0 && ...
        direct.ops.rescaling.referenceWallGradientCentroid==C0 && ...
        replay.ops.rescaling.referenceWallGradientCentroid==C0);
    qDirect=direct.rhsCache.rhoRate*direct.ops.Dx';
    qReplay=replay.rhsCache.rhoRate*replay.ops.Dx';
    densityError=max(abs(direct.rho-replay.rho),[],'all')/max(1,max(abs(direct.rho),[],'all'));
    qRateError=max(abs(qDirect-qReplay),[],'all');
    names=fieldnames(direct.scale);scaleError=0;
    for k=1:numel(names)
        scaleError=max(scaleError,max(abs(direct.scale.(names{k})-replay.scale.(names{k})),[],'all'));
    end
    assert(densityError<=1e-8 && qRateError<=1e-7*max(1,max(abs(qDirect),[],'all')) && scaleError<=1e-10);
    assert(strcmp(local_sha(checkpointFile),report.checkpointSha256));
    report.genuineEvolution=struct('actualTau',direct.scale.canonicalTime, ...
        'actualPhysicalTime',direct.scale.physicalTime,'actualSteps',direct.step, ...
        'firstAcceptedDt',first.timeStep.dt,'secondAcceptedDt',direct.timeStep.dt, ...
        'firstStageDiagnostics',first.runMetadata.lastImplicitStep, ...
        'secondStageDiagnostics',direct.runMetadata.lastImplicitStep, ...
        'firstAudit',audit(first),'secondAudit',audit(direct), ...
        'replayAudit',audit(replay),'relativeDensityReplayError',densityError, ...
        'completeQRateReplayError',qRateError,'scaleReplayError',scaleError);
    report.status='genuine_steps_and_checkpoint_replay_passed';persist(outputRoot,report);

    % Separate static morphology fixture. Stretch Y by 1.35; wall rho is exact
    % unchanged, so I0 and C0 are retained rather than reset to fit the fixture.
    [~,baselineOps,baselineInfo]=ipm.remesh.adapt(fresh.rho,fresh.ops,fresh.flow,fresh.config,1);
    assert(~baselineInfo.applied);
    shape=fresh;shape.ops=baselineOps;
    analyticOps=shape.ops;analyticOps.Y=shape.ops.Y/1.35;
    shape.rho=ipm.field.initialDensity(analyticOps,shape.config.physics);
    assert(isequal(shape.rho(1,:),fresh.rho(1,:)));
    [shapeRate,shape.flow]=ipm.evolve.flow(shape.rho,shape.ops,shape.scale);
    shape.rhsCache=ipm.evolve.makeRhsCache(shape.rho,shapeRate,shape.flow,shape.ops,shape.scale);
    before=audit(shape);
    [rhoNew,opsNew,info,rateNew,flowNew]=ipm.remesh.adapt( ...
        shape.rho,shape.ops,shape.flow,shape.config,1);
    save(fullfile(outputRoot,'synthetic_area_transfer.mat'),'info','-v7.3');
    assert(info.areaRequested && info.maxAreaChange>=.2 && info.applied);
    assert(opsNew.remeshCount==shape.ops.remeshCount+1 && ...
        (~isequal(opsNew.x,shape.ops.x) || ~isequal(opsNew.y,shape.ops.y)));
    assert(opsNew.rescaling.referenceWallIntegral==I0 && ...
        opsNew.rescaling.referenceWallGradientCentroid==C0);
    transferred=shape;transferred.rho=rhoNew;transferred.ops=opsNew;
    if isempty(rateNew) || isempty(flowNew)
        [rateNew,flowNew]=ipm.evolve.flow(rhoNew,opsNew,transferred.scale);
    end
    transferred.flow=flowNew;
    transferred.rhsCache=ipm.evolve.makeRhsCache(rhoNew,rateNew,flowNew,opsNew,transferred.scale);
    report.areaFixture=struct('kind','independent_static_Y_stretch_not_real_time', ...
        'stretchFactor',1.35,'expectedAreaRelativeChange',.35, ...
        'actualAreaRelativeChange',info.maxAreaChange,'info',info, ...
        'originalI0C0Preserved',true,'beforeAudit',before,'afterAudit',audit(transferred));
    verify_runtime(packageRoot,manifest.runtimeRows);
    assert(isempty(which('ipm.field.weno5LocalMex')));
    report.runtimeBytesUnchanged=true;report.passed=true;
    report.status='portable_two_real_steps_checkpoint_replay_area_trigger_passed_not_goal_hold';
catch e
    report.status='failed';report.errorIdentifier=e.identifier;
    report.errorMessage=getReport(e,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputRoot,report);rethrow(e)
end
report.elapsedSeconds=toc(timer);persist(outputRoot,report);
fprintf('PORTABLE_SDIRK26 passed=%d elapsed=%.3fs\n',report.passed,report.elapsedSeconds);
end

function state=advance_one(state)
previous=state.step;[state,reason,plan]=ipm.evolve.advance(state);
assert(isempty(reason) && isempty(plan) && state.step==previous+1 && ipm.evolve.isFinite(state));
assert(state.runMetadata.lastImplicitStep.stageCount==2);
end
function a=audit(state)
q=state.rho*state.ops.Dx';qd=state.rhsCache.rhoRate*state.ops.Dx';
assert(all(isfinite(q),'all') && all(isfinite(qd),'all'));
peak=ipm.evolve.globalSplineMaximum(q(1,:),state.ops);
centroid=ipm.evolve.wallGradientCentroid(q(1,:),state.ops);
integral=ipm.evolve.wallGradientIntegral(q(1,:),state.ops,state.ops.rescaling);
Idot=integral.weights*qd(1,:)';
Pdot=peak.weights*qd(1,peak.indices)';
Cdot=centroid.weights*qd(1,centroid.indices)';
C0=state.ops.rescaling.referenceWallGradientCentroid;
phase=Cdot+state.ops.rescaling.peakPositionRestoringRate*(centroid.value-C0);
assert(max(abs([Idot,Pdot,phase]))<=1e-5*max(1,abs(peak.value)));
[full,index]=max(abs(qd),[],'all');[iy,ix]=ind2sub(size(qd),index);
core=abs(state.ops.X)<=2 & state.ops.Y<=2;
a=struct('tau',state.scale.canonicalTime,'physicalTime',state.scale.physicalTime, ...
    'step',state.step,'remeshCount',state.ops.remeshCount,'fullQRate',full, ...
    'coreQRate',max(abs(qd(core))),'fullQRateMaximumPosition',[state.ops.x(ix),state.ops.y(iy)], ...
    'I',integral.value,'I0',state.ops.rescaling.referenceWallIntegral,'I_dot',Idot, ...
    'P',peak.value,'P_dot',Pdot,'C',centroid.value,'C0',C0,'phaseResidual',phase);
end
function s=settings()
s.grid=struct('nodesX',513,'nodesY',321,'adaptiveRemesh',true,'spatialOrder','high_order', ...
    'transferScheme','high_order','adaptiveLevels',[.1,.2,.5],'maximumCellRatio',6);
s.box=struct('left',1e10,'right',1e10,'top',1e10);
s.initial=struct('kind','localized_gaussian_primitive','power',8,'cutoffB',8, ...
    'gaussianAmplitude',.7,'gaussianCenterX',0,'gaussianWidthX',.75, ...
    'gaussianWidthY',1,'gaussianSkewWeight',.2);
s.gauge=struct('assembly','net_velocity_density_centered_advective','length','wall_gradient_integral', ...
    'amplitude','wall_omega_global_spline_maximum','position','wall_gradient_centroid', ...
    'contract','wall_integral_centroid_practical_v1','peakPositionRestoringRate',.5,'wallIntegralRadius',1);
s.time=struct('finalTau',1,'integrator','sdirk2_alexander_practical_v2','cfl',.5, ...
    'maxDt',.2,'physicalLimit',Inf,'maxSteps',100,'outputEvery',.1,'checkpointEvery',.5);
s.stop=struct('fullAbsoluteRhoXRate',1e-5,'holdTau',.5);
s.transport=struct('inflowDensityDegree',1,'scheme','weno5_fd','alphaMode','local');
s.boundary=struct('farMode','green','transportMode','open');
s.poisson=struct('tolerance',1e-5,'solver','separable','greenSourceTolerance',1e-12,'greenMaxSources',48000);
s.runtime=struct('threads',1,'sourceCflRel',[]);s.output=struct('storeSnapshots',false,'livePlot',false);
end
function verify_runtime(root,rows)
for k=1:numel(rows),assert(strcmp(local_sha(fullfile(root,rows(k).path)),rows(k).sha256));end
end
function value=local_sha(path)
fid=fopen(path,'rb');assert(fid>0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
function persist(folder,report)
fid=fopen(fullfile(folder,'report.json'),'w');assert(fid>0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
save(fullfile(folder,'report.mat'),'report','-v7.3');
end
