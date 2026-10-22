function [state,log,cursor] = restoreCheckpoint(inputCheckpoint,overrides)
%IPM.OUTPUT.RESTORECHECKPOINT Rebuild one same-grid trusted solver state.
%   Only terminal horizons and output controls may be overridden. Numerical,
%   physical, scaling, diagnostic, and grid choices remain those checkpointed.
%   Runtime references are allowlisted onto freshly derived schema-4 ops.

if nargin < 2 || isempty(overrides)
    overrides = struct();
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('ipm:CheckpointOverrides', ...
        'Checkpoint resume overrides must be a flat scalar structure.');
end
allowed = {'finalTime','physicalFinalTime','maxSteps','rhoXStop', ...
    'cfl','timeIntegrator', ...
    'saveResults','resultFile','makePlots','livePlot','writeVideo', ...
    'videoFile','videoFrameRate','videoQuality','verbose','caseMetadata', ...
    'checkpoint'};
unexpected = setdiff(fieldnames(overrides),allowed,'stable');
if ~isempty(unexpected)
    error('ipm:CheckpointOverrideNotAllowed', ...
        'Checkpoint resume cannot override: %s.',strjoin(unexpected,', '));
end

checkpoint = ipm.output.readCheckpoint(inputCheckpoint);
payload = checkpoint.payload;
saved = payload.state;
flat = flatten_config(saved.config);
names = fieldnames(overrides);
for index = 1:numel(names)
    flat.(names{index}) = overrides.(names{index});
end
config = ipm.config.resolve(flat);
assert_open_horizons(config.time,saved);

if isfield(config.remesh,'autonomousMesh') && config.remesh.autonomousMesh.enabled
    % The validated controller records custom-axis construction from t=0.
    % Preserve that operator representation even when axes equal the analytic
    % configured coordinates; this also avoids an unused initial-grid LU.
    grid = config.grid;
    grid.nx = numel(saved.x);
    grid.ny = numel(saved.y);
    grid.customX = saved.x;
    grid.customY = saved.y;
    ops = ipm.mesh.build(config,grid);
else
    ops = ipm.mesh.build(config);
    if ~isequal(ops.x,saved.x) || ~isequal(ops.y,saved.y)
        grid = config.grid;
        grid.nx = numel(saved.x);
        grid.ny = numel(saved.y);
        grid.customX = saved.x;
        grid.customY = saved.y;
        ops = ipm.mesh.build(config,grid);
    end
end
ops.baseX = saved.baseX;
ops.baseY = saved.baseY;
ops.rescaling = ipm.output.restoreRuntimeReferences( ...
    ops.rescaling,saved.rescaling,ops.x, ...
    strcmp(config.scaling.rescalingMode,'dynamic'));
ops.remeshCount = saved.remeshCount;
if strcmp(config.scaling.cOmegaGauge,'wall_omega_tracked_spline_peak')
    window=ops.rescaling.pinX+[-1,1]*ops.rescaling.peakTrackingHalfWidth;
    ipm.evolve.validateTrackedPeakAnchor(ops.rescaling.trackedPeakAnchor, ...
        ops.x,window,saved.scale.canonicalTime);
    ops.rescaling.peakEvaluationTime=saved.scale.canonicalTime;
end
[rhoRate,flow] = ipm.evolve.flow(saved.rho,ops,saved.scale);
if isfield(flow,'peakTrackingAnchor')
    ops.rescaling.trackedPeakAnchor=flow.peakTrackingAnchor;
end
rhsCache = ipm.evolve.makeRhsCache( ...
    saved.rho,rhoRate,flow,ops,saved.scale);

metadata = update_metadata(saved.runMetadata,config,overrides,checkpoint);
if isfield(ops,'timeIntegrationContract')
    if strcmp(saved.config.time.timeIntegrator,config.time.timeIntegrator)
        assert(isfield(saved.runMetadata,'timeIntegrationContract') && ...
            isequaln(saved.runMetadata.timeIntegrationContract,ops.timeIntegrationContract), ...
            'ipm:ImplicitCheckpointContract', ...
            'The saved versioned implicit solver contract differs from this runtime.');
    else
        metadata.timeIntegratorExperiment=struct('source',saved.config.time.timeIntegrator, ...
            'selected',config.time.timeIntegrator,'explicitOverride',true);
    end
    metadata.timeIntegrationContract=ops.timeIntegrationContract;
    metadata.transportCflPolicy='complete_rate_recorded; implicit_stage_convergence_and_time_accuracy_control';
    metadata.implicitTimeMethodExperimental=true;
end
state = struct('config',config,'runMetadata',metadata, ...
    'ops',ops,'rho',saved.rho,'flow',flow,'scale',saved.scale, ...
    'rhsCache',rhsCache, ...
    'normalizedTime',saved.normalizedTime,'step',saved.step, ...
    'timeStep',ipm.evolve.emptyTimestep( ...
    config.time,saved.scale,'restored'), ...
    'mass0',saved.mass0,'rhoRange0',saved.rhoRange0);
if ~ipm.evolve.isFinite(state)
    error('ipm:CheckpointRestoredState', ...
        'The reconstructed checkpoint state is not finite.');
end
log = payload.log;
if isfield(config.remesh,'autonomousMesh') || isfield(metadata,'autonomousMesh')
    current = struct();
    if isfield(config.remesh,'autonomousMesh') && config.remesh.autonomousMesh.enabled
        current = struct('step',state.step,'canonicalTime',state.scale.canonicalTime, ...
            'physicalTime',state.scale.physicalTime,'normalizedTime',state.normalizedTime, ...
            'remeshCount',ops.remeshCount,'coreCells',[flow.coreGridPoints,flow.verticalCoreGridPoints], ...
            'safety',flow.safetyFactor,'x',ops.x,'y',ops.y, ...
            'baseX',ops.baseX,'baseY',ops.baseY,'history',log.history);
        if any(config.remesh.autonomousMesh.version == [2,3,4,5])
            current.levelId=metadata.autonomousMesh.currentLevelId;
            current.nodeCount=[ops.nx,ops.ny];
            current.snapshots=struct('rho',{log.snapshotRho},'x',{log.snapshotX},'y',{log.snapshotY});
        end
    end
    ipm.remesh.validateController(metadata,config,current);
end
cursor = payload.cursor;
cursor = update_checkpoint_cursor( ...
    cursor,saved.config.output,config.output,saved.scale.canonicalTime,overrides);
end

function flat = flatten_config(config)
schema = ipm.config.schema();
flat = struct();
for domainIndex = 1:numel(schema.domainNames)
    domain = config.(schema.domainNames{domainIndex});
    names = fieldnames(domain);
    for nameIndex = 1:numel(names)
        flat.(names{nameIndex}) = domain.(names{nameIndex});
    end
end
end

function assert_open_horizons(time,saved)
tolerance = 64*eps(max([1,saved.normalizedTime, ...
    saved.scale.canonicalTime,saved.scale.physicalTime]));
if time.finalTime+ tolerance < saved.scale.canonicalTime || ...
        time.physicalFinalTime+tolerance < saved.scale.physicalTime || ...
        time.maxSteps < saved.step
    error('ipm:CheckpointHorizon', ...
        ['Resume horizons must not precede the checkpoint canonical time, ' ...
        'physical time, or step.']);
end
end

function metadata = update_metadata(metadata,config,overrides,checkpoint)
if ~isfield(metadata,'resumeCount')
    metadata.resumeCount = 0;
end
metadata.resumeCount = metadata.resumeCount+1;
metadata.resumedFromStep = checkpoint.payload.state.step;
metadata.resumedFromCanonicalTime = ...
    checkpoint.payload.state.scale.canonicalTime;
metadata.resumedFromPhysicalTime = ...
    checkpoint.payload.state.scale.physicalTime;
if isfield(checkpoint,'storageFile')
    metadata.resumedFromCheckpoint = checkpoint.storageFile;
else
    metadata.resumedFromCheckpoint = "in_memory_checkpoint";
end
metadata.caseMetadata = config.output.caseMetadata;

saveChanged = isfield(overrides,'saveResults') || ...
    isfield(overrides,'resultFile');
if ~config.output.saveResults
    metadata.resultFile = "";
elseif saveChanged
    metadata.resultFile = string(ipm.output.uniquePath( ...
        config.output.resultFile,metadata.caseId));
end
videoChanged = isfield(overrides,'writeVideo') || ...
    isfield(overrides,'videoFile');
if ~config.output.writeVideo
    metadata.videoFile = "";
elseif videoChanged
    metadata.videoFile = string(ipm.output.uniquePath( ...
        config.output.videoFile,metadata.caseId));
end
end

function cursor = update_checkpoint_cursor( ...
        cursor,oldOutput,newOutput,tau,overrides)
if ~isfield(overrides,'checkpoint')
    return;
end
oldPolicy = disabled_policy();
newPolicy = disabled_policy();
if isfield(oldOutput,'checkpoint')
    oldPolicy = oldOutput.checkpoint;
end
if isfield(newOutput,'checkpoint')
    newPolicy = newOutput.checkpoint;
end
cadenceChanged = oldPolicy.enabled ~= newPolicy.enabled || ...
    oldPolicy.every ~= newPolicy.every;
if ~cadenceChanged
    return;
end
if newPolicy.enabled && isfinite(newPolicy.every)
    cursor.nextCheckpoint = tau+newPolicy.every;
else
    cursor.nextCheckpoint = Inf;
end
end

function policy = disabled_policy()
policy = struct('enabled',false,'file','','every',Inf,'atExit',false);
end
