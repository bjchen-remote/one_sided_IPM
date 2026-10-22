function result = solve(userOpts,restartCheckpoint,datasetPolicy)
%IPM.SOLVE Solve the physical 2-D incompressible porous media equation.
%   RESULT = IPM.SOLVE(OPTIONS) is the only solver entry point. It resolves
%   configuration once, advances a shared numerical path, and returns the
%   grouped version-2 result contract.
%   RESULT = IPM.SOLVE(OVERRIDES,CHECKPOINT) resumes a validated same-grid
%   accepted-state checkpoint. Only terminal horizons and output controls
%   may be overridden; every numerical and physical choice remains frozen.
%   A third optional DATASETPOLICY controls native plot output independently
%   of the checkpoint configuration (enabled, file, frameEvery). By default
%   saved runs produce data.mat plus immutable native frame files.
%
%   The maintained equation is
%       rho_t + u dot grad(rho) = 0,
%       u = grad^perp (-Delta)^(-1) d_x rho
%   on a nonperiodic truncation of R x R_+. Isotropic scaling is the exact
%   C_x=C_y specialization of the two-scale geometry.

if nargin < 1
    userOpts = [];
end
if nargin < 2, restartCheckpoint = []; end
if nargin < 3, datasetPolicy = struct(); end
if nargin < 2 || isempty(restartCheckpoint)
    state = ipm.evolve.initialize(userOpts);
    log = ipm.output.initializeLog();
    [log,~] = ipm.output.record(log,state);
    cursor = initial_cursor(state.config);
else
    [state,log,cursor] = ipm.output.restoreCheckpoint( ...
        restartCheckpoint,userOpts);
end
time = state.config.time;
output = state.config.output;
dataset = ipm.output.dataset('open',state,log,datasetPolicy,restartCheckpoint);
dataset = publish_dataset(dataset,state,log,cursor,true,'initial_or_resumed_state');

viz = [];
if output.livePlot || output.writeVideo
    viz = ipm.output.visualize(viz,state.rho,state.flow,state.ops, ...
        state.scale,log.history,output,state.runMetadata);
    if output.writeVideo && isfield(viz,'videoFile')
        state.runMetadata.videoFile = string(viz.videoFile);
    end
end
videoCleanup = onCleanup(@()ipm.output.closeLive(viz));

stopReason = rho_x_stop_reason(state);
rateTarget = full_rate_target(state);
belowTargetSince = NaN;
if ~isempty(rateTarget)
    if observe_full_rate(state,rateTarget)
        belowTargetSince = restored_hold_start(state,rateTarget);
    end
    state = remember_hold(state,rateTarget,belowTargetSince);
    if isempty(stopReason) && isfinite(belowTargetSince) && ...
            state.scale.canonicalTime-belowTargetSince >= ...
            rateTarget.holdTau-10*eps(max(1,state.scale.canonicalTime))
        stopReason = 'rho_x_time_derivative_full_target';
    end
end
while isempty(stopReason) && ipm.evolve.isActive(state)
    try
        [state,stepStopReason,meshPlan] = ipm.evolve.advance(state);
    catch advanceError
        % Preserve the last accepted state if a trial stage fails. Keep
        % the original numerical error even if checkpoint I/O also fails.
        try
            [log,~] = ipm.output.recordIfNew(log,state);
            [state,cursor] = ipm.output.maybeCheckpoint( ...
                state,log,cursor,true);
            dataset = publish_dataset(dataset,state,log,cursor,true,'advance_error');
        catch checkpointError
            warning('ipm:FailureCheckpoint', ...
                'Could not persist all accepted-state outputs: %s', ...
                checkpointError.message);
        end
        rethrow(advanceError);
    end
    if ~isempty(meshPlan)
        % Own the factor lifetime here: advance has returned, so neither its
        % pre-step rollback state nor this state can retain a second old LU.
        state.ops = rmfield(state.ops,'poisson');
        [state,applied,attempts] = ipm.evolve.applyAutonomousMesh(state,meshPlan);
        if applied
            state = ipm.evolve.publishRemeshWarnings(state);
        else
            if isfield(state.ops,'separablePoisson')
                state.ops.poisson = [];
            else
                state.ops.poisson = decomposition(state.ops.A,'lu');
            end
            state.runMetadata.autonomousMesh.lastFailure = struct( ...
                'kind','candidate_transaction_failure','attempts',attempts, ...
                'step',state.step,'canonicalTime',state.scale.canonicalTime);
            if any(state.config.remesh.autonomousMesh.version == [2,3,4,5])
                state.runMetadata.autonomousMesh.lastFailure.controllerDecision= ...
                    rmfield(meshPlan,{'candidates','axisReport'});
                state.runMetadata.autonomousMesh.lastFailure.axisReport=meshPlan.axisReport;
            end
            stepStopReason = 'autonomous_mesh_candidates_rejected';
            if ~isempty(attempts) && ~isempty(attempts(end).errorIdentifier)
                stepStopReason = 'autonomous_mesh_error';
            end
        end
    end
    if isempty(stepStopReason)
        stepStopReason = rho_x_stop_reason(state);
    end
    if ~isempty(stepStopReason)
        stopReason = stepStopReason;
        [log,recorded] = ipm.output.recordIfNew(log,state);
        [state,cursor] = ipm.output.maybeCheckpoint( ...
            state,log,cursor,true);
        dataset = publish_dataset(dataset,state,log,cursor,true,stopReason);
        if recorded && ~isempty(viz)
            viz = ipm.output.visualize(viz,state.rho,state.flow, ...
                state.ops,state.scale,log.history,output,state.runMetadata);
        end
        break;
    end

    belowTarget = false;
    if ~isempty(rateTarget)
        belowTarget = full_rate_probe(state,rateTarget);
        if belowTarget && isnan(belowTargetSince)
            belowTargetSince = state.scale.canonicalTime;
        elseif ~belowTarget
            belowTargetSince = NaN;
        end
        state = remember_hold(state,rateTarget,belowTargetSince);
    end
    if belowTarget && state.scale.canonicalTime-belowTargetSince >= ...
            rateTarget.holdTau-10*eps(state.scale.canonicalTime)
        stopReason = 'rho_x_time_derivative_full_target';
        [log,recorded] = ipm.output.recordIfNew(log,state);
        observe_full_rate(state,rateTarget);
        [state,cursor] = ipm.output.maybeCheckpoint( ...
            state,log,cursor,true);
        dataset = publish_dataset(dataset,state,log,cursor,true,stopReason);
        if recorded && ~isempty(viz)
            viz = ipm.output.visualize(viz,state.rho,state.flow, ...
                state.ops,state.scale,log.history,output,state.runMetadata);
        end
        break;
    end

    if ~ipm.output.shouldRecord(state,cursor.nextOutput)
        continue;
    end
    [log,~] = ipm.output.record(log,state);
    cursor.nextOutput = cursor.nextOutput+time.outputEvery;
    dataset = publish_dataset(dataset,state,log,cursor,false,'recorded_output');
    if ~isempty(rateTarget)
        observe_full_rate(state,rateTarget);
    end
    [state,cursor] = ipm.output.maybeCheckpoint( ...
        state,log,cursor,false);
    if ~isempty(viz)
        viz = ipm.output.visualize(viz,state.rho,state.flow,state.ops, ...
            state.scale,log.history,output,state.runMetadata);
    end
end

if isempty(stopReason)
    stopReason = 'final_time';
    if state.scale.physicalTime >= time.physicalFinalTime
        stopReason = 'physical_final_time';
    elseif state.step >= time.maxSteps && ...
            state.scale.canonicalTime < time.finalTime
        stopReason = 'maximum_steps';
    end
end

% Always preserve the last finite field, even when it falls between output
% intervals. Snapshot grids remain paired with their corresponding fields.
[log,recorded] = ipm.output.recordIfNew(log,state);
[state,~] = ipm.output.maybeCheckpoint(state,log,cursor,true);
if recorded && ~isempty(viz)
    viz = ipm.output.visualize(viz,state.rho,state.flow,state.ops, ...
        state.scale,log.history,output,state.runMetadata);
end

result = ipm.output.finalize(state,log,stopReason);
if dataset.enabled, result.metadata.datasetFile=string(dataset.file); end
analysis = struct('fit',result.fit,'quality',result.quality);
publish_dataset(dataset,state,log,cursor,true,stopReason,analysis);
result = ipm.output.write(result,output);
ipm.output.report(result);
if output.makePlots && isempty(viz)
    ipm.output.plotResult(result);
end
end

function stream = publish_dataset(stream,state,log,cursor,force,reason,analysis)
% An output failure preserves the accepted checkpoint, then fails visibly.
% Export state never enters the PDE state, signed checkpoint or RK cache.
if nargin<7, analysis=struct(); end
try
    stream=ipm.output.dataset('record',stream,state,log,force,reason,analysis);
catch outputError
    try
        ipm.output.maybeCheckpoint(state,log,cursor,true);
    catch checkpointError
        outputError=addCause(outputError,checkpointError);
    end
    rethrow(outputError);
end
end

function stopReason = rho_x_stop_reason(state)
maximum = ipm.diagnostics.physicalRhoXInf(state.flow,state.scale);
threshold = state.config.diagnostics.rhoXStop;
if isfinite(threshold) && maximum >= threshold
    stopReason = 'rho_x_threshold';
else
    stopReason = '';
end
end

function cursor = initial_cursor(config)
nextCheckpoint = Inf;
if isfield(config.output,'checkpoint') && ...
        config.output.checkpoint.enabled
    nextCheckpoint = config.output.checkpoint.every;
end
cursor = struct('nextOutput',config.time.outputEvery, ...
    'nextCheckpoint',nextCheckpoint,'lastCheckpointStep',-1, ...
    'lastCheckpointCanonicalTime',NaN);
end

function target = full_rate_target(state)
% Explicit archived-run metadata enables this observation-only stopping test.
target = [];
metadata = state.config.output.caseMetadata;
if ~isstruct(metadata) || ...
        ~isfield(metadata,'fullRhoXRateTarget')
    return;
end
threshold = metadata.fullRhoXRateTarget;
validateattributes(threshold,{'double'}, ...
    {'scalar','real','finite','positive'});
assert(isfield(metadata,'fullRhoXRateMonitorFile'), ...
    'A full-rate target requires a monitor CSV path.');
file = char(metadata.fullRhoXRateMonitorFile);
assert(~isempty(file),'The monitor CSV path cannot be empty.');
holdTau = 0.5;
if isfield(metadata,'fullRhoXRateHoldTau')
    holdTau = metadata.fullRhoXRateHoldTau;
end
validateattributes(holdTau,{'double'}, ...
    {'scalar','real','finite','nonnegative'});
target = struct('threshold',threshold,'holdTau',holdTau,'file',file);
end

function reached = observe_full_rate(state,target)
% The accepted-state RHS cache is the same semi-discrete RHS used by RK.
[reached,rate] = full_rate_probe(state,target);
source = state.flow.source;
[fullRate,index] = max(abs(rate),[],'all');
[rateY,rateX] = ind2sub(size(rate),index);
coreX = state.ops.x>=-2 & state.ops.x<=2;
coreY = state.ops.y>=0 & state.ops.y<=2;
assert(any(coreX) && any(coreY),'The fixed core has no grid nodes.');
coreRate = max(abs(rate(coreY,coreX)),[],'all');
leftRate = max(abs(rate(:,1)));
rightRate = max(abs(rate(:,end)));
topRate = max(abs(rate(end,:)));
boundaryRate = max([leftRate,rightRate,topRate]);
leftSource = max(abs(source(:,1)));
rightSource = max(abs(source(:,end)));
topSource = max(abs(source(end,:)));
boundarySource = max([leftSource,rightSource,topSource]);
fullSource = max(abs(source),[],'all');
denominator = state.flow.c_l-state.flow.c_omega;
cs = NaN;
if abs(denominator) > 100*eps(max([1,abs(state.flow.c_l),abs(state.flow.c_omega)]))
    cs = state.flow.c_l/denominator;
end
if strcmp(state.config.scaling.scalingContract,'wall_integral_centroid_practical_v1')
    peak=state.flow.globalWallPeak;
else
    peak=state.flow.omegaGaugeQuadraticPeakValue;
end
assert(isfinite(peak) && peak>0,'Invalid wall amplitude in the full-rate monitor.');
columnMax = max(abs(source),[],1)/peak;
frontIndex = find(columnMax>=1e-5,1,'first');
assert(~isempty(frontIndex),'The 1e-5 front is not present.');
frontX = state.ops.x(frontIndex);
leftInflow = nnz(state.flow.transportU1(:,1)>0);
rightInflow = nnz(state.flow.transportU1(:,end)<0);
topInflow = nnz(state.flow.transportU2(end,:)<0);
assert(all(isfinite([fullRate,coreRate,boundaryRate, ...
    boundarySource,fullSource,frontX])), ...
    'Nonfinite full-rate monitor value.');
newFile = ~isfile(target.file);
[fid,message] = fopen(target.file,'a');
assert(fid>0,'Cannot open full-rate monitor: %s',message);
cleanup = onCleanup(@()fclose(fid));
if newFile
    fprintf(fid,['tau,step,fullRateNorm,coreRateNorm,fullSourceNorm,' ...
        'leftRateRatio,rightRateRatio,topRateRatio,boundaryRateRatio,' ...
        'leftSourceRatio,rightSourceRatio,topSourceRatio,' ...
        'boundarySourceRatio,Cs,front1e5X,frontDistanceToLeft,' ...
        'leftInflowNodes,rightInflowNodes,topInflowNodes,' ...
        'maxRateX,maxRateY,belowTarget,targetHoldTau\n']);
end
fprintf(fid,['%.17g,%d,%.17g,%.17g,%.17g,' ...
    '%.17g,%.17g,%.17g,%.17g,' ...
    '%.17g,%.17g,%.17g,%.17g,' ...
    '%.17g,%.17g,%.17g,%d,%d,%d,' ...
    '%.17g,%.17g,%d,%.17g\n'], ...
    state.scale.canonicalTime,state.step,fullRate,coreRate,fullSource, ...
    leftRate/max(fullRate,realmin), ...
    rightRate/max(fullRate,realmin), ...
    topRate/max(fullRate,realmin), ...
    boundaryRate/max(fullRate,realmin), ...
    leftSource/max(fullSource,realmin), ...
    rightSource/max(fullSource,realmin), ...
    topSource/max(fullSource,realmin), ...
    boundarySource/max(fullSource,realmin),cs, ...
    frontX,frontX-state.ops.x(1), ...
    leftInflow,rightInflow,topInflow, ...
    state.ops.x(rateX),state.ops.y(rateY),reached,target.holdTau);
end

function [reached,rate] = full_rate_probe(state,target)
% This touches no PDE state, Poisson solve, or Runge--Kutta stage.
cache = state.rhsCache;
assert(isequaln(cache.rho,state.rho) && ...
    isequaln(cache.scale,state.scale) && ...
    isequaln(cache.x,state.ops.x) && isequaln(cache.y,state.ops.y) && ...
    cache.remeshCount==state.ops.remeshCount, ...
    'The accepted-state RHS cache does not match the observed state.');
rate = cache.rhoRate*state.ops.Dx';
assert(isequal(size(rate),[state.ops.ny,state.ops.nx]) && ...
    all(isfinite(rate),'all'),'Invalid full-rate field.');
reached = max(abs(rate),[],'all')<=target.threshold;
if reached
    % A small residual on an unresolved or oscillatory field is not a
    % convergence certificate. Reuse the existing shared quality contract
    % only near the target; no extra Poisson solve or stored snapshot.
    probe = state;
    probe.config.output.verbose = false;
    probe.config.output.storeSnapshots = false;
    audit = ipm.output.initializeLog();
    [audit,~] = ipm.output.record(audit,probe);
    reached = all(ipm.output.trustedMask(audit.history,state.config));
end
end

function start = restored_hold_start(state,target)
% Preserve a continuous convergence-observation interval across segments.
start = state.scale.canonicalTime;
if ~isfield(state.runMetadata,'convergenceObservation')
    return
end
saved = state.runMetadata.convergenceObservation;
names = {'threshold','holdTau','belowTargetSince','lastCanonicalTime'};
if ~isstruct(saved) || ~isscalar(saved) || ~all(isfield(saved,names))
    return
end
for index = 1:numel(names)
    value = saved.(names{index});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
        return
    end
end
if saved.threshold ~= target.threshold || saved.holdTau ~= target.holdTau || ...
        saved.lastCanonicalTime ~= state.scale.canonicalTime || ...
        saved.belowTargetSince < 0 || saved.belowTargetSince > state.scale.canonicalTime
    return
end
start = saved.belowTargetSince;
end

function state = remember_hold(state,target,start)
state.runMetadata.convergenceObservation = struct( ...
    'threshold',target.threshold,'holdTau',target.holdTau, ...
    'belowTargetSince',start,'lastCanonicalTime',state.scale.canonicalTime);
end
