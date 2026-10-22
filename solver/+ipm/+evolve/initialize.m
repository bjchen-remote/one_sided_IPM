function state = initialize(userOpts)
%IPM.EVOLVE.INITIALIZE Resolve configuration and build the initial solver state.

if nargin < 1 || isempty(userOpts)
    userOpts = ipm.config.activeCase();
end
config = ipm.config.resolve(userOpts);
remesh = config.remesh;
output = config.output;
autonomous = isfield(remesh,'autonomousMesh') && remesh.autonomousMesh.enabled;
ops = ipm.mesh.build(config);
runMetadata = ipm.output.initializeMetadata(config,ops);
rho = ipm.field.initialDensity(ops,config.physics);
ops = ipm.evolve.initializeScaling(rho,ops);
scale = ipm.evolve.initialScale(ops);
[rhoRate,flow] = ipm.evolve.flow(rho,ops,scale);

if ipm.mesh.isQuadrant(ops)
    initialRemeshReady=all(isfinite([flow.coreGridPoints, ...
        flow.verticalCoreGridPoints]));
else
    initialRemeshReady=isfinite(flow.safetyFactor);
end
if ~autonomous && remesh.initialAnalyticRemesh && remesh.adaptiveRemesh && ...
        initialRemeshReady
    for pass = 1:remesh.initialAnalyticRemeshPasses
        try
            [~,candidateOps,info] = ...
                ipm.remesh.adapt( ...
                rho,ops,flow,config,exp(scale.logC_omega));
        catch exception
            runMetadata = initial_remesh_warning( ...
                runMetadata,ops,flow,scale, ...
                'remesh_exception',sprintf('%s: %s', ...
                exception.identifier,exception.message));
            break;
        end
        if ~info.applied
            if isfield(info,'status') && ...
                    ~strcmp(info.status,'not_triggered')
                runMetadata = initial_remesh_warning( ...
                    runMetadata,ops,flow,scale, ...
                    'remesh_candidate_rejected',sprintf( ...
                    'initial direct level-set remesh returned %s', ...
                    info.status));
            end
            break;
        end
        % Sample the analytic datum on the accepted grid; interpolating the
        % old samples would preserve their derivative error.
        ops = candidateOps;
        rho = ipm.field.initialDensity(ops,config.physics);
        ops = ipm.evolve.initializeScaling(rho,ops);
        [rhoRate,flow] = ipm.evolve.flow(rho,ops,scale);
        if output.verbose
            fprintf(['  analytic initial remesh %d: peak dx %.3e -> %.3e, ' ...
                'wall dy %.3e -> %.3e\n'],pass,info.oldPeakSpacing, ...
                info.newPeakSpacing,info.oldWallSpacing,info.newWallSpacing);
        end
    end
end

rhsCache = ipm.evolve.makeRhsCache(rho,rhoRate,flow,ops,scale);
state = struct('config',config,'runMetadata',runMetadata, ...
    'ops',ops,'rho',rho,'flow',flow, ...
    'scale',scale,'rhsCache',rhsCache,'normalizedTime',0,'step',0, ...
    'timeStep',ipm.evolve.emptyTimestep(config.time,scale,'initial'), ...
    'mass0',sum(rho.*ops.integrationWeights,'all'), ...
    'rhoRange0',[min(rho,[],'all'),max(rho,[],'all')]);
if autonomous
    [state,plan] = ipm.evolve.planAutonomousMesh(state,true);
    if ~isempty(plan.stopReason) && initial_observer_enabled(config)
        [state,plan] = initial_observation_fallback(state,plan,struct([]));
    end
    if ~isempty(plan.stopReason)
        error('ipm:AutonomousMeshInitialization', ...
            'Zero-time automatic mesh planning failed: %s.',plan.stopReason);
    end
    % This scope owns the initial factor. Remove its local alias as well as
    % the state reference before constructing any candidate factor.
    clear ops
    state.ops = rmfield(state.ops,'poisson');
    [state,success,attempts] = ipm.evolve.applyAutonomousMesh(state,plan);
    if ~success && initial_observer_enabled(config) && ~isfield(plan,'observationFallback')
        [state,plan] = initial_observation_fallback(state,plan,attempts);
        [state,success,attempts] = ipm.evolve.applyAutonomousMesh(state,plan);
    end
    if ~success
        detail = jsonencode(attempts);
        if isfield(plan,'observationFallback')
            evidence=plan.observationFallback;evidence.nativeAttempts=attempts;
            evidence.status='native_candidates_rejected';
            exception=MException('ipm:AutonomousMeshInitialization', ...
                'The observed initial candidates failed their unchanged native gates.');
            throw(addCause(exception,MException('ipm:InitialObservationEvidence','%s', ...
                initial_failure_message(state,evidence))));
        end
        error('ipm:AutonomousMeshInitialization', ...
            'Zero-time automatic mesh candidates failed: %s',detail);
    end
end
if ~ipm.evolve.isFinite(state)
    error('ipm:NonFiniteInitialState', ...
        'The resolved initial state contains a non-finite evolution field.');
end
end

function metadata = initial_remesh_warning( ...
        metadata,ops,flow,scale,code,message)
message = strrep(char(message),newline,' ');
row = struct('code',char(code),'metric','remesh_attempt', ...
    'value',NaN,'threshold',NaN,'message',message, ...
    'step',0,'canonicalTime',0,'physicalTime',0, ...
    'remeshCount',ops.remeshCount);
if isfield(metadata,'remeshWarnings') && ~isempty(metadata.remeshWarnings)
    metadata.remeshWarnings(end+1) = row;
else
    metadata.remeshWarnings = row;
end
metadata.quadrantRemeshRetry = ipm.remesh.recordLevelSetFailure( ...
    flow,ops,code,[],0,scale);
fprintf(2,['WARNING [initial-remesh/%s] remesh=%d: %s. ' ...
    'Initial grid retained; evolution continues.\n'], ...
    row.code,row.remeshCount,row.message);
end

function yes = initial_observer_enabled(config)
yes=isfield(config.remesh,'initialMeshObservationFallback') && ...
    config.remesh.initialMeshObservationFallback.enabled;
end

function [state,plan] = initial_observation_fallback(state,plan,attempts)
assert(state.step==0 && state.normalizedTime==0 && state.scale.canonicalTime==0 && ...
    state.scale.physicalTime==0 && state.ops.remeshCount==0, ...
    'ipm:InitialObservationContext','Observation fallback cannot replace an evolved field.');
original=struct('plan',plan,'attempts',attempts);
try
    [pairs,evidence]=ipm.remesh.planInitialAnalyticFallback(state.config, ...
        struct('x',state.ops.baseX,'y',state.ops.baseY),original);
catch exception
    cause=MException('ipm:InitialObservationOriginalFailure','%s',jsonencode(original));
    throw(addCause(exception,cause));
end
if isempty(pairs)
    exception=MException('ipm:InitialObservationCapacity', ...
        'The one registered initial observation fallback has no qualified candidates.');
    cause=MException('ipm:InitialObservationEvidence','%s', ...
        initial_failure_message(state,evidence));
    throw(addCause(exception,cause));
end
plan.candidates=pairs;plan.axisReport=evidence.axisReport;plan.stopReason='';
plan.observationFallback=evidence;
state.runMetadata.autonomousMesh.lastDecision=rmfield(plan, ...
    {'candidates','axisReport','observationFallback'});
end

function message=initial_failure_message(state,evidence)
% Version five keeps a bounded console diagnostic. The original analytic
% datum and policy reproduce the complete trial arrays without logging them.
if state.config.remesh.autonomousMesh.version~=5
    message=jsonencode(evidence);return
end
original=evidence.originalAttempt;
summary=struct('version',evidence.version,'status',evidence.status, ...
    'sourceNodeCount',evidence.sourceNodeCount, ...
    'observationNodeCount',evidence.observation.nodeCount, ...
    'originalFailureClass',evidence.originalFailureClass, ...
    'originalStopReason',original.plan.stopReason, ...
    'xAdmitted',sum([evidence.axisReport.xTrials.admissible]), ...
    'yAdmitted',sum([evidence.axisReport.yTrials.admissible]), ...
    'xRejections',rejection_counts(evidence.axisReport.xTrials), ...
    'yRejections',rejection_counts(evidence.axisReport.yTrials), ...
    'candidateAudits',struct([]));
for k=1:numel(evidence.candidateAudits)
    a=evidence.candidateAudits(k);
    row=struct('originalPairIndex',a.originalPairIndex, ...
        'geometryPassed',a.geometryPassed,'familyPassed',a.familyPassed, ...
        'actualCoreCells',a.features.actualCoreCells, ...
        'actualFrontCells',a.features.leftFrontCells);
    if isempty(summary.candidateAudits),summary.candidateAudits=row;
    else,summary.candidateAudits(end+1)=row;end
end
if isfield(evidence,'nativeAttempts')
    summary.nativeAttemptCount=numel(evidence.nativeAttempts);
end
message=jsonencode(summary);
end

function rows=rejection_counts(trials)
names={};
for k=1:numel(trials)
    names=[names,trials(k).reasons]; %#ok<AGROW>
end
uniqueNames=unique(names);rows=struct([]);
for k=1:numel(uniqueNames)
    row=struct('reason',uniqueNames{k},'count',sum(strcmp(names,uniqueNames{k})));
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
end
