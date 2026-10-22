function [state,stopReason,meshPlan] = advance(state)
%IPM.EVOLVE.ADVANCE Select dt, advance the configured RK scheme, and remesh.

meshPlan = [];
[dt,stopReason,timeStep] = ipm.evolve.selectTimestep( ...
    state.flow,state.scale,state.ops,state.config.time);
if ~isempty(stopReason)
    return;
end
recoveringGauge = false;
previousDt = NaN;
if isfield(state,'timeStep') && ...
        isfield(state.timeStep,'gaugeRecovery') && state.timeStep.gaugeRecovery && ...
        isfinite(state.timeStep.dt) && state.timeStep.dt > 0
    previousDt = state.timeStep.dt;
elseif isfield(state.runMetadata,'gaugeStepRecovery')
    saved = state.runMetadata.gaugeStepRecovery;
    if isstruct(saved) && isscalar(saved) && ...
            all(isfield(saved,{'dt','active','lastCanonicalTime'})) && ...
            islogical(saved.active) && isscalar(saved.active) && saved.active && ...
            isnumeric(saved.dt) && isreal(saved.dt) && isscalar(saved.dt) && ...
            isfinite(saved.dt) && saved.dt > 0 && ...
            isequal(saved.lastCanonicalTime,state.scale.canonicalTime)
        previousDt = saved.dt;
    end
end
if isfinite(previousDt)
    % After a rejected stage, regrow the accepted dt gradually. Repeating
    % the same failed large step on every call wastes complete Poisson/RK
    % evaluations. This ceiling never exceeds the newly computed CFL limit.
    recoveryLimit = 1.25*previousDt;
    if recoveryLimit < dt
        dt = recoveryLimit;
        timeStep.dt = dt;
        timeStep.realizedCfl = dt*timeStep.rate;
        timeStep.activeLimiter = 'wall_integral_gauge_recovery';
        recoveringGauge = true;
    end
end

% Commit a step only after its complete state is known to be finite.  This
% keeps the last accepted field available for diagnostics and persistence
% when an attempted Runge--Kutta stage overflows.
acceptedState = state;
initialRhsCache = [];
if isfield(state,'rhsCache')
    initialRhsCache = state.rhsCache;
end
initialDt = dt;
retryLimit=12;
implicit=any(strcmp(state.config.time.timeIntegrator, ...
    {'sdirk3_alexander_practical_v1','sdirk3_alexander_practical_v2','sdirk3_alexander_practical_v3','sdirk2_alexander_practical_v1','sdirk2_alexander_practical_v2'}));
if implicit,retryLimit=state.ops.timeIntegrationContract.maximumStepRetries;end
for gaugeRetry = 0:retryLimit
    try
        switch state.config.time.timeIntegrator
            case 'ssprk3'
                [state.rho,state.flow,state.scale,state.rhsCache] = ...
                    ipm.evolve.stepSsprk3( ...
                    acceptedState.rho,dt,acceptedState.ops, ...
                    acceptedState.scale,initialRhsCache);
            case 'ssprk54'
                [state.rho,state.flow,state.scale,state.rhsCache] = ...
                    ipm.evolve.stepSsprk54( ...
                    acceptedState.rho,dt,acceptedState.ops, ...
                    acceptedState.scale,initialRhsCache);
            case 'rk6'
                [state.rho,state.flow,state.scale,state.rhsCache] = ...
                    ipm.evolve.stepRk6( ...
                    acceptedState.rho,dt,acceptedState.ops, ...
                    acceptedState.scale,initialRhsCache);
            case {'sdirk3_alexander_practical_v1','sdirk3_alexander_practical_v2','sdirk3_alexander_practical_v3','sdirk2_alexander_practical_v1','sdirk2_alexander_practical_v2'}
                [state.rho,state.flow,state.scale,state.rhsCache,implicitEvaluation] = ...
                    ipm.evolve.stepRk(acceptedState.rho,dt,acceptedState.ops, ...
                    state.config.time.timeIntegrator,acceptedState.scale,initialRhsCache);
                if state.scale.physicalTime<=acceptedState.scale.physicalTime || ...
                        (isfinite(state.config.time.physicalFinalTime) && ...
                        state.scale.physicalTime>state.config.time.physicalFinalTime+ ...
                        64*eps(max(1,state.config.time.physicalFinalTime)))
                    error('ipm:ImplicitStageClock','The implicit step must advance within its physical horizon.');
                end
            otherwise
                error('ipm:TimeIntegrator','Unknown time integrator ''%s''.', ...
                    state.config.time.timeIntegrator);
        end
        break
    catch exception
        % CFL controls transport, but an intermediate field can still leave
        % the integral gauge's solvable neighbourhood. Retry the entire RK
        % step from the same accepted field with unchanged gauge tolerances.
        state = acceptedState;
        recoverable = state.ops.rescaling.enabled && ...
            strcmp(state.ops.rescaling.lengthGauge,'wall_gradient_integral') && ...
            any(strcmp(exception.identifier, ...
            {'ipm:UnresolvedWallIntegralGauge','ipm:DegenerateWallIntegralGauge', ...
             'ipm:TrackedPeakBranchLost'}));
        if implicit && any(strcmp(exception.identifier, ...
                {'ipm:ImplicitStageConvergence','ipm:ImplicitStageClock'}))
            recoverable=true;
        end
        nextDt = dt/2;
        if ~recoverable || gaugeRetry == retryLimit || ...
                nextDt < state.config.time.minDt || ...
                state.scale.canonicalTime+nextDt == state.scale.canonicalTime
            rethrow(exception);
        end
        if state.config.output.verbose
            fprintf(['  STEP_GAUGE_RETRY step=%d tau=%.8g retry=%d ' ...
                'dt=%.6e -> %.6e reason=%s\n'], ...
                state.step+1,state.scale.canonicalTime,gaugeRetry+1, ...
                dt,nextDt,exception.identifier);
        end
        dt = nextDt;
        timeStep.dt = dt;
        timeStep.realizedCfl = dt*timeStep.rate;
        timeStep.activeLimiter = 'wall_integral_gauge_retry';
        timeStep.initialDt = initialDt;
        timeStep.gaugeRetries = gaugeRetry+1;
    end
end
timeStep.gaugeRecovery = recoveringGauge || gaugeRetry > 0;
if implicit
    if any(strcmp(state.config.time.timeIntegrator,{'sdirk3_alexander_practical_v2','sdirk3_alexander_practical_v3','sdirk2_alexander_practical_v1','sdirk2_alexander_practical_v2'}))
        initialCfl=dt*timeStep.rate;
        stageCfl=[implicitEvaluation.implicit.stages.fullTransportCflDiagnostic];
        terminalCfl=dt*state.flow.generatorCflRate;
        timeStep.initialCfl=initialCfl;
        timeStep.actualStageMaximumCfl=max(stageCfl);
        timeStep.terminalCfl=terminalCfl;
        timeStep.realizedCfl=max([initialCfl,stageCfl,terminalCfl]);
        timeStep.realizedCflScope='maximum_initial_converged_implicit_stages_and_terminal_full_rates';
        implicitEvaluation.initialCfl=initialCfl;
        implicitEvaluation.actualStageMaximumCfl=timeStep.actualStageMaximumCfl;
        implicitEvaluation.terminalCfl=terminalCfl;
        implicitEvaluation.reportedRealizedCfl=timeStep.realizedCfl;
    end
    timeStep.implicit=implicitEvaluation;
    state.runMetadata.lastImplicitStep=implicitEvaluation;
end
% normalizedTime is retained as a result/checkpoint compatibility alias.
state.normalizedTime = state.scale.canonicalTime;
state.step = state.step+1;
if isfield(state.flow,'peakTrackingAnchor')
    state.ops.rescaling.trackedPeakAnchor = state.flow.peakTrackingAnchor;
    state.ops.rescaling.peakEvaluationTime = state.scale.canonicalTime;
end
if state.ops.rescaling.enabled && ...
        strcmp(state.ops.dynamicScaleGeometry,'isotropic') && ...
        strcmp(state.ops.rescaling.lengthGauge,'wall_gradient_integral')
    state.ops.rescaling.gaugeRateGuess = [state.flow.c_l,state.flow.c_r];
end

if ~ipm.evolve.isFinite(state)
    state = acceptedState;
    stopReason = 'non_finite_solution';
    return;
end
% Checkpoints omit the derived RHS cache/timeStep. Preserve the accepted
% controller ceiling in signed run metadata so segmented runs do not keep
% retrying a large dt that was already rejected.
state.runMetadata.gaugeStepRecovery = struct('dt',dt, ...
    'active',logical(timeStep.gaugeRecovery), ...
    'lastCanonicalTime',state.scale.canonicalTime);

if isfield(state.ops,'spatialOperators')
    % The RK, nonlinear stage solve, retry and recovery above are shared.
    % Mesh geometry owns proposal/transfer/validation and preserves the
    % accepted time endpoint when a candidate mesh is rejected.
    if isfield(state.ops.spatialOperators,'remesh')
        [state,remeshInfo]=state.ops.spatialOperators.remesh(state);
        state.runMetadata.lastRemesh=remeshInfo;
        if remeshInfo.applied
            if ~isfield(state,'rhsCache')||~isstruct(state.rhsCache)||isempty(fieldnames(state.rhsCache))
                [rate,state.flow]=ipm.evolve.flow(state.rho,state.ops,state.scale);
                state.rhsCache=ipm.evolve.makeRhsCache(state.rho,rate,state.flow,state.ops,state.scale);
            end
        end
    end
    state.timeStep=timeStep;
    if ~ipm.evolve.isFinite(state)
        state=acceptedState;stopReason='non_finite_solution';
    end
    return;
end

if isfield(state.config.remesh,'autonomousMesh') && ...
        state.config.remesh.autonomousMesh.enabled
    % Return the pure plan to solve, which owns the last external reference
    % to the old factor. No new LU may be constructed in this call scope.
    try
        [state,meshPlan] = ipm.evolve.planAutonomousMesh(state,false);
    catch exception
        failedStep = state.step;
        failureClock = [state.scale.canonicalTime,state.scale.physicalTime,state.normalizedTime];
        state = acceptedState;
        state.runMetadata.autonomousMesh.lastFailure = struct( ...
            'kind','planning_exception','attemptedStep',failedStep, ...
            'identifier',exception.identifier,'message',exception.message);
        if any(state.config.remesh.autonomousMesh.version == [2,3,4,5])
            state.runMetadata.autonomousMesh.lastFailure.sourceCanonicalTime=failureClock(1);
            state.runMetadata.autonomousMesh.lastFailure.sourcePhysicalTime=failureClock(2);
            state.runMetadata.autonomousMesh.lastFailure.sourceNormalizedTime=failureClock(3);
        end
        stopReason = 'autonomous_mesh_error';
        meshPlan = [];
        return;
    end
    if ~isempty(meshPlan.stopReason)
        failedPlan = rmfield(meshPlan,{'candidates','axisReport'});
        if any(state.config.remesh.autonomousMesh.version == [2,3,4,5])
            failedPlan.axisReport=meshPlan.axisReport;
        end
        state = acceptedState;
        state.runMetadata.autonomousMesh.lastFailure = failedPlan;
        stopReason = meshPlan.stopReason;
        meshPlan = [];
        return;
    end
    policy = state.config.remesh.autonomousMesh;
    if any([state.flow.coreGridPoints,state.flow.verticalCoreGridPoints] < ...
            policy.endpointCoreFloor)
        failureDecision=rmfield(meshPlan,{'candidates','axisReport'});
        state = acceptedState;
        if any(policy.version == [2,3,4,5])
            state.runMetadata.autonomousMesh.lastFailure=failureDecision;
            state.runMetadata.autonomousMesh.lastFailure.stopReason='autonomous_mesh_resolution_failure';
        end
        stopReason = 'autonomous_mesh_resolution_failure';
        meshPlan = [];
        return;
    end
    state.timeStep = timeStep;
    if ~meshPlan.requested
        meshPlan = [];
    end
    return;
end

retryMemory = [];
if ipm.mesh.isQuadrant(state.ops) && ...
        isfield(state.runMetadata,'quadrantRemeshRetry')
    retryMemory = state.runMetadata.quadrantRemeshRetry;
    retry = ipm.remesh.levelSetRetryDecision( ...
        state.flow,state.ops,retryMemory);
    if ~retry
        state.timeStep = timeStep;
        return;
    end
end
try
    [state.rho,state.ops,state.flow,remeshInfo,remeshedRhoRate] = ...
        ipm.evolve.remeshIfNeeded( ...
        state.rho,state.ops,state.flow,state.scale,state.config);
catch exception
    state = ipm.evolve.publishRemeshAttemptWarning(state, ...
        'remesh_exception',sprintf('%s: %s', ...
        exception.identifier,exception.message));
    state.runMetadata.quadrantRemeshRetry = ...
        ipm.remesh.recordLevelSetFailure(state.flow,state.ops, ...
        'remesh_exception',retryMemory,state.step,state.scale);
    state.timeStep = timeStep;
    return;
end
if remeshInfo.applied
    state.rhsCache = ipm.evolve.makeRhsCache( ...
        state.rho,remeshedRhoRate,state.flow,state.ops,state.scale);
end
if ~remeshInfo.applied && isfield(remeshInfo,'status') && ...
        ~strcmp(remeshInfo.status,'not_triggered')
    state = ipm.evolve.publishRemeshAttemptWarning(state, ...
        'remesh_candidate_rejected',sprintf( ...
        'direct level-set remesh returned %s',remeshInfo.status));
    state.runMetadata.quadrantRemeshRetry = ...
        ipm.remesh.recordLevelSetFailure(state.flow,state.ops, ...
        remeshInfo.status,retryMemory,state.step,state.scale);
elseif ipm.mesh.isQuadrant(state.ops) && ...
        isfield(state.runMetadata,'quadrantRemeshRetry')
    state.runMetadata = rmfield(state.runMetadata,'quadrantRemeshRetry');
end
if ~ipm.evolve.isFinite(state)
    state = acceptedState;
    stopReason = 'non_finite_solution';
    return;
end
if remeshInfo.applied
    if isfield(state.runMetadata,'quadrantRemeshRetry')
        state.runMetadata = rmfield(state.runMetadata, ...
            'quadrantRemeshRetry');
    end
    state = ipm.evolve.publishRemeshWarnings(state);
end
state.timeStep = timeStep;
end
