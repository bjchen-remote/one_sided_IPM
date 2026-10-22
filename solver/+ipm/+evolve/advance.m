function [state,stopReason,meshPlan] = advance(state)
%IPM.EVOLVE.ADVANCE Select dt, advance the configured RK scheme, and remesh.

meshPlan = [];
[dt,stopReason,timeStep] = ipm.evolve.selectTimestep( ...
    state.flow,state.scale,state.ops,state.config.time);
if ~isempty(stopReason)
    return;
end

% Commit a step only after its complete state is known to be finite.  This
% keeps the last accepted field available for diagnostics and persistence
% when an attempted Runge--Kutta stage overflows.
acceptedState = state;
initialRhsCache = [];
if isfield(state,'rhsCache')
    initialRhsCache = state.rhsCache;
end
switch state.config.time.timeIntegrator
    case 'ssprk3'
        [state.rho,state.flow,state.scale,state.rhsCache] = ...
            ipm.evolve.stepSsprk3( ...
            state.rho,dt,state.ops,state.scale,initialRhsCache);
    case 'ssprk54'
        [state.rho,state.flow,state.scale,state.rhsCache] = ...
            ipm.evolve.stepSsprk54( ...
            state.rho,dt,state.ops,state.scale,initialRhsCache);
    case 'rk6'
        [state.rho,state.flow,state.scale,state.rhsCache] = ...
            ipm.evolve.stepRk6( ...
            state.rho,dt,state.ops,state.scale,initialRhsCache);
    otherwise
        error('ipm:TimeIntegrator','Unknown time integrator ''%s''.', ...
            state.config.time.timeIntegrator);
end
% normalizedTime is retained as a result/checkpoint compatibility alias.
state.normalizedTime = state.scale.canonicalTime;
state.step = state.step+1;

if ~ipm.evolve.isFinite(state)
    state = acceptedState;
    stopReason = 'non_finite_solution';
    return;
end

% Isolated masked-CFL continuation guard. The full RHS is still evaluated
% everywhere; this only stops the experiment if excluded-region dynamics or
% one-step norms cease to be small relative to the resolved field.
if ~isempty(getenv('IPM_SOURCE_CFL_REL'))
    sourceMagnitude = abs(state.flow.source);
    sourceMaximum = max(sourceMagnitude,[],'all');
    previousSourceMaximum = max(abs(acceptedState.flow.source),[],'all');
    rhoMaximum = max(abs(state.rho),[],'all');
    previousRhoMaximum = max(abs(acceptedState.rho),[],'all');
    relativeThreshold = str2double(getenv('IPM_SOURCE_CFL_REL'));
    previousSourceMagnitude = abs(acceptedState.flow.source);
    previousActive = conv2(double(previousSourceMagnitude >= ...
        relativeThreshold*previousSourceMaximum), ...
        ones(7,7),'same') > 0;
    active = conv2(double(sourceMagnitude >= ...
        relativeThreshold*sourceMaximum), ...
        ones(7,7),'same') > 0;
    rhsMagnitude = abs(state.rhsCache.rhoRate);
    rhsMaximum = max(rhsMagnitude,[],'all');
    if any(~active,'all')
        omittedRhsRatio = max(rhsMagnitude(~active))/max(rhsMaximum,realmin);
    else
        omittedRhsRatio = 0;
    end
    previousOmitted = find(~previousActive);
    if isempty(previousOmitted)
        previousOmittedSourceRatio = 0;
        previousOmittedRhsRatio = 0;
        omittedSourceX = NaN;
        omittedSourceY = NaN;
    else
        [previousOmittedSource,peakLocalIndex] = ...
            max(sourceMagnitude(previousOmitted));
        previousOmittedSourceRatio = previousOmittedSource/ ...
            max(sourceMaximum,realmin);
        previousOmittedRhsRatio = max(rhsMagnitude(previousOmitted))/ ...
            max(rhsMaximum,realmin);
        [peakRow,peakColumn] = ind2sub(size(sourceMagnitude), ...
            previousOmitted(peakLocalIndex));
        omittedSourceX = state.ops.x(peakColumn);
        omittedSourceY = state.ops.y(peakRow);
    end
    farCore = abs(state.ops.X) <= 1 & state.ops.Y >= 1e5;
    farCoreSourceRatio = max(sourceMagnitude(farCore))/ ...
        max(sourceMaximum,realmin);
    state.runMetadata.experimentalCflGuard = struct( ...
        'previousOmittedSourceRatio',previousOmittedSourceRatio, ...
        'previousOmittedRhsRatio',previousOmittedRhsRatio, ...
        'previousOmittedPeakX',omittedSourceX, ...
        'previousOmittedPeakY',omittedSourceY, ...
        'currentOmittedRhsRatio',omittedRhsRatio, ...
        'farCoreSourceRatio',farCoreSourceRatio);
    unsafe = omittedRhsRatio > 1e-6 || ...
        previousOmittedSourceRatio > 1e-6 || ...
        previousOmittedRhsRatio > 1e-6 || ...
        farCoreSourceRatio > 1e-6 || ...
        sourceMaximum > 1.02*previousSourceMaximum || ...
        rhoMaximum > 1.02*previousRhoMaximum;
    if unsafe
        fprintf(2,['MASK_GUARD tau=%.12g step=%d ' ...
            'omittedRhsRatio=%.3e previousOmittedSourceRatio=%.3e ' ...
            'previousOmittedRhsRatio=%.3e atX=%.9g atY=%.9g ' ...
            'farCoreSourceRatio=%.3e sourceGrowth=%.6g ' ...
            'rhoGrowth=%.6g\n'], ...
            state.scale.canonicalTime,state.step,omittedRhsRatio, ...
            previousOmittedSourceRatio,previousOmittedRhsRatio, ...
            omittedSourceX,omittedSourceY,farCoreSourceRatio, ...
            sourceMaximum/max(previousSourceMaximum,realmin), ...
            rhoMaximum/max(previousRhoMaximum,realmin));
        state = acceptedState;
        stopReason = 'experimental_cfl_tail_guard';
        return;
    end
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
