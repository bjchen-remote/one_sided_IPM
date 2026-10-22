function [rho,ops,flow,info,rhoRate] = ...
    remeshIfNeeded(rho,ops,flow,scale,config)
%IPM.EVOLVE.REMESHIFNEEDED Apply one accepted adaptive-grid proposal when due.

remesh = config.remesh;
output = config.output;
info = struct('applied',false);
rhoRate = [];
if ~remesh.adaptiveRemesh || ops.remeshCount>=remesh.maxRemeshes
    return
end
if ipm.mesh.isQuadrant(ops)
    % The positive-only planner uses only its 90% core level-set counts.
    % Reuse flow's counts only when its highest tracked level is 90%.
    % Other user-selected levels must be measured by the direct planner.
    core=[flow.coreGridPoints,flow.verticalCoreGridPoints];
    requested=all(isfinite(core));
    if requested && abs(ops.rescaling.adaptiveLevels(end)-.9)<=100*eps(1)
        [~,triggerCells]=ipm.remesh.directLevelSetTargets(ops);
        requested=any(core<triggerCells);
    end
elseif strcmp(ops.symmetryMode,'half_plane') && ...
        strcmp(ops.dynamicScaleGeometry,'isotropic')
    % The existing remesh transaction owns the geometric-area trigger and
    % its accepted reference. Receive its ops even when no grid is applied.
    requested = true;
else
    requested = isfinite(flow.safetyFactor) && ...
        flow.safetyFactor > remesh.remeshSafetyTrigger;
end
if ~requested
    return;
end

[rho,ops,info,rhoRate,validatedFlow] = ipm.remesh.adapt( ...
    rho,ops,flow,config,exp(scale.logC_omega));
if ~info.applied
    return;
end
if ~isempty(rhoRate) && ~isempty(validatedFlow)
    % The isotropic area transaction already audited this exact accepted
    % field and grid. Reuse its RHS and root instead of solving twice.
    flow = validatedFlow;
else
    [rhoRate,flow] = ipm.evolve.flow(rho,ops,scale);
end
if strcmp(ops.dynamicScaleGeometry,'isotropic') && ...
        ops.rescaling.enabled && ...
        strcmp(ops.rescaling.lengthGauge,'wall_gradient_integral')
    ops.rescaling.gaugeRateGuess = [flow.c_l,flow.c_r];
end
if isfield(flow,'peakTrackingAnchor')
    ops.rescaling.trackedPeakAnchor = flow.peakTrackingAnchor;
    ops.rescaling.peakEvaluationTime = flow.peakTrackingAnchor.tau;
end
if output.verbose
    detail = '';
    if isfield(info,'damping') && isfield(info,'targetReached')
        detail = sprintf(' damping=%.5g targetReached=%d', ...
            info.damping,info.targetReached);
    end
    fprintf(['  adaptive remesh %d at physical t=%.6g: ' ...
        'peak dx %.3e -> %.3e, wall dy %.3e -> %.3e%s\n'], ...
        info.count,scale.physicalTime,info.oldPeakSpacing, ...
        info.newPeakSpacing,info.oldWallSpacing,info.newWallSpacing,detail);
end
end
