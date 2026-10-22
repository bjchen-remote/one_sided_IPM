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
else
    requested = isfinite(flow.safetyFactor) && ...
        flow.safetyFactor > remesh.remeshSafetyTrigger;
end
if ~requested
    return;
end

[rho,ops,info] = ipm.remesh.adapt( ...
    rho,ops,flow,config,exp(scale.logC_omega));
if ~info.applied
    return;
end
[rhoRate,flow] = ipm.evolve.flow(rho,ops,scale);
if output.verbose
    fprintf(['  adaptive remesh %d at physical t=%.6g: ' ...
        'peak dx %.3e -> %.3e, wall dy %.3e -> %.3e\n'], ...
        info.count,scale.physicalTime,info.oldPeakSpacing, ...
        info.newPeakSpacing,info.oldWallSpacing,info.newWallSpacing);
end
end
