function [targets,triggerCells]=directLevelSetTargets(ops)
%IPM.REMESH.DIRECTLEVELSETTARGETS Core goals and 80% trigger cell counts.
targets=[target_cells(ops,'adaptiveTargetLevelPoints',ops.nx), ...
    target_cells(ops,'adaptiveTargetVerticalLevelPoints',ops.ny)];
triggerCells=.8*targets;
end

function target=target_cells(ops,fieldName,nodeCount)
if isfield(ops.rescaling,fieldName)
    values=ops.rescaling.(fieldName);
    target=ceil(values(end));
else
    target=16;
end
target=min([32,max(12,target),max(8,floor((nodeCount-1)/3))]);
end
