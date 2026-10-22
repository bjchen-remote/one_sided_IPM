function [retry,evidence] = levelSetRetryDecision(flow,ops,memory)
%IPM.REMESH.LEVELSETRETRYDECISION Retry only after resolved feature motion.
%   The comparison uses the existing 90% wall-core center/width and vertical
%   core width. One target-cell change is the sole retry threshold; elapsed
%   steps and time do not enter the decision.

if nargin < 3
    memory = [];
end
[targets,~] = ipm.remesh.directLevelSetTargets(ops);
requiredChange = 1/min(targets);
current = [flow.trackedPeakX,flow.trackedWallCoreWidth, ...
    flow.trackedVerticalCoreWidth];
currentValid = all(isfinite(current)) && all(current(2:3) > 0);

reference = [NaN,NaN,NaN];
if isstruct(memory) && isscalar(memory) && ...
        isfield(memory,'reference') && numel(memory.reference) == 3
    reference = memory.reference(:)';
end
referenceValid = all(isfinite(reference)) && all(reference(2:3) > 0);
relativeChange = [NaN,NaN,NaN];
if currentValid && referenceValid
    horizontalScale = max([current(2),reference(2),realmin]);
    relativeChange = [abs(current(1)-reference(1))/horizontalScale, ...
        abs(log(current(2)/reference(2))), ...
        abs(log(current(3)/reference(3)))];
    retry = max(relativeChange) >= requiredChange;
elseif currentValid
    retry = true;
else
    retry = false;
end
if isempty(memory)
    retry = true;
end

evidence = struct('current',current,'reference',reference, ...
    'relativeChange',relativeChange,'maximumRelativeChange', ...
    max(relativeChange,[],'all','omitnan'), ...
    'requiredRelativeChange',requiredChange, ...
    'targetCells',targets,'currentValid',currentValid, ...
    'referenceValid',referenceValid,'stepOrTimeUsed',false);
end
