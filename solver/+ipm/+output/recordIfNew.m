function [log,recorded] = recordIfNew(log,state)
%IPM.OUTPUT.RECORDIFNEW Append a state unless its solver time is already stored.

recorded = false;
if isfield(log,'history') && isfield(log.history,'common') && ...
        isfield(log.history.common,'t') && ...
        ~isempty(log.history.common.t)
    previousTime = log.history.common.t(end);
    tolerance = 10*eps(max([1,abs(previousTime), ...
        abs(state.normalizedTime)]));
    if abs(state.normalizedTime-previousTime) <= tolerance
        return;
    end
end
[log,~] = ipm.output.record(log,state);
recorded = true;
end
