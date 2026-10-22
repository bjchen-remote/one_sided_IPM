function [state,cursor,fileName] = maybeCheckpoint(state,log,cursor,force)
%IPM.OUTPUT.MAYBECHECKPOINT Persist a due trusted recorded state.

if nargin < 4
    force = false;
end
fileName = "";
output = state.config.output;
if ~isfield(output,'checkpoint') || ~output.checkpoint.enabled
    return;
end
policy = output.checkpoint;
tau = state.scale.canonicalTime;
tolerance = 10*eps(max(1,abs(tau)));
due = tau+tolerance >= cursor.nextCheckpoint;
exitDue = force && policy.atExit;
if ~(due || exitDue) || cursor.lastCheckpointStep == state.step
    return;
end
trusted = ipm.output.continuousTrustedPrefix( ...
    ipm.output.trustedMask(log.history,state.config));
if ~trusted(end)
    return;
end
if due && isfinite(policy.every)
    while tau+tolerance >= cursor.nextCheckpoint
        cursor.nextCheckpoint = cursor.nextCheckpoint+policy.every;
    end
end
cursor.lastCheckpointStep = state.step;
cursor.lastCheckpointCanonicalTime = tau;
checkpoint = ipm.output.makeCheckpoint(state,log,cursor);
[~,fileName] = ipm.output.writeCheckpoint(checkpoint,policy.file);
state.runMetadata.latestCheckpointFile = string(fileName);
end
