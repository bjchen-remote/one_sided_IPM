function report(result)
%IPM.OUTPUT.REPORT Print only the terminal stop reason for verbose runs.

result = ipm.output.validate(result);
if ~result.config.output.verbose
    return;
end
fprintf('stop=%s tau=%.5f steps=%d\n', ...
    result.state.stopReason,result.state.canonicalTime,result.state.steps);
end
