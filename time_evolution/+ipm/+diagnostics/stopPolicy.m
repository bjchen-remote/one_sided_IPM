function stopReason = stopPolicy( ...
        history,physicalGradientInf,flow,ops,diagnostics,remesh)
%IPM.DIAGNOSTICS.STOPPOLICY Reproduce the retired terminal classification.
%   The solver does not call this historical audit helper. Runtime quality
%   thresholds are emitted by postRemeshWarnings only after an accepted
%   replacement grid and never become a stop reason.

common = history.common;
mesh = history.mesh;
stopReason = '';
if physicalGradientInf >= diagnostics.gradientStop
    stopReason = 'gradient_threshold';
elseif common.physicalRangeViolation(end) >= diagnostics.rangeStopTolerance
    stopReason = 'maximum_principle_violation';
elseif (common.wallPeakLocalMaxima(end) > 1 && ...
        common.wallPeakTVRatio(end) > diagnostics.oscillationTVTolerance) || ...
        common.positiveWallNegativeRatio(end) > ...
        diagnostics.positiveWallNegativeTolerance
    stopReason = 'wall_profile_oscillation';
elseif mesh.maximumCellRatioX(end) > ...
        remesh.remeshMaximumCellRatio*(1+1e-10) || ...
        mesh.maximumCellRatioY(end) > ...
        remesh.remeshMaximumCellRatio*(1+1e-10)
    stopReason = 'grid_smoothness_failure';
elseif ipm.diagnostics.resolutionFailed(flow,ops)
    stopReason = 'grid_resolution_failure';
end
end
