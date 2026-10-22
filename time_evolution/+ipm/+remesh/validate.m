function validation = validate( ...
    metrics,ops,flow,config,proposal,context,retry,amplitudeScale)
%IPM.REMESH.VALIDATE Check core counts and physical peak preservation.

if nargin < 8
    amplitudeScale = 1;
end
remesh = config.remesh;

if retry == 0
    if proposal.xImproved
        xCoreFloor = 0.98*flow.coreGridPoints;
    else
        xCoreFloor = min( ...
            0.98*flow.coreGridPoints,context.xTargets(end));
    end
    if ~context.trackY
        yCoreFloor = -Inf;
    elseif proposal.yImproved
        yCoreFloor = 0.98*flow.verticalCoreGridPoints;
    else
        yCoreFloor = min(0.98*flow.verticalCoreGridPoints, ...
            context.yTargets(end));
    end
else
    xCoreFloor = 0.98*flow.coreGridPoints;
    yCoreFloor = 0.98*flow.verticalCoreGridPoints;
end
xCoreSafe = ~isfinite(flow.coreGridPoints) || ...
    flow.coreGridPoints <= 0 || metrics.xCore >= xCoreFloor;
yCoreSafe = ~context.trackY || ...
    ~isfinite(flow.verticalCoreGridPoints) || ...
    flow.verticalCoreGridPoints <= 0 || metrics.yCore >= yCoreFloor;
peakSafe = ~strcmp(ops.rescalingMode,'physical') || ...
    abs(metrics.peak-flow.trackedWallPeak) <= ...
    remesh.remeshPeakChangeTolerance*max(flow.trackedWallPeak,eps);
physicalRangeViolation = metrics.rangeViolation/amplitudeScale;
rangeSafe = physicalRangeViolation <= ...
    config.diagnostics.rangeStopTolerance;
validation = struct( ...
    'accepted',xCoreSafe && yCoreSafe && peakSafe && rangeSafe, ...
    'xCoreSafe',xCoreSafe,'yCoreSafe',yCoreSafe,'peakSafe',peakSafe, ...
    'rangeSafe',rangeSafe,'physicalRangeViolation',physicalRangeViolation, ...
    'xCoreFloor',xCoreFloor,'yCoreFloor',yCoreFloor);
end
