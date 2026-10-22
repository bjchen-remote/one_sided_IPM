function [proposal,context] = ...
    propose(ops,flow,config,retry,context)
%IPM.REMESH.PROPOSE Propose one initial or damped remesh grid.

grid = config.grid;
remesh = config.remesh;
if nargin < 5 || retry == 0
    [proposal,context] = initial_proposal(ops,flow,grid,remesh);
    return;
end

retryFactor = 0.8^retry;
if context.initialXImproved
    [xNew,xAxisInfo] = ipm.remesh.axis(context.xReference, ...
        context.wallMetrics.bounds,context.peakX,context.xTargets, ...
        context.minimumSpacing(1),remesh,context.paired, ...
        context.initialXAmplitude*retryFactor);
    xAmplitude = xAxisInfo.amplitude;
else
    xNew = context.xOld;
    xAmplitude = 0;
    xAxisInfo = context.initialXAxisInfo;
end
if context.initialYImproved
    [yRow,yAxisInfo] = ipm.remesh.axis(context.yReference, ...
        context.verticalMetrics.bounds,0,context.yTargets, ...
        context.minimumSpacing(2),remesh,false, ...
        context.initialYAmplitude*retryFactor);
    yNew = yRow';
    yAmplitude = yAxisInfo.amplitude;
else
    yNew = context.yOld;
    yAmplitude = 0;
    yAxisInfo = context.initialYAxisInfo;
end
newPeakSpacing = local_spacing(xNew,context.peakX);
newWallSpacing = yNew(2)-yNew(1);
xImproved = context.oldPeakSpacing-newPeakSpacing >= ...
    remesh.remeshMinimumImprovement*context.oldPeakSpacing;
yImproved = context.oldWallSpacing-newWallSpacing >= ...
    remesh.remeshMinimumImprovement*context.oldWallSpacing;
proposal = make_proposal(xNew,yNew,xAxisInfo,yAxisInfo, ...
    xAmplitude,yAmplitude,newPeakSpacing,newWallSpacing, ...
    xImproved,yImproved);
end

function [proposal,context] = initial_proposal(ops,flow,grid,remesh)
xOld = ops.x;
yOld = ops.y;
if strcmpi(string(remesh.remeshReferenceMode),'current')
    xReference = xOld;
    yReference = yOld';
else
    xReference = ops.baseX;
    yReference = ops.baseY';
end
if isfinite(flow.omegaGaugePeakX)
    peakX = flow.omegaGaugePeakX;
else
    peakX = flow.trackedPeakX;
end
targetSpacing = grid.targetCenterSpacing(:)';
if isscalar(targetSpacing)
    targetSpacing = [targetSpacing,targetSpacing];
end
minimumSpacing = remesh.remeshMinimumSpacingFactor*targetSpacing;
levels = ops.rescaling.adaptiveLevels;
wallMetrics = ipm.diagnostics.peakResolution(max(flow.source(1,:),0),ops.x,levels);
xSafety = max([ops.rescaling.safetyPeakPoints / ...
    max(wallMetrics.areaPoints(1),1), ...
    ops.rescaling.safetyLevelPoints./max(wallMetrics.gridPoints,1)]);
xTargets = min(ops.rescaling.adaptiveTargetLevelPoints, ...
    max(ops.rescaling.safetyLevelPoints/remesh.remeshTargetSafety, ...
    1.15*wallMetrics.gridPoints));
xRequested = xSafety > remesh.remeshSafetyTrigger;
paired = strcmp(ops.symmetryMode,'double_odd_omega');
if xRequested
    [xNew,xAxisInfo] = ipm.remesh.axis(xReference, ...
        wallMetrics.bounds,peakX,xTargets,minimumSpacing(1),remesh,paired);
else
    xNew = xOld;
    xAxisInfo = unchanged_axis_info(xOld,wallMetrics.gridPoints,xTargets);
end
xAmplitude = xAxisInfo.amplitude;
xWindow = max(wallMetrics.widths);

trackY = isfinite(flow.trackedVerticalCoreWidth) && ...
    flow.trackedVerticalCoreWidth > 0;
if trackY
    verticalFeature = abs(flow.source(:,wallMetrics.peakIndex))';
    verticalMetrics = ipm.diagnostics.peakResolution(verticalFeature,ops.y',levels);
    ySafety = max([ops.rescaling.safetyVerticalPeakPoints / ...
        max(verticalMetrics.areaPoints(1),1), ...
        ops.rescaling.safetyVerticalLevelPoints ./ ...
        max(verticalMetrics.gridPoints,1)]);
    yTargets = min(ops.rescaling.adaptiveTargetVerticalLevelPoints, ...
        max(ops.rescaling.safetyVerticalLevelPoints / ...
        remesh.remeshTargetSafety,1.15*verticalMetrics.gridPoints));
    yRequested = ySafety > remesh.remeshSafetyTrigger;
    if yRequested
        [yRow,yAxisInfo] = ipm.remesh.axis(yReference, ...
            verticalMetrics.bounds,0,yTargets,minimumSpacing(2),remesh,false);
        yNew = yRow';
    else
        yNew = yOld;
        yAxisInfo = unchanged_axis_info( ...
            yOld',verticalMetrics.gridPoints,yTargets);
    end
    yAmplitude = yAxisInfo.amplitude;
    yWindow = max(verticalMetrics.widths);
else
    verticalMetrics = [];
    yTargets = NaN(size(levels));
    yNew = yOld;
    yAmplitude = 0;
    yWindow = NaN;
    yRequested = false;
    yAxisInfo = unchanged_axis_info( ...
        yOld',NaN(size(levels)),NaN(size(levels)));
end

oldPeakSpacing = local_spacing(xOld,peakX);
newPeakSpacing = local_spacing(xNew,peakX);
oldWallSpacing = yOld(2)-yOld(1);
newWallSpacing = yNew(2)-yNew(1);
xImproved = xRequested && oldPeakSpacing-newPeakSpacing >= ...
    remesh.remeshMinimumImprovement*oldPeakSpacing;
yImproved = yRequested && oldWallSpacing-newWallSpacing >= ...
    remesh.remeshMinimumImprovement*oldWallSpacing;
if ~xImproved
    xNew = xOld;
    newPeakSpacing = oldPeakSpacing;
end
if ~yImproved
    yNew = yOld;
    newWallSpacing = oldWallSpacing;
end
proposal = make_proposal(xNew,yNew,xAxisInfo,yAxisInfo, ...
    xAmplitude,yAmplitude,newPeakSpacing,newWallSpacing, ...
    xImproved,yImproved);
context = struct('xOld',xOld,'yOld',yOld, ...
    'xReference',xReference,'yReference',yReference, ...
    'peakX',peakX,'minimumSpacing',minimumSpacing, ...
    'wallMetrics',wallMetrics,'xTargets',xTargets,'paired',paired, ...
    'xWindow',xWindow,'trackY',trackY, ...
    'verticalMetrics',verticalMetrics,'yTargets',yTargets, ...
    'yWindow',yWindow,'oldPeakSpacing',oldPeakSpacing, ...
    'oldWallSpacing',oldWallSpacing, ...
    'initialXImproved',xImproved,'initialYImproved',yImproved, ...
    'initialXAmplitude',xAmplitude,'initialYAmplitude',yAmplitude, ...
    'initialXAxisInfo',xAxisInfo,'initialYAxisInfo',yAxisInfo);
end

function proposal = make_proposal(xNew,yNew,xAxisInfo,yAxisInfo, ...
        xAmplitude,yAmplitude,newPeakSpacing,newWallSpacing, ...
        xImproved,yImproved)
proposal = struct('x',xNew,'y',yNew, ...
    'xAxisInfo',xAxisInfo,'yAxisInfo',yAxisInfo, ...
    'xAmplitude',xAmplitude,'yAmplitude',yAmplitude, ...
    'newPeakSpacing',newPeakSpacing,'newWallSpacing',newWallSpacing, ...
    'xImproved',xImproved,'yImproved',yImproved, ...
    'anyImproved',xImproved || yImproved);
end

function spacing = local_spacing(axis,target)
[~,index] = min(abs(axis-target));
left = max(index-1,1);
right = min(index+1,numel(axis));
spacing = (axis(right)-axis(left))/(right-left);
end

function info = unchanged_axis_info(axis,counts,targets)
faces = diff(axis);
if numel(faces) > 1
    maximumCellRatio = max(max(faces(2:end)./faces(1:end-1), ...
        faces(1:end-1)./faces(2:end)));
else
    maximumCellRatio = 1;
end
info = struct('amplitude',0,'minimumSpacing',min(faces), ...
    'maximumCellRatio',maximumCellRatio,'counts',counts, ...
    'targetCounts',targets,'targetReached',true);
end
