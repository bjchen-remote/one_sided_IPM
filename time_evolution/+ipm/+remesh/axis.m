function [newAxis,info] = axis( ...
        baseAxis,featureBounds,peakCenter,targetCounts,minimumSpacing, ...
        remesh,paired,amplitudeCap)
%IPM.REMESH.AXIS Build a smooth, target-driven nonuniform axis.
%   Nested connected-level intervals supply the monitor. For a paired odd
%   feature, a graded bridge prevents a coarse hole between the origin and
%   the approaching positive peak. The weakest amplitude meeting the target
%   counts is used, subject to spacing and adjacent-cell-ratio constraints.
%   No global max(dx)/min(dx) bound is imposed. The optional lattice profile
%   adapts the flat-core/logistic-transition idea in the legacy lat* maps.

baseAxis = baseAxis(:)';
targetCounts = targetCounts(:)';
if size(featureBounds,1) ~= numel(targetCounts) || ...
        size(featureBounds,2) ~= 2
    error('ipm:AdaptiveAxisBounds', ...
        'Feature bounds and target counts must have compatible sizes.');
end
reference = linspace(0,1,numel(baseAxis));
if nargin < 8
    amplitudeCap = 1e10;
end

[baseCandidate,baseStats] = candidate_axis(0);
if all(baseStats.counts >= targetCounts)
    newAxis = baseCandidate;
    info = make_info(0,baseStats,true);
    return;
end

best = struct('amplitude',0,'axis',baseCandidate,'stats',baseStats, ...
    'score',min(baseStats.counts./targetCounts));
previousAmplitude = 0;
lastValidAmplitude = 0;
lastValidAxis = baseCandidate;
lastValidStats = baseStats;
trialAmplitude = min(max(remesh.remeshConcentration,1),amplitudeCap);
targetReached = false;
while trialAmplitude > previousAmplitude
    [trialAxis,trialStats] = candidate_axis(trialAmplitude);
    if trialStats.valid
        lastValidAmplitude = trialAmplitude;
        lastValidAxis = trialAxis;
        lastValidStats = trialStats;
        trialScore = min(trialStats.counts./targetCounts);
        if trialScore > best.score
            best = struct('amplitude',trialAmplitude,'axis',trialAxis, ...
                'stats',trialStats,'score',trialScore);
        end
        if all(trialStats.counts >= targetCounts)
            targetReached = true;
            break;
        end
    elseif all(trialStats.counts >= targetCounts)
        % Cell-ratio admissibility is not monotone in amplitude: a
        % partially formed transition can be rougher than the stronger,
        % fully resolved monitor. Locate the first count-satisfying axis
        % inside this bracket and test that axis directly.
        bracketLower = previousAmplitude;
        bracketUpper = trialAmplitude;
        bracketAxis = trialAxis;
        bracketStats = trialStats;
        for iteration = 1:32
            bracketAmplitude = 0.5*(bracketLower+bracketUpper);
            [bracketCandidate,bracketCandidateStats] = ...
                candidate_axis(bracketAmplitude);
            if all(bracketCandidateStats.counts >= targetCounts)
                bracketUpper = bracketAmplitude;
                bracketAxis = bracketCandidate;
                bracketStats = bracketCandidateStats;
            else
                bracketLower = bracketAmplitude;
            end
        end
        if bracketStats.valid
            trialAmplitude = bracketUpper;
            trialAxis = bracketAxis;
            trialStats = bracketStats;
            targetReached = true;
            previousAmplitude = bracketLower;
            break;
        else
            % With a strict adjacent-ratio cap, even the first
            % count-satisfying amplitude can already be inadmissible. Keep
            % the strongest valid sub-target deformation below it instead
            % of falling back to the unchanged grid.
            best = retain_strongest_valid(best,lastValidAmplitude, ...
                bracketUpper,lastValidAxis,lastValidStats);
        end
    elseif trialStats.minimumSpacing < minimumSpacing
        % Increasing monitor strength cannot restore the spacing floor.
        break;
    else
        % A ratio-invalid trial can still have a useful admissible amplitude
        % below it even when no single remesh reaches every target. Retain the
        % strongest valid incremental improvement; current-reference mode can
        % then resolve nested scales through several accepted small changes.
        best = retain_strongest_valid(best,lastValidAmplitude, ...
            trialAmplitude,lastValidAxis,lastValidStats);
    end
    previousAmplitude = trialAmplitude;
    if trialAmplitude >= amplitudeCap
        break;
    end
    trialAmplitude = min(2*trialAmplitude,amplitudeCap);
end

if targetReached
    % Find the least deformation that supplies the requested reserve.
    targetUpper = trialAmplitude;
    targetAxis = trialAxis;
    targetStats = trialStats;
    targetLower = previousAmplitude;
    for iteration = 1:32
        trialAmplitude = 0.5*(targetLower+targetUpper);
        [trialAxis,trialStats] = candidate_axis(trialAmplitude);
        if trialStats.valid && all(trialStats.counts >= targetCounts)
            targetUpper = trialAmplitude;
            targetAxis = trialAxis;
            targetStats = trialStats;
        else
            targetLower = trialAmplitude;
        end
    end
    newAxis = targetAxis;
    info = make_info(targetUpper,targetStats,true);
else
    newAxis = best.axis;
    info = make_info(best.amplitude,best.stats,false);
end

    function best = retain_strongest_valid( ...
            best,validAmplitude,invalidAmplitude,validAxis,validStats)
        for boundaryIteration = 1:32
            amplitude = 0.5*(validAmplitude+invalidAmplitude);
            [axis,stats] = candidate_axis(amplitude);
            if stats.valid
                validAmplitude = amplitude;
                validAxis = axis;
                validStats = stats;
            else
                invalidAmplitude = amplitude;
            end
        end
        score = min(validStats.counts./targetCounts);
        if score > best.score
            best = struct('amplitude',validAmplitude,'axis',validAxis, ...
                'stats',validStats,'score',score);
        end
    end

    function [axis,stats] = candidate_axis(amplitude)
        axis = redistribute(baseAxis,reference,featureBounds,peakCenter, ...
            targetCounts,amplitude,minimumSpacing,remesh,paired);
        faces = diff(axis);
        if numel(faces) > 1
            adjacentRatio = max(max(faces(2:end)./faces(1:end-1), ...
                faces(1:end-1)./faces(2:end)));
        else
            adjacentRatio = 1;
        end
        counts = interval_counts(axis,featureBounds);
        stats = struct('minimumSpacing',min(faces), ...
            'maximumCellRatio',adjacentRatio,'counts',counts, ...
            'valid',min(faces) >= minimumSpacing && ...
            adjacentRatio <= remesh.remeshMaximumCellRatio);
    end

    function output = make_info(amplitude,stats,reached)
        output = struct('amplitude',amplitude, ...
            'minimumSpacing',stats.minimumSpacing, ...
            'maximumCellRatio',stats.maximumCellRatio, ...
            'counts',stats.counts,'targetCounts',targetCounts, ...
            'targetReached',reached);
    end
end

function newAxis = redistribute(baseAxis,reference,bounds,peakCenter, ...
        targetCounts,amplitude,minimumSpacing,remesh,paired)
if isfinite(remesh.remeshOuterAnchorFactor)
    featureExtent = max(abs(bounds),[],'all');
    featureWidth = max(bounds(:,2)-bounds(:,1));
    anchorInner = featureExtent + ...
        remesh.remeshOuterAnchorFactor*featureWidth;
    anchorOuter = min(max(abs(baseAxis)),2.5*anchorInner);
    baseFaces = diff(baseAxis);
    baseRatios = max(baseFaces(2:end)./baseFaces(1:end-1), ...
        baseFaces(1:end-1)./baseFaces(2:end));
    ratioNearlySpent = max(baseRatios) > ...
        0.9*remesh.remeshMaximumCellRatio;
    if paired
        rightAnchor = find(baseAxis <= anchorOuter,1,'last');
        % Include two fixed-side halo nodes.  Without this buffer, repeated
        % current-grid maps can spend the entire cell-ratio budget at the
        % splice while the feature itself is still under-resolved.
        if ratioNearlySpent
            rightAnchor = min(rightAnchor+2,numel(baseAxis));
        end
        leftAnchor = numel(baseAxis)-rightAnchor+1;
        localIndices = leftAnchor:rightAnchor;
    elseif baseAxis(1) >= 0
        rightAnchor = find(baseAxis <= anchorOuter,1,'last');
        if ratioNearlySpent
            rightAnchor = min(rightAnchor+2,numel(baseAxis));
        end
        localIndices = 1:rightAnchor;
    else
        error('ipm:OuterAnchorGeometry', ...
            'A finite outer anchor requires a paired axis or a nonnegative axis.');
    end
    if numel(localIndices) >= 9 && numel(localIndices) < numel(baseAxis)
        localRemesh = remesh;
        localRemesh.remeshOuterAnchorFactor = Inf;
        localBase = baseAxis(localIndices);
        if strcmpi(string(remesh.remeshReferenceMode),'current') && ...
                ratioNearlySpent
            % Recover ratio headroom before composing another local map.
            % The hard admissibility cap is unchanged; sqrt(cap) is only a
            % soft transition target and leaves the feature monitor free to
            % use the remaining ratio budget.
            localBase = relax_outer_cell_ratios(localBase, ...
                sqrt(remesh.remeshMaximumCellRatio),paired,featureExtent);
        end
        localReference = linspace(0,1,numel(localBase));
        localCandidate = redistribute(localBase,localReference,bounds, ...
            peakCenter,targetCounts,amplitude,minimumSpacing, ...
            localRemesh,paired);
        newAxis = baseAxis;
        newAxis(localIndices) = localCandidate;
        return;
    end
end
fineReference = linspace(0,1,16*(numel(reference)-1)+1);
localAxis = [];
for level = 1:size(bounds,1)
    width = bounds(level,2)-bounds(level,1);
    transition = transition_width(width,minimumSpacing,remesh);
    localAxis = [localAxis,linspace(bounds(level,1)-2*transition, ...
        bounds(level,2)+2*transition,129)]; %#ok<AGROW>
    if paired
        localAxis = [localAxis,linspace(-bounds(level,2)-2*transition, ...
            -bounds(level,1)+2*transition,129)]; %#ok<AGROW>
    end
end
if paired && peakCenter > 0
    localAxis = [localAxis,linspace(-peakCenter,peakCenter,257)];
end
localAxis = localAxis(localAxis > baseAxis(1) & ...
    localAxis < baseAxis(end));
localReference = interp1(baseAxis,reference,localAxis,'pchip');
fineReference = unique([fineReference,localReference]);
fineAxis = interp1(reference,baseAxis,fineReference,'pchip');

composite = zeros(size(fineAxis));
numberOfLevels = size(bounds,1);
featureWidths = bounds(:,2)-bounds(:,1);
densityDemand = targetCounts(:)./max(featureWidths,minimumSpacing);
levelWeights = densityDemand/max(densityDemand);
for level = 1:numberOfLevels
    width = bounds(level,2)-bounds(level,1);
    transition = transition_width(width,minimumSpacing,remesh);
    kernel = interval_kernel(fineAxis,bounds(level,:),transition,remesh);
    if paired
        kernel = max(kernel,interval_kernel(fineAxis, ...
            -fliplr(bounds(level,:)),transition,remesh));
    end
    weightedKernel = levelWeights(level)*kernel;
    if strcmpi(string(remesh.remeshAxisProfile),'lattice')
        % The old lat maps use a nearly constant fine metric in the core and
        % a localized smooth transition. Taking the maximum keeps nested
        % plateaus from broadening one another by addition.
        composite = max(composite,weightedKernel);
    else
        composite = composite+weightedKernel;
    end
end
composite = composite/max(composite);
if paired && peakCenter > 0
    absoluteAxis = abs(fineAxis);
    inside = absoluteAxis <= peakCenter;
    bridge = zeros(size(fineAxis));
    normalizedDistance = absoluteAxis(inside)/peakCenter;
    bridge(inside) = remesh.remeshBridgeFloor + ...
        (1-remesh.remeshBridgeFloor)* ...
        normalizedDistance.^remesh.remeshBridgePower;
    composite = max(composite,bridge);
end
monitor = 1+amplitude*composite;
massCoordinate = zeros(size(fineReference));
massCoordinate(2:end) = cumsum(0.5*(monitor(1:end-1)+ ...
    monitor(2:end)).*diff(fineReference));
massCoordinate = massCoordinate/massCoordinate(end);
[massCoordinate,uniqueMassIndices] = unique(massCoordinate,'stable');
fineAxis = fineAxis(uniqueMassIndices);
if numel(massCoordinate) < 2
    error('ipm:AdaptiveMonitorDegenerate', ...
        'Adaptive monitor lost all distinct cumulative coordinates.');
end
newAxis = interp1(massCoordinate,fineAxis,reference,'pchip');
newAxis([1,end]) = baseAxis([1,end]);
if paired
    newAxis = 0.5*(newAxis-fliplr(newAxis));
end
if min(abs(baseAxis)) <= 100*eps(max(abs(baseAxis)))
    [~,originIndex] = min(abs(newAxis));
    newAxis(originIndex) = 0;
end
if any(diff(newAxis) <= 0)
    error('ipm:AdaptiveGridMonotonicity', ...
        'Adaptive redistribution produced a non-monotone axis.');
end
end

function axis = relax_outer_cell_ratios(axis,ratioLimit,paired,featureExtent)
% Smooth only the exterior transition; do not sacrifice feature cells.
widths = diff(axis);
centers = 0.5*(axis(1:end-1)+axis(2:end));
if paired
    positive = find(centers >= 0);
    positiveWidths = relax_one_side(widths(positive), ...
        centers(positive),ratioLimit,featureExtent);
    widths(positive) = positiveWidths;
    widths(numel(widths)-positive+1) = positiveWidths;
else
    widths = relax_one_side(widths,centers,ratioLimit,featureExtent);
end
axis = axis(1)+[0,cumsum(widths)];
axis(end) = axis(1)+sum(widths);
if paired
    axis = 0.5*(axis-fliplr(axis));
end
end

function widths = relax_one_side(widths,centers,ratioLimit,featureExtent)
protected = find(centers <= featureExtent,1,'last');
if isempty(protected)
    protected = 1;
end
if protected >= numel(widths)
    return;
end
outerOriginal = widths(protected+1:end);
outerSpan = sum(outerOriginal);
innerLogWidth = log(widths(protected));
logLimit = log(ratioLimit);
lowerShift = -60;
upperShift = 60;
for iteration = 1:80
    middleShift = 0.5*(lowerShift+upperShift);
    trial = constrained_outer_widths(outerOriginal,middleShift, ...
        innerLogWidth,logLimit);
    if sum(trial) < outerSpan
        lowerShift = middleShift;
    else
        upperShift = middleShift;
    end
end
widths(protected+1:end) = constrained_outer_widths( ...
    outerOriginal,0.5*(lowerShift+upperShift),innerLogWidth,logLimit);
widths(protected+1:end) = outerSpan*widths(protected+1:end) / ...
    sum(widths(protected+1:end));
end

function widths = constrained_outer_widths(original,shift,previous,logLimit)
logWidths = zeros(size(original));
for index = 1:numel(original)
    requested = log(original(index))+shift;
    logWidths(index) = min(max(requested,previous-logLimit), ...
        previous+logLimit);
    previous = logWidths(index);
end
widths = exp(logWidths);
end

function transition = transition_width(width,minimumSpacing,remesh)
if strcmpi(string(remesh.remeshAxisProfile),'lattice')
    transition = max(remesh.remeshLatticeTransitionFraction*width, ...
        8*minimumSpacing);
else
    transition = max(0.5*remesh.remeshWindowFactor*width, ...
        8*minimumSpacing);
end
end

function kernel = interval_kernel(axis,bounds,transition,remesh)
if strcmpi(string(remesh.remeshAxisProfile),'lattice')
    % Stable difference of two logistic steps: a flat fine-grid plateau
    % inside the connected level set and a narrow C-infinity transition.
    kernel = 0.5*(tanh((axis-bounds(1))/transition)- ...
        tanh((axis-bounds(2))/transition));
    kernel = max(kernel,0);
    return;
end
distance = max([bounds(1)-axis;axis-bounds(2);zeros(size(axis))],[],1);
kernel = exp(-(distance/transition).^4);
end

function counts = interval_counts(axis,bounds)
faces = [axis(1:end-1);axis(2:end)];
cellWidths = diff(axis);
counts = zeros(1,size(bounds,1));
for level = 1:size(bounds,1)
    overlap = max(0,min(faces(2,:),bounds(level,2))- ...
        max(faces(1,:),bounds(level,1)));
    counts(level) = sum(overlap./cellWidths);
end
end
