function [rhoNew,opsNew,info,rhoRate,validatedFlow] = ...
    adapt(rho,ops,flow,config,amplitudeScale)
%IPM.REMESH.ADAPT Orchestrate propose-transfer-validate remeshing.

rhoRate = [];
validatedFlow = [];
if nargin < 5
    amplitudeScale = 1;
end
validateattributes(amplitudeScale,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename,'amplitudeScale');
if ipm.mesh.isQuadrant(ops)
    [rhoNew,opsNew,info] = ipm.remesh.directLevelSetAdapt( ...
        rho,ops,flow,config,amplitudeScale);
    return
end

if strcmp(ops.symmetryMode,'half_plane') && ...
        strcmp(ops.dynamicScaleGeometry,'isotropic')
    [rhoNew,opsNew,info,rhoRate,validatedFlow] = ...
        area_adapt(rho,ops,flow,config,amplitudeScale);
    return
end

context = [];
accepted = false;
firstProposal = [];
transferFailureIdentifier = '';
firstMetrics = struct('xCore',NaN,'yCore',NaN,'peak',NaN, ...
    'rangeViolation',NaN,'relativeRangeViolation',NaN);
for retry = 0:8
    [proposal,context] = ...
        ipm.remesh.propose(ops,flow,config,retry,context);
    if retry == 0
        firstProposal = proposal;
        if ~proposal.anyImproved
            rhoNew = rho;
            opsNew = ops;
            info = remesh_info(false,ops.remeshCount,context, ...
                proposal,flow,firstMetrics,transferFailureIdentifier);
            return;
        end
    elseif ~proposal.anyImproved
        break;
    end

    try
        [candidateRho,candidateOps,metrics] = ...
            ipm.remesh.transfer(rho,ops,config,proposal);
    catch exception
        if recoverable_high_order_transfer_failure(exception.identifier)
            % Reject a numerically inadmissible candidate and let the
            % existing retry damping reduce its deformation. Programming
            % and configuration errors still escape immediately.
            transferFailureIdentifier = exception.identifier;
            continue;
        end
        rethrow(exception);
    end
    if retry == 0
        firstMetrics = metrics;
    end
    validation = ipm.remesh.validate( ...
        metrics,ops,flow,config,proposal,context,retry,amplitudeScale);
    if validation.accepted
        rhoNew = candidateRho;
        opsNew = candidateOps;
        acceptedProposal = proposal;
        acceptedMetrics = metrics;
        accepted = true;
        break;
    end
end

if ~accepted
    rhoNew = rho;
    opsNew = ops;
    info = remesh_info(false,ops.remeshCount,context, ...
        firstProposal,flow,firstMetrics,transferFailureIdentifier);
    return;
end
opsNew.remeshCount = ops.remeshCount+1;
info = remesh_info(true,opsNew.remeshCount,context, ...
    acceptedProposal,flow,acceptedMetrics,transferFailureIdentifier);
end

function info = remesh_info( ...
        applied,count,context,proposal,flow,metrics,transferFailureIdentifier)
info = struct('applied',applied,'count',count,'peakX',context.peakX, ...
    'window',context.xWindow,'xWindow',context.xWindow, ...
    'yWindow',context.yWindow, ...
    'requestedAmplitude',proposal.xAmplitude, ...
    'amplitude',proposal.xAmplitude, ...
    'xAmplitude',proposal.xAmplitude,'yAmplitude',proposal.yAmplitude, ...
    'minimumSpacing',context.minimumSpacing(1), ...
    'minimumYSpacing',context.minimumSpacing(2), ...
    'oldPeakSpacing',context.oldPeakSpacing, ...
    'newPeakSpacing',proposal.newPeakSpacing, ...
    'oldWallSpacing',context.oldWallSpacing, ...
    'newWallSpacing',proposal.newWallSpacing, ...
    'oldCoreGridPoints',flow.coreGridPoints, ...
    'newCoreGridPoints',metrics.xCore, ...
    'oldVerticalCoreGridPoints',flow.verticalCoreGridPoints, ...
    'newVerticalCoreGridPoints',metrics.yCore, ...
    'rangeViolation',metrics.rangeViolation, ...
    'relativeRangeViolation',metrics.relativeRangeViolation, ...
    'xMaximumCellRatio',proposal.xAxisInfo.maximumCellRatio, ...
    'yMaximumCellRatio',proposal.yAxisInfo.maximumCellRatio, ...
    'xTargetCounts',proposal.xAxisInfo.targetCounts, ...
    'yTargetCounts',proposal.yAxisInfo.targetCounts, ...
    'xPredictedCounts',proposal.xAxisInfo.counts, ...
    'yPredictedCounts',proposal.yAxisInfo.counts, ...
    'xTargetReached',proposal.xAxisInfo.targetReached, ...
    'yTargetReached',proposal.yAxisInfo.targetReached, ...
    'transferFailureIdentifier',transferFailureIdentifier);
end

function answer = recoverable_high_order_transfer_failure(identifier)
recoverable = {'ipm:AreaRemeshGeometry','ipm:AreaRemeshOrigin', ...
    'ipm:AreaRemeshNodeCount','ipm:LevelSetAreaInput', ...
    'ipm:AdaptiveAxisBounds','ipm:AdaptiveMonitorDegenerate', ...
    'ipm:AdaptiveGridMonotonicity','ipm:HighOrderQuadratureWeights', ...
    'ipm:RemeshWallIntegralDrift','ipm:WallIntegralSupport', ...
    'ipm:UnresolvedWallIntegralGauge','ipm:DegenerateWallIntegralGauge', ...
    'ipm:DegenerateWallIntegral','ipm:TrackedPeakBranchLost', ...
    'ipm:HighOrderGridMetric','ipm:FiniteDifferenceScale', ...
    'ipm:HighOrderRemeshBubble','ipm:SixthOrderGridMetric', ...
    'ipm:SixthOrderQuadratureWeights', ...
    'ipm:SixthOrderQuadratureMetric', ...
    'ipm:SixthOrderGridCellRatio', ...
    'ipm:SixthOrderGridStencilCondition', ...
    'ipm:SixthOrderGridQuadratureMargin', ...
    'ipm:SixthOrderGridMetricConsistency', ...
    'ipm:SeparablePoissonResidual','ipm:AnisotropicPoissonPcg'};
answer = any(strcmp(identifier,recoverable));
end


function [rhoNew,opsNew,info,rhoRate,validatedFlow] = ...
        area_adapt(rho,ops,flow,config,amplitudeScale)
% Full two-dimensional superlevel sets; all disconnected components count.
practical = strcmp(config.scaling.scalingContract, ...
    'wall_integral_centroid_practical_v1');
rhoNew = rho;
opsNew = ops;
rhoRate = [];
validatedFlow = [];
features = level_set_areas(flow.source,ops);
info = struct('applied',false,'count',ops.remeshCount, ...
    'areas',features.areas,'referenceAreas',[],'maxAreaChange',NaN, ...
    'xBounds',features.xBounds,'yBounds',features.yBounds, ...
    'referenceXBounds',[],'referenceYBounds',[], ...
    'areaRequested',false,'shapeRequested',false,'requested',false, ...
    'maxShapeChange',NaN,'deferred',false,'pending',false, ...
    'damping',NaN,'targetReached',false, ...
    'countsDiagnosticOnly',practical,'rangeDiagnosticOnly',practical, ...
    'projectionShapeDiagnosticOnly',practical);
if features.peak <= 0 && ...
        ~isfield(ops.rescaling,'referenceLevelSetAreas')
    % No positive feature exists from which to establish a baseline.
    return
end
if ~isfield(ops.rescaling,'referenceLevelSetAreas')
    if ops.remeshCount == 0
        initialRho = ipm.field.initialDensity(ops,config.physics);
        initial = level_set_areas(initialRho*ops.Dx',ops);
        reference = initial.areas;
        referenceXBounds = initial.xBounds;
        referenceYBounds = initial.yBounds;
        clear initialRho initial
    else
        % The last historical remesh field is unavailable in old checkpoints.
        reference = features.areas;
        referenceXBounds = features.xBounds;
        referenceYBounds = features.yBounds;
        if config.output.verbose
            fprintf('  REMESH_AREA_BASELINE legacy count=%d areas=%s\n', ...
                ops.remeshCount,mat2str(reference,6));
        end
    end
    ops.rescaling.referenceLevelSetAreas = reference;
    ops.rescaling.referenceLevelSetXBounds = referenceXBounds;
    ops.rescaling.referenceLevelSetYBounds = referenceYBounds;
elseif (~isfield(ops.rescaling,'referenceLevelSetXBounds') || ...
        ~isfield(ops.rescaling,'referenceLevelSetYBounds')) && valid_geometry(features)
    % Old area-only checkpoints keep their existing area baseline. Their
    % historical shape is unavailable, so begin shape tracking here.
    ops.rescaling.referenceLevelSetXBounds = features.xBounds;
    ops.rescaling.referenceLevelSetYBounds = features.yBounds;
end
reference = ops.rescaling.referenceLevelSetAreas;
[areaRequested,maxChange] = area_change(features.areas,reference);
if isfield(ops.rescaling,'referenceLevelSetXBounds') && ...
        isfield(ops.rescaling,'referenceLevelSetYBounds')
    [shapeRequested,maxShapeChange] = shape_change(features, ...
        ops.rescaling.referenceLevelSetXBounds, ...
        ops.rescaling.referenceLevelSetYBounds);
    info.referenceXBounds = ops.rescaling.referenceLevelSetXBounds;
    info.referenceYBounds = ops.rescaling.referenceLevelSetYBounds;
else
    % An old area-only reference cannot establish bounds from a vanished
    % positive set. Preserve it and reject the requested degenerate shape.
    shapeRequested = true;
    maxShapeChange = Inf;
end
requested = areaRequested || shapeRequested;
if practical
    % The practical trigger is the requested two-dimensional R_X level-set
    % area change. Projection motion remains a map input and a diagnostic.
    requested = areaRequested;
end
observer = remesh_observer(ops);
observer.checks = observer.checks+1;
ops.rescaling.areaRemeshObserver = observer;
opsNew = ops; % Retain the first baseline even when no remesh is requested.
info.referenceAreas = reference;
info.maxAreaChange = maxChange;
info.maxShapeChange = maxShapeChange;
info.areaRequested = areaRequested;
info.shapeRequested = shapeRequested;
info.requested = requested;
info.pending = requested;
if ~requested
    observer.failures = 0;
    observer.retryAt = 0;
    observer.features = [];
    opsNew.rescaling.areaRemeshObserver = observer;
    return
end
if observer.checks < observer.retryAt && ...
        failed_shape_change(features,observer.features) < .02-64*eps(.02)
    % A request stays pending. Retry after bounded check-count backoff, or
    % immediately when this rejected shape changes by another two percent.
    info.deferred = true;
    info.retryAfterChecks = observer.retryAt-observer.checks;
    return
end
if config.output.verbose
    fprintf(['  REMESH_REQUEST count=%d areaChange=%.3f%% shapeChange=%.3f%% ' ...
        'trigger=20%% areas=%s reference=%s\n'], ...
        ops.remeshCount,100*maxChange,100*maxShapeChange, ...
        mat2str(features.areas,6),mat2str(reference,6));
end
if ~valid_geometry(features)
    [opsNew,info] = area_reject(opsNew,info,features,config,'ipm:AreaRemeshGeometry', ...
        'The requested level sets have degenerate projection bounds.');
    return
end
xOld = ops.x(:)';
yOld = ops.y(:)';
[~,peakColumn] = ind2sub(size(flow.source),features.peakIndex);
peakX = xOld(peakColumn);
info.oldPeakSpacing = local_spacing(xOld,peakX);
info.newPeakSpacing = info.oldPeakSpacing;
info.oldWallSpacing = yOld(2)-yOld(1);
info.newWallSpacing = info.oldWallSpacing;
xTargets = ops.rescaling.adaptiveTargetLevelPoints(:)';
yTargets = ops.rescaling.adaptiveTargetVerticalLevelPoints(:)';
if numel(xTargets) ~= numel(features.areas) || ...
        numel(yTargets) ~= numel(features.areas) || ...
        any(~isfinite([xTargets,yTargets])) || any([xTargets,yTargets] <= 0)
    error('ipm:AreaRemeshTargets','Every level needs a positive fixed X/Y target.');
end
xFloors = min(.98*interval_counts(xOld,features.xBounds),xTargets);
yFloors = min(.98*interval_counts(yOld,features.yBounds),yTargets);
info.xValidationFloors = xFloors;
info.yValidationFloors = yFloors;
minimumSpacing = config.grid.targetCenterSpacing(:)';
if isscalar(minimumSpacing)
    minimumSpacing = [minimumSpacing,minimumSpacing];
end
minimumSpacing = config.remesh.remeshMinimumSpacingFactor*minimumSpacing;
axisSettings = config.remesh;
axisSettings.remeshReferenceMode = 'original';
axisSettings.remeshOuterAnchorFactor = Inf;
zeroIndex = find(xOld == 0,1);
if isempty(zeroIndex) || zeroIndex == 1 || zeroIndex == numel(xOld)
    [opsNew,info] = area_reject(opsNew,info,features,config,'ipm:AreaRemeshOrigin', ...
        'Area remeshing requires an exact interior origin.');
    return
end
try
    % Geometry is the primary proposal, even when the current mesh already
    % exceeds its minimum point targets. A pure count monitor can return
    % the original mesh for a 20%-changed, equally resolved morphology.
    xReference = [];
    yReference = [];
    if isfield(ops.rescaling,'referenceLevelSetXBounds')
        xReference = ops.rescaling.referenceLevelSetXBounds;
    end
    if isfield(ops.rescaling,'referenceLevelSetYBounds')
        yReference = ops.rescaling.referenceLevelSetYBounds;
    end
    if practical
        % Rebuild from a positive analytic monitor. Projection endpoints
        % may cross the fixed origin: no old/current landmark pairing is
        % required. Fixed Gaussian masses set a smooth resolution budget,
        % while counts remain diagnostics under the existing policy.
        [targetX,xMorph] = smooth_positive_feature_axis(xOld,features.xBounds,6);
        [targetY,yMorph] = smooth_positive_feature_axis(yOld,features.yBounds,6);
        info.xMorphologyProposal = xMorph;
        info.yMorphologyProposal = yMorph;
        info.xProposalMode = 'analytic_positive_monitor';
        info.yProposalMode = 'analytic_positive_monitor';
        info.xCountCorrection = false;
        info.yCountCorrection = false;
    else
        [targetX,xMapped,xMorph] = levelset_morph_axis( ...
            xOld,xReference,features.xBounds);
        [targetY,yMapped,yMorph] = levelset_morph_axis( ...
            yOld,yReference,features.yBounds);
        info.xMorphologyProposal = xMorph;
        info.yMorphologyProposal = yMorph;
        info.xProposalMode = 'geometry_map';
        info.yProposalMode = 'geometry_map';
        info.xCountCorrection = false;
        info.yCountCorrection = false;
        if ~xMapped
            if practical
                % A missing/conflicting geometry map cannot justify a count-only
                % replacement. Keep the accepted grid and retry new geometry.
                error('ipm:AreaRemeshGeometry', ...
                    'The X landmark map is unavailable: %s.',xMorph.reason);
            end
            [targetX,~] = ipm.remesh.axis(ops.baseX,features.xBounds,peakX, ...
                xTargets,minimumSpacing(1),axisSettings,false);
            info.xProposalMode = 'monitor_fallback';
        elseif ~practical && any(interval_counts(targetX,features.xBounds) < xFloors)
            % Preserve the existing transaction's minimum count floors. If a
            % coarsening warp spends too many cells, correct counts from the
            % *warped* mesh; never erase morphology by restarting at baseX.
            [targetX,~] = ipm.remesh.axis(targetX,features.xBounds,peakX, ...
                xTargets,minimumSpacing(1),axisSettings,false);
            info.xCountCorrection = true;
        end
        if ~yMapped
            if practical
                error('ipm:AreaRemeshGeometry', ...
                    'The Y landmark map is unavailable: %s.',yMorph.reason);
            end
            [targetY,~] = ipm.remesh.axis(ops.baseY',features.yBounds,0, ...
                yTargets,minimumSpacing(2),axisSettings,false);
            info.yProposalMode = 'monitor_fallback';
        elseif ~practical && any(interval_counts(targetY,features.yBounds) < yFloors)
            [targetY,~] = ipm.remesh.axis(targetY,features.yBounds,0, ...
                yTargets,minimumSpacing(2),axisSettings,false);
            info.yCountCorrection = true;
        end
    end
    targetX = lock_origin_index(targetX,zeroIndex,xOld);
    targetY = targetY(:)';
    if numel(targetX) ~= numel(xOld) || numel(targetY) ~= numel(yOld)
        error('ipm:AreaRemeshNodeCount','Area remeshing must preserve node counts.');
    end
catch exception
    if ~recoverable_high_order_transfer_failure(exception.identifier)
        rethrow(exception);
    end
    [opsNew,info] = area_reject(opsNew,info,features,config,exception.identifier,exception.message);
    return
end
reason = 'no_geometry_change';
detail = 'The requested projection grid is unchanged.';
for retry = 0:12
    damping = 2^(-retry);
    if practical && damping < config.remesh.remeshMinimumImprovement
        % A nearly identity transfer is not an effective remesh. Preserve
        % the accepted reference and leave its area request pending under
        % the existing retry backoff, even if point targets are already met.
        reason = 'insufficient_mesh_change';
        detail = sprintf(['retry=%d damping=%.3e is below remeshMinimumImprovement=%.3e; ' ...
            'the area request remains pending.'], ...
            retry,damping,config.remesh.remeshMinimumImprovement);
        break
    end
    x = xOld+damping*(targetX-xOld);
    y = yOld+damping*(targetY-yOld);
    x([1,zeroIndex,end]) = xOld([1,zeroIndex,end]);
    y([1,end]) = yOld([1,end]);
    xValid = axis_quality(x,minimumSpacing(1),config.remesh.remeshMaximumCellRatio);
    yValid = axis_quality(y,minimumSpacing(2),config.remesh.remeshMaximumCellRatio);
    if ~xValid || ~yValid
        reason = 'axis_quality_failed';
        detail = sprintf('retry=%d checks[X/Y]=%d/%d (spacing, monotonicity, adjacent ratio).', ...
            retry,xValid,yValid);
        continue
    end
    if ~axis_changed(x,xOld) && ~axis_changed(y,yOld)
        break
    end
    if ~practical && damping < config.remesh.remeshMinimumImprovement && ...
            ~(all(interval_counts(x,features.xBounds) >= xTargets) && ...
            all(interval_counts(y,features.yBounds) >= yTargets))
        reason = 'insufficient_mesh_change';
        detail = sprintf(['retry=%d damping=%.3e is below remeshMinimumImprovement=%.3e ' ...
            'and the level-set targets remain unmet.'], ...
            retry,damping,config.remesh.remeshMinimumImprovement);
        break
    end
    clear candidateRho candidateOps metrics candidateFeatures candidateRhoRate candidateFlow
    try
        [candidateRho,candidateOps,metrics] = ...
            ipm.remesh.transfer(rho,ops,config,struct('x',x,'y',y'));
        candidateFeatures = level_set_areas(candidateRho*candidateOps.Dx',candidateOps);
        if ~valid_geometry(candidateFeatures)
            error('ipm:AreaRemeshGeometry','Transferred level-set geometry is degenerate.');
        end
    catch exception
        if ~recoverable_high_order_transfer_failure(exception.identifier)
            rethrow(exception);
        end
        reason = exception.identifier;
        detail = exception.message;
        continue
    end
    countTolerance = 256*eps(max([1,xTargets,yTargets]));
    xCounts = interval_counts(x,candidateFeatures.xBounds);
    yCounts = interval_counts(y,candidateFeatures.yBounds);
    xSafe = all(xCounts+countTolerance >= xFloors);
    ySafe = all(yCounts+countTolerance >= yFloors);
    targetReached = all(xCounts+countTolerance >= xTargets) && ...
        all(yCounts+countTolerance >= yTargets);
    peakSafe = ~strcmp(ops.rescalingMode,'physical') || ...
        abs(candidateFeatures.peak-features.peak) <= ...
        config.remesh.remeshPeakChangeTolerance*max(features.peak,eps);
    physicalViolation = metrics.rangeViolation/amplitudeScale;
    rangeSafe = isfinite(physicalViolation) && ...
        physicalViolation <= config.diagnostics.rangeStopTolerance;
    % Counts, peak changes and small range overshoots remain observations
    % in practical mode. Do not replace a smooth landmark warp with a CDF
    % count correction or reject it only for an unmet resolution target.
    info.xCounts = xCounts;info.yCounts = yCounts;
    info.xFloorReached = xSafe;info.yFloorReached = ySafe;
    info.peakCheckPassed = peakSafe;info.rangeCheckPassed = rangeSafe;
    info.physicalRangeViolation = physicalViolation;
    if ~isfinite(physicalViolation) || ...
            (~practical && ~(xSafe && ySafe && peakSafe && rangeSafe))
        reason = 'validation_failed';
        detail = sprintf('retry=%d checks[X/Y/peak/range]=%d/%d/%d/%d.', ...
            retry,xSafe,ySafe,peakSafe,rangeSafe);
        continue
    end
    if ~practical && damping < config.remesh.remeshMinimumImprovement && ~targetReached
        reason = 'insufficient_mesh_change';
        detail = sprintf('retry=%d damping=%.3e does not restore the level-set targets.', ...
            retry,damping);
        break
    end
    try
        % Validate the complete discrete gauge before committing this grid.
        % A failed solve must try the next damped candidate in this loop.
        % This branch is isotropic, so flow needs no scale-state argument.
        [candidateRhoRate,candidateFlow] = ipm.evolve.flow(candidateRho,candidateOps);
    catch exception
        if ~recoverable_high_order_transfer_failure(exception.identifier)
            rethrow(exception);
        end
        reason = exception.identifier;
        detail = sprintf('retry=%d: %s',retry,exception.message);
        continue
    end
    rhoNew = candidateRho;
    opsNew = candidateOps;
    rhoRate = candidateRhoRate;
    validatedFlow = candidateFlow;
    if opsNew.rescaling.enabled && ...
            strcmp(opsNew.rescaling.lengthGauge,'wall_gradient_integral')
        opsNew.rescaling.gaugeRateGuess = [candidateFlow.c_l,candidateFlow.c_r];
    end
    if isfield(candidateFlow,'peakTrackingAnchor')
        opsNew.rescaling.trackedPeakAnchor = candidateFlow.peakTrackingAnchor;
        opsNew.rescaling.peakEvaluationTime = candidateFlow.peakTrackingAnchor.tau;
    end
    opsNew.remeshCount = ops.remeshCount+1;
    opsNew.rescaling.referenceLevelSetAreas = candidateFeatures.areas;
    opsNew.rescaling.referenceLevelSetXBounds = candidateFeatures.xBounds;
    opsNew.rescaling.referenceLevelSetYBounds = candidateFeatures.yBounds;
    observer.failures = 0;
    observer.retryAt = 0;
    observer.features = [];
    opsNew.rescaling.areaRemeshObserver = observer;
    info.applied = true;
    info.count = opsNew.remeshCount;
    info.areas = candidateFeatures.areas;
    info.referenceAreas = candidateFeatures.areas;
    info.xBounds = candidateFeatures.xBounds;
    info.yBounds = candidateFeatures.yBounds;
    info.referenceXBounds = candidateFeatures.xBounds;
    info.referenceYBounds = candidateFeatures.yBounds;
    info.pending = false;
    info.damping = damping;
    info.targetReached = targetReached;
    info.newPeakSpacing = local_spacing(x,peakX);
    info.newWallSpacing = y(2)-y(1);
    return
end
[opsNew,info] = area_reject(opsNew,info,features,config,reason,detail);
end

function [mapped,available,detail] = ...
        levelset_morph_axis(oldAxis,referenceBounds,currentBounds)
%LEVELSET_MORPH_AXIS Follow accepted level-set landmarks on one tensor axis.
% A monotone coordinate map sends every previous projection endpoint to
% its current endpoint. It transports the accepted grid slots, so expansion
% can coarsen, contraction can refine, and equal-area changes can act on
% both tensor axes. It is a proposal, not an acceptance certificate.
% Conflicting repeated landmarks or a crossing of the fixed origin cannot
% define a monotone bijection. Report that conflict instead of silently
% dropping a level or changing a fixed physical coordinate.
oldAxis = oldAxis(:)';
mapped = oldAxis;
available = false;
detail = struct('available',false,'reason','missing_reference', ...
    'sourceAnchors',[],'targetAnchors',[], ...
    'featureExtent',NaN,'outerIdentityRadius',NaN, ...
    'minimumAnchorSlope',NaN,'maximumAnchorSlope',NaN, ...
    'interpolant','monotone_beta7_C6', ...
    'minimumNodeSlope',NaN,'maximumNodeSlope',NaN);
if isempty(referenceBounds)
    return
end
if ~isnumeric(referenceBounds) || ~isreal(referenceBounds) || ...
        ~isnumeric(currentBounds) || ~isreal(currentBounds) || ...
        ~isequal(size(referenceBounds),size(currentBounds)) || ...
        size(referenceBounds,2)~=2 || isempty(referenceBounds) || ...
        any(~isfinite(referenceBounds),'all') || ...
        any(~isfinite(currentBounds),'all') || ...
        any(diff(referenceBounds,1,2)<=0) || any(diff(currentBounds,1,2)<=0)
    detail.reason = 'invalid_reference_geometry';
    return
end
if any(referenceBounds(:)<mapped(1) | referenceBounds(:)>mapped(end)) || ...
        any(currentBounds(:)<mapped(1) | currentBounds(:)>mapped(end))
    detail.reason = 'landmark_outside_axis';
    return
end
extent = max(abs([referenceBounds(:);currentBounds(:)]));
if ~isfinite(extent) || extent<=0
    detail.reason = 'invalid_feature_extent';
    return
end
detail.featureExtent = extent;
detail.outerIdentityRadius = 4*extent;
source = [mapped([1,end]),referenceBounds(:)'];
target = [mapped([1,end]),currentBounds(:)'];
fixed = [true,true,false(1,numel(referenceBounds))];
if mapped(1)<=0 && mapped(end)>=0
    source(end+1)=0;target(end+1)=0;fixed(end+1)=true;
end
% Two identity landmarks on each outer side give the C6 map unit derivative at
% the 4E splice. A lone 4E landmark followed by a 1e10 endpoint can otherwise
% spread an O(1) feature displacement across that enormous outer interval.
guards = extent*[-4,-3,3,4];
guards = guards(guards>mapped(1) & guards<mapped(end));
source = [source,guards];
target = [target,guards];
fixed = [fixed,true(size(guards))];
[source,order] = sort(source);
target = target(order);fixed = fixed(order);
detail.sourceAnchors = source;detail.targetAnchors = target;
mergedSource = zeros(size(source));
mergedTarget = zeros(size(target));
count = 0;first = 1;
while first<=numel(source)
    last = first;
    tolerance = 128*eps(max(1,abs(source(first))));
    while last<numel(source) && source(last+1)-source(first)<=tolerance
        last=last+1;
    end
    indices = first:last;
    targetTolerance = 128*eps(max(1,max(abs(target(indices)))));
    if max(target(indices))-min(target(indices))>targetTolerance
        detail.reason = 'duplicate_landmark_conflict';
        return
    end
    count=count+1;
    hard = indices(find(fixed(indices),1));
    if isempty(hard)
        % Center relative to the first sample to avoid a large-coordinate
        % sum. Compatible duplicate threshold endpoints carry equal weight.
        mergedSource(count)=source(first)+mean(source(indices)-source(first));
        mergedTarget(count)=target(first)+mean(target(indices)-target(first));
    else
        mergedSource(count)=source(hard);
        mergedTarget(count)=target(hard);
    end
    first=last+1;
end
source = mergedSource(1:count);
target = mergedTarget(1:count);
detail.sourceAnchors = source;detail.targetAnchors = target;
for k=1:numel(source)-1
    tolerance = 128*eps(max([1,abs(target(k)),abs(target(k+1))]));
    if target(k+1)-target(k)<=tolerance
        detail.reason = 'landmark_order_conflict';
        return
    end
end
slopes = diff(target)./diff(source);
if any(~isfinite(slopes)) || any(slopes<=0)
    detail.reason = 'invalid_c6_anchor_slopes';
    return
end
detail.minimumAnchorSlope=min(slopes);
detail.maximumAnchorSlope=max(slopes);
[mapped,detail.minimumNodeSlope,detail.maximumNodeSlope] = ...
    monotone_c6_landmarks(source,target,mapped);
% Keep the exterior identically fixed, including floating evaluation error.
outside = abs(oldAxis(:)')>=4*extent;
mapped(outside)=oldAxis(outside);
mapped([1,end])=oldAxis([1,end]);
origin = find(oldAxis==0,1);
if ~isempty(origin),mapped(origin)=0;end
if any(~isfinite(mapped)) || any(diff(mapped)<=0)
    mapped=oldAxis(:)';
    detail.reason = 'nonmonotone_evaluated_map';
    return
end
available=true;
detail.available=true;detail.reason='mapped';
end

function [mapped,minimumNodeSlope,maximumNodeSlope] = ...
        monotone_c6_landmarks(source,target,queries)
%MONOTONE_C6_LANDMARKS Positive-derivative interpolation through landmarks.
% For interval H, positive secant d and shared node slopes m0,m1<=d,
%   F = y0+H*[m0*(t-T)+m1*T+b*S], b=d-(m0+m1)/2>=0,
% where S is the Beta(7,7) CDF and T is its integral from zero to t.
%   F' = (1-S)*m0+S*m1+b*12012*t^6*(1-t)^6 > 0.
% S vanishes to order seven at t=0 and 1-S at t=1. The second through
% sixth derivatives of F vanish at both endpoints; matching the shared
% slopes therefore gives C6 continuity. Seventh derivatives need not match.
% Bernstein sums use nonnegative coefficients. This does not relax any
% discrete metric, spacing, adjacent-cell-ratio or transfer acceptance gate.
source=source(:)';target=target(:)';queries=queries(:)';
widths=diff(source);secants=diff(target)./widths;
nodeSlopes=[secants(1),min(secants(1:end-1),secants(2:end)),secants(end)];
minimumNodeSlope=min(nodeSlopes);maximumNodeSlope=max(nodeSlopes);
mapped=NaN(size(queries));
for interval=1:numel(widths)
    if interval<numel(widths)
        selected=queries>=source(interval) & queries<source(interval+1);
    else
        selected=queries>=source(interval) & queries<=source(interval+1);
    end
    if ~any(selected),continue;end
    localQueries=queries(selected);
    m0=nodeSlopes(interval);m1=nodeSlopes(interval+1);
    d=secants(interval);H=widths(interval);
    if source(interval)==target(interval) && ...
            source(interval+1)==target(interval+1) && m0==1 && m1==1
        mapped(selected)=localQueries;
        continue
    end
    if m0==d && m1==d
        values=target(interval)+d*(localQueries-source(interval));
    else
        t=(localQueries-source(interval))/H;
        % Queries are interval-contained. Clipping removes only arithmetic
        % overshoot in the normalized parameter, never changes a PDE field.
        t=min(max(t,0),1);
        [S,T]=beta7_bernstein(t);
        b=max(0,d-(.5*m0+.5*m1));
        values=target(interval)+H*(m0*(t-T)+m1*T+b*S);
    end
    values(localQueries==source(interval))=target(interval);
    values(localQueries==source(interval+1))=target(interval+1);
    mapped(selected)=values;
end
end

function [S,T] = beta7_bernstein(t)
% S=sum_{k=7}^{13} B_{k,13};
% T=(1/14)*sum_{k=8}^{14}(k-7)*B_{k,14}, with T(1)=1/2.
choose13=[1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,1];
choose14=[1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14,1];
S=zeros(size(t));T=zeros(size(t));oneMinusT=1-t;
for k=7:13
    S=S+choose13(k+1)*(t.^k).*(oneMinusT.^(13-k));
end
for k=8:14
    T=T+((k-7)/14)*choose14(k+1)*(t.^k).*(oneMinusT.^(14-k));
end
% These exact range bounds ensure the nonnegative form of F and F' is
% retained when a sum rounds just outside its analytic range near t=1.
S=min(max(S,0),1);T=min(max(T,0),t);
S(t==0)=0;S(t==1)=1;T(t==0)=0;T(t==1)=.5;
end

function [requested,maxChange] = area_change(current,reference)
if ~isnumeric(current) || ~isnumeric(reference)
    error('ipm:LevelSetAreaReference','Level-set areas must be numeric.');
end
current = current(:)';
reference = reference(:)';
if isempty(current) || numel(current) ~= numel(reference) || ...
        ~isreal(reference) || ...
        any(~isfinite(reference)) || any(reference <= 0)
    error('ipm:LevelSetAreaReference','Reference areas must be positive, finite, and match the levels.');
end
maxChange = max(abs(current./reference-1));
requested = maxChange >= .2-64*eps(.2);
end

function [requested,maxChange] = shape_change(features,xReference,yReference)
references = {xReference,yReference};
current = {features.xBounds,features.yBounds};
changes = zeros(1,4);
for axis = 1:2
    reference = references{axis};
    if ~isnumeric(reference) || ~isreal(reference) || ...
            ~isequal(size(reference),size(current{axis})) || ...
            any(~isfinite(reference),'all') || any(diff(reference,1,2) <= 0)
        error('ipm:LevelSetShapeReference', ...
            'Reference level-set bounds must be finite, positive-width, and match the levels.');
    end
    width = diff(reference,1,2);
    now = current{axis};
    if any(~isfinite(now),'all') || any(diff(now,1,2) <= 0)
        maxChange = Inf;
        requested = true;
        return
    end
    changes(2*axis-1) = max(abs(diff(now,1,2)./width-1));
    changes(2*axis) = max(abs(mean(now,2)-mean(reference,2))./width);
end
maxChange = max(changes);
requested = maxChange >= .2-64*eps(.2);
end

function observer = remesh_observer(ops)
% This is a transient retry scheduler, not a frozen physical reference.
% Restoring a checkpoint deliberately starts with an immediate retry.
observer = struct('checks',0,'failures',0,'retryAt',0,'features',[]);
if ~isfield(ops.rescaling,'areaRemeshObserver')
    return
end
stored = ops.rescaling.areaRemeshObserver;
if ~isstruct(stored) || ~isscalar(stored) || ...
        ~all(isfield(stored,{'checks','failures','retryAt','features'}))
    error('ipm:AreaRemeshObserver','The remesh retry observer is malformed.');
end
values = [stored.checks,stored.failures,stored.retryAt];
if ~isnumeric(values) || ~isreal(values) || numel(values) ~= 3 || ...
        any(~isfinite(values)) || any(values < 0) || any(values ~= fix(values))
    error('ipm:AreaRemeshObserver','Remesh retry counters must be finite nonnegative integers.');
end
observer = stored;
end

function change = failed_shape_change(features,previous)
if isempty(previous)
    change = Inf;
elseif ~valid_geometry(features) || ~valid_geometry(previous)
    % Repeated disappearance of a positive feature should not print once
    % per step. Any newly distinct degenerate geometry retries immediately.
    if isequaln(features.areas,previous.areas) && ...
            isequaln(features.xBounds,previous.xBounds) && ...
            isequaln(features.yBounds,previous.yBounds)
        change = 0;
    else
        change = Inf;
    end
else
    [~,areaChange] = area_change(features.areas,previous.areas);
    [~,shapeChange] = shape_change(features,previous.xBounds,previous.yBounds);
    change = max(areaChange,shapeChange);
end
end

function valid = valid_geometry(features)
valid = all(isfinite([features.areas(:);features.xBounds(:);features.yBounds(:)])) && ...
    all(features.areas > 0) && all(diff(features.xBounds,1,2) > 0) && ...
    all(diff(features.yBounds,1,2) > 0);
end

function target = lock_origin_index(target,oldZeroIndex,oldAxis)
target = target(:)';
newZeroIndex = find(target == 0,1);
if isempty(newZeroIndex) || newZeroIndex == 1 || newZeroIndex == numel(target)
    error('ipm:AreaRemeshOrigin','The requested grid lost its exact interior origin.');
end
if newZeroIndex ~= oldZeroIndex
    left = interp1(linspace(0,1,newZeroIndex),target(1:newZeroIndex), ...
        linspace(0,1,oldZeroIndex),'pchip');
    right = interp1(linspace(0,1,numel(target)-newZeroIndex+1),target(newZeroIndex:end), ...
        linspace(0,1,numel(oldAxis)-oldZeroIndex+1),'pchip');
    target = [left,right(2:end)];
end
target([1,oldZeroIndex,end]) = oldAxis([1,oldZeroIndex,end]);
end

function passed = axis_quality(axis,minimumSpacing,maximumRatio)
faces = diff(axis);
passed = all(isfinite(axis)) && all(faces > 0) && min(faces) >= minimumSpacing;
if passed && numel(faces) > 1
    passed = max([faces(2:end)./faces(1:end-1),faces(1:end-1)./faces(2:end)]) <= maximumRatio;
end
end

function changed = axis_changed(candidate,old)
changed = any(abs(candidate-old) > 128*eps(max(1,abs(old))));
end

function counts = interval_counts(axis,bounds)
faces = [axis(1:end-1);axis(2:end)];
counts = zeros(1,size(bounds,1));
for level = 1:size(bounds,1)
    overlap = max(0,min(faces(2,:),bounds(level,2))-max(faces(1,:),bounds(level,1)));
    counts(level) = sum(overlap./diff(axis));
end
end

function spacing = local_spacing(axis,point)
[~,index] = min(abs(axis-point));
left = max(index-1,1);
right = min(index+1,numel(axis));
spacing = (axis(right)-axis(left))/(right-left);
end

function [ops,info] = area_reject(ops,info,features,config,identifier,message)
info.rejectionReason = identifier;
info.rejectionDetail = strtrim(regexprep(message,'[\r\n]+',' '));
observer = ops.rescaling.areaRemeshObserver;
observer.failures = min(observer.failures+1,5);
delay = 8*2^(observer.failures-1);
observer.retryAt = observer.checks+delay;
observer.features = struct('areas',features.areas, ...
    'xBounds',features.xBounds,'yBounds',features.yBounds);
ops.rescaling.areaRemeshObserver = observer;
info.retryAfterChecks = delay;
if config.output.verbose
    fprintf(['  REMESH_REJECT count=%d reason=%s: %s ' ...
        'pending=1 retryAfterChecks=%d earlyRetryChange=2%%\n'], ...
        info.count,identifier,info.rejectionDetail,delay);
end
end

function measurement = level_set_areas(source,ops)
%LEVEL_SET_AREAS Geometric areas of two-dimensional R_X sets.
%   A_l = area{R_X >= l*max(R_X)} over the computational domain. All
%   components contribute. Each rectangle is split along its lower-left to
%   upper-right diagonal; linear interpolation on each triangle gives exact
%   clipped areas for that discrete field. Counts of grid cells are not used.

x = ops.x(:)';
y = ops.y(:);
levels = ops.rescaling.adaptiveLevels(:)';
if ~isnumeric(source) || ~isreal(source) || ...
        ~isequal(size(source),[numel(y),numel(x)]) || ...
        any(~isfinite(source),'all') || numel(x)<2 || numel(y)<2 || ...
        any(~isfinite(x)) || any(~isfinite(y)) || ...
        any(diff(x)<=0) || any(diff(y)<=0) || isempty(levels) || ...
        any(~isfinite(levels)) || any(levels<=0 | levels>=1) || ...
        any(diff(levels)<=0)
    error('ipm:LevelSetAreaInput', ...
        'Level-set areas need a finite R_X field on increasing X/Y axes and levels in (0,1).');
end
[peak,peakIndex] = max(source(:));
measurement = struct('levels',levels,'areas',zeros(size(levels)), ...
    'xBounds',NaN(numel(levels),2),'yBounds',NaN(numel(levels),2), ...
    'peak',peak,'peakIndex',peakIndex);
if peak <= 0
    return
end

% A piecewise linear triangle cannot cross a positive threshold unless a
% vertex reaches it. Include all lowest-level active vertices and one node
% of halo on each side; every contributing cell and crossing stays inside.
% Keep peakIndex in the original field for the subsequent peak-column lookup.
active = source >= levels(1)*peak;
activeRows = find(any(active,2));
activeColumns = find(any(active,1));
rows = max(1,activeRows(1)-1):min(numel(y),activeRows(end)+1);
columns = max(1,activeColumns(1)-1):min(numel(x),activeColumns(end)+1);
source = source(rows,columns);
x = x(columns);
y = y(rows);

q00 = source(1:end-1,1:end-1);
q10 = source(1:end-1,2:end);
q01 = source(2:end,1:end-1);
q11 = source(2:end,2:end);
triangleArea = 0.5*(diff(y)*diff(x));
first = sort([q00(:),q10(:),q11(:)],2);
second = sort([q00(:),q11(:),q01(:)],2);
[X,Y] = meshgrid(x,y);
for index = 1:numel(levels)
    threshold = levels(index)*peak;
    fractions = triangle_fraction(first,threshold)+ ...
        triangle_fraction(second,threshold);
    measurement.areas(index) = sum(triangleArea(:).*fractions);

    % Vertices and edge crossings are all extrema of the clipped polygons.
    inside = source >= threshold;
    bounds = [min(X(inside)),max(X(inside)), ...
        min(Y(inside)),max(Y(inside))];
    bounds = crossing_bounds(bounds,threshold, ...
        source(:,1:end-1),source(:,2:end), ...
        X(:,1:end-1),X(:,2:end),Y(:,1:end-1),Y(:,2:end));
    bounds = crossing_bounds(bounds,threshold, ...
        source(1:end-1,:),source(2:end,:), ...
        X(1:end-1,:),X(2:end,:),Y(1:end-1,:),Y(2:end,:));
    bounds = crossing_bounds(bounds,threshold,q00,q11, ...
        X(1:end-1,1:end-1),X(2:end,2:end), ...
        Y(1:end-1,1:end-1),Y(2:end,2:end));
    measurement.xBounds(index,:) = bounds(1:2);
    measurement.yBounds(index,:) = bounds(3:4);
end
end

function fractions = triangle_fraction(sorted,threshold)
a = sorted(:,1);
b = sorted(:,2);
c = sorted(:,3);
fractions = zeros(size(a));
full = threshold <= a;
fractions(full) = 1;
% Full takes precedence for an exactly threshold-valued plateau.
partial = ~full & threshold < c;
lower = partial & threshold <= b;
upper = partial & threshold > b;
fractions(lower) = 1- ...
    ((threshold-a(lower))./(b(lower)-a(lower))).* ...
    ((threshold-a(lower))./(c(lower)-a(lower)));
fractions(upper) = ...
    ((c(upper)-threshold)./(c(upper)-a(upper))).* ...
    ((c(upper)-threshold)./(c(upper)-b(upper)));
fractions = min(max(fractions,0),1);
end

function bounds = crossing_bounds(bounds,threshold,qA,qB,xA,xB,yA,yB)
crossing = (qA>=threshold) ~= (qB>=threshold);
if ~any(crossing,'all')
    return
end
fraction = (threshold-qA(crossing))./(qB(crossing)-qA(crossing));
crossX = xA(crossing)+fraction.*(xB(crossing)-xA(crossing));
crossY = yA(crossing)+fraction.*(yB(crossing)-yA(crossing));
bounds = [min(bounds(1),min(crossX)),max(bounds(2),max(crossX)), ...
    min(bounds(3),min(crossY)),max(bounds(4),max(crossY))];
end

function [axis,info] = smooth_positive_feature_axis(oldAxis,bounds,featureMass)
%SMOOTH_POSITIVE_FEATURE_AXIS Positive analytic-CDF practical proposal.
% The surrounding transaction keeps its existing acceptance policy.
% Monitor = asinh-tail density + fixed origin Gaussian + fixed-mass
% feature Gaussians. Its integral is analytic (asinh + erf). Features can
% cross X=0 because no old/current projection endpoints are paired.
% The existing origin INDEX is retained using one smooth far-field
% balancing Gaussian; both sides use the SAME global CDF mass step.

if nargin<3 || isempty(featureMass),featureMass=6;end
oldAxis=oldAxis(:)';
assert(isnumeric(oldAxis) && isreal(oldAxis) && all(isfinite(oldAxis)) && ...
    numel(oldAxis)>=9 && all(diff(oldAxis)>0), ...
    'ipm:AreaRemeshGeometry','The analytic monitor needs an increasing finite axis.');
assert(isnumeric(bounds) && isreal(bounds) && size(bounds,2)==2 && ...
    size(bounds,1)>=1 && all(isfinite(bounds),'all') && ...
    all(diff(bounds,1,2)>0) && all(bounds(:)>=oldAxis(1)) && ...
    all(bounds(:)<=oldAxis(end)), ...
    'ipm:AreaRemeshGeometry','The analytic monitor needs nondegenerate in-domain projections.');
assert(isscalar(featureMass) && isfinite(featureMass) && featureMass>0, ...
    'ipm:AreaRemeshGeometry','The fixed analytic monitor mass must be positive.');
left=oldAxis(1);right=oldAxis(end);number=numel(oldAxis);
zero=find(oldAxis==0,1);
horizontal=left<0 && right>0;
assert((horizontal && ~isempty(zero)) || left==0, ...
    'ipm:AreaRemeshOrigin','The analytic monitor must retain its exact origin.');
if horizontal
    coreScale=.4;
    components=[.5*.55*sqrt(2*pi),0,.55];
else
    coreScale=.375;
    components=zeros(0,3);
end
if size(bounds,1)==3
    weights=featureMass*[.2,.3,.5];
else
    weights=featureMass*ones(1,size(bounds,1))/size(bounds,1);
end
components=[components;weights(:),mean(bounds,2),.5*diff(bounds,1,2)];
balanceMass=0;balanceMean=NaN;balanceSigma=NaN;
if horizontal
    q0=(zero-1)/(number-1);
    total=mass(right);atZero=mass(0);
    deficit=q0*total-atZero;
    if abs(deficit)>8*eps*total
        extent=max([coreScale;abs(components(:,2))+3*components(:,3)]);
        if deficit>0
            available=-left;signSide=-1;
        else
            available=right;signSide=1;
        end
        balanceSigma=min(extent,available/12);
        balanceMean=signSide*min(4*extent,available/2);
        h0=normal_cdf(-balanceMean/balanceSigma)- ...
            normal_cdf((left-balanceMean)/balanceSigma);
        ht=normal_cdf((right-balanceMean)/balanceSigma)- ...
            normal_cdf((left-balanceMean)/balanceSigma);
        balanceMass=deficit/(h0-q0*ht);
        assert(isfinite(balanceMass) && balanceMass>=0, ...
            'ipm:AreaRemeshGeometry','The analytic origin-balance mass must be nonnegative.');
        components(end+1,:)=[balanceMass,balanceMean,balanceSigma];
    end
end
total=mass(right);
targets=(0:number-1)*(total/(number-1));
lower=left*ones(1,number);upper=right*ones(1,number);
for iteration=1:128
    middle=lower+(upper-lower)/2;
    if all(middle==lower | middle==upper),break;end
    below=mass(middle)<targets;
    lower(below)=middle(below);upper(~below)=middle(~below);
end
axis=lower+(upper-lower)/2;
axis([1,end])=[left,right];
if horizontal,axis(zero)=0;end
assert(all(isfinite(axis)) && all(diff(axis)>0), ...
    'ipm:AreaRemeshGeometry','The analytic CDF quantiles must be finite and increasing.');
widths=diff(axis);
ratio=max([widths(2:end)./widths(1:end-1), ...
    widths(1:end-1)./widths(2:end)]);
info=struct('kind','analytic_positive_Cinf_feature_monitor', ...
    'featureMass',featureMass,'components',components, ...
    'coreScale',coreScale,'originIndex',zero, ...
    'balanceMass',balanceMass,'balanceMean',balanceMean, ...
    'balanceSigma',balanceSigma,'totalMonitorMass',total, ...
    'maximumAdjacentRatio',ratio,'minimumSpacing',min(widths), ...
    'maximumCdfQuantileError',max(abs(mass(axis)-targets))/total, ...
    'pointCountsAreNotOptimizationTargets',true);
if horizontal,info.originMassFraction=mass(0)/total;end

    function value=mass(location)
        value=asinh(location/coreScale)-asinh(left/coreScale);
        for k=1:size(components,1)
            weight=components(k,1);mu=components(k,2);sigma=components(k,3);
            value=value+weight*(normal_cdf((location-mu)/sigma)- ...
                normal_cdf((left-mu)/sigma));
        end
    end
end

function value=normal_cdf(location)
value=.5*(1+erf(location/sqrt(2)));
end
