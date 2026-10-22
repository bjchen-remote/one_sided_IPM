function ops = initializeScaling(rho,ops)
%IPM.EVOLVE.INITIALIZESCALING Fix dynamic-rescaling normalization locations.
%   Dynamic mode requires a resolved positive wall feature and an exact
%   origin. Physical mode records the same feature references only when they
%   exist; an unscaled PDE evolution never depends on a gauge normalization.

enableModulation = strcmp(ops.rescalingMode,'dynamic');
settings = clear_derived_fields(ops.rescaling);
settings.enabled = enableModulation;
ops.rescaling = settings;

[u1,~,~,source] = ipm.field.velocity(rho,ops);
if enableModulation
    if strcmp(ops.symmetryMode,'half_plane')
        eligible = find(ops.x > ops.x(1) & ops.x < ops.x(end));
    else
        eligible = find(ops.x > 0);
    end
    if isempty(eligible)
        error('ipm:RescalingDomain', ...
            'Dynamic rescaling requires interior peak candidates.');
    end
end
initialWall = ipm.diagnostics.measureFeatures( ...
    source(1,:),ops.x,settings.adaptiveLevels,[]);
if initialWall.peak <= 0 || ~isfinite(initialWall.widths(1))
    if enableModulation
        error('ipm:DegeneratePeak', ...
            'Dynamic rescaling requires a resolved positive wall-omega peak.');
    end
    return;
end

if enableModulation
    [~,localPeak] = max(source(1,eligible));
    pinIndex = eligible(localPeak);
else
    pinIndex = initialWall.peakIndex;
end
if strcmp(settings.scalingContract,'wall_integral_centroid_practical_v1')
    initialPeak=ipm.evolve.globalSplineMaximum(source(1,:),ops);
    pinX=initialPeak.x;[~,pinIndex]=min(abs(ops.x-pinX));
else
    pinX = ipm.diagnostics.peakLocation(initialWall.feature,ops.x,pinIndex);
end

if enableModulation
    [~,originIndex] = min(abs(ops.x));
    if abs(ops.x(originIndex)) > 100*eps(max(abs(ops.x)))
        error('ipm:RescalingOrigin', ...
            ['Dynamic rescaling requires a grid point at x=0. ' ...
            'Use odd nx and a grid containing the origin.']);
    end
    u1X = u1*ipm.mesh.oddDx(ops)';
    strainTarget = u1X(1,originIndex);
    if strcmp(settings.cOmegaGauge,'strain_point') && ...
            abs(strainTarget) < 1e-12
        error('ipm:DegenerateNormalization', ...
            'd_x u_1(0,0) is too small for the point-strain c_omega gauge.');
    end
end

horizontal = reference_counts(initialWall,settings);
peakTrackingHalfWidth = settings.peakTrackingWidthFactor * ...
    initialWall.widths(1);
initialVertical = ipm.diagnostics.measureFeatures( ...
    abs(source(:,pinIndex))',ops.y',settings.adaptiveLevels,[]);
if initialVertical.peak <= 0 || ~isfinite(initialVertical.widths(1))
    if enableModulation
        error('ipm:DegenerateVerticalPeak', ...
            'The wall-normal omega profile must be resolved at the initial peak.');
    end
    return;
end
vertical = reference_counts(initialVertical,settings);
ops.rescaling.pinIndex = pinIndex;
ops.rescaling.pinX = pinX;
if enableModulation && strcmp(settings.lengthGauge,'wall_gradient_integral')
    integral = ipm.evolve.wallGradientIntegral(source(1,:),ops,settings);
    ops.rescaling.referenceWallIntegral = integral.value;
    if strcmp(settings.scalingContract,'wall_integral_centroid_practical_v1')
        centroid=ipm.evolve.wallGradientCentroid(source(1,:),ops);
        ops.rescaling.referenceWallGradientCentroid=centroid.value;
    end
end
if enableModulation
    ops.rescaling.originIndex = originIndex;
    ops.rescaling.strainTarget = strainTarget;
end
if enableModulation && strcmp(settings.cOmegaGauge,'anchor_wall_slope')
    anchorWallSlope = interp1(ops.x,source(1,:), ...
        settings.transportAnchorX,'linear');
    wallSlopeScale = max(abs(source(1,:)));
    anchorWallSlopeCondition = ...
        abs(anchorWallSlope)/max(wallSlopeScale,realmin);
    if abs(anchorWallSlope) < settings.cOmegaGaugeFloor || ...
            anchorWallSlopeCondition < ...
            settings.pointStrainConditionFloor
        error('ipm:DegenerateNormalization', ...
            ['The anchor-wall c_omega gauge has a negligible initial ' ...
            'wall slope at transportAnchorX.']);
    end
    ops.rescaling.referenceAnchorWallSlope = anchorWallSlope;
end
if enableModulation && any(strcmp(settings.cOmegaGauge, ...
        {'anchor_wall_window_l2','anchor_wall_window_l4'}))
    if strcmp(settings.cOmegaGauge,'anchor_wall_window_l4')
        momentOrder = 4;
    else
        momentOrder = 2;
    end
    window = ipm.evolve.wallWindowMoment( ...
        source(1,:),ops,settings,momentOrder);
    ops.rescaling.referenceAnchorWallWindowValue = window.value;
end
if enableModulation && ...
        strcmp(settings.cOmegaGauge,'outer_wall_density_window_l2')
    window = ipm.evolve.outerWallDensityWindow( ...
        rho(1,:),ops,settings);
    ops.rescaling.referenceOuterWallDensityWindowValue = window.value;
end
if enableModulation && ...
        strcmp(settings.cOmegaGauge,'anchor_wall_template_projection')
    geometry = ipm.evolve.wallWindowGeometry(ops,settings);
    templateNorm = sqrt(sum(geometry.normalizedWeights.*source(1,:).^2));
    if ~isfinite(templateNorm) || templateNorm < settings.cOmegaGaugeFloor
        error('ipm:DegenerateNormalization', ...
            'The initial anchor-wall template has negligible norm.');
    end
    unitTemplate = source(1,:)/templateNorm;
    referenceProjection = sum(geometry.normalizedWeights .* ...
        source(1,:).*unitTemplate);
    ops.rescaling.referenceAnchorWallTemplate = unitTemplate;
    ops.rescaling.referenceAnchorWallTemplateProjection = ...
        referenceProjection;
end
if enableModulation && ...
        strcmp(settings.cOmegaGauge,'anchor_bulk_gradient_l2')
    rhoY = ops.Dy*rho;
    bulk = ipm.evolve.bulkGradientWindow(source,rhoY,ops,settings);
    ops.rescaling.referenceAnchorBulkGradientL2 = bulk.value;
end
ops.rescaling.referencePeakPoints = horizontal.peak;
ops.rescaling.referenceLevelPoints = horizontal.levels;
ops.rescaling.referenceVerticalPeakPoints = vertical.peak;
ops.rescaling.referenceVerticalLevelPoints = vertical.levels;
ops.rescaling.adaptiveTargetPeakPoints = horizontal.targetPeak;
ops.rescaling.adaptiveTargetLevelPoints = horizontal.targetLevels;
ops.rescaling.adaptiveTargetVerticalPeakPoints = ...
    vertical.targetPeak;
ops.rescaling.adaptiveTargetVerticalLevelPoints = ...
    vertical.targetLevels;
ops.rescaling.safetyPeakPoints = horizontal.safetyPeak;
ops.rescaling.safetyLevelPoints = horizontal.safetyLevels;
ops.rescaling.safetyVerticalPeakPoints = vertical.safetyPeak;
ops.rescaling.safetyVerticalLevelPoints = vertical.safetyLevels;
ops.rescaling.peakTrackingHalfWidth = peakTrackingHalfWidth;
if enableModulation && strcmp(settings.cOmegaGauge,'wall_omega_tracked_spline_peak')
    ops.rescaling.peakEvaluationTime = 0;
    peak = ipm.evolve.trackedSplinePeakFunctional(source(1,:),ops);
    ops.rescaling.trackedPeakAnchor = peak.anchor;
end
end

function counts = reference_counts(metrics,settings)
counts.peak = metrics.areaPoints(1);
counts.levels = metrics.gridPoints;
counts.targetPeak = min(settings.targetPeakPoints, ...
    settings.targetPeakExpansion*counts.peak);
counts.targetLevels = min(settings.targetLevelPoints, ...
    settings.targetPeakExpansion*counts.levels);
counts.safetyLevels = min(settings.minimumLevelPoints,counts.levels);
counts.safetyPeak = min(settings.minimumPeakPoints,counts.peak);
end

function settings = clear_derived_fields(settings)
% Preserve frozen configuration and discard references from an older field.
derived = {'pinIndex','pinX','originIndex','strainTarget','trackedPeakAnchor', ...
    'referenceWallIntegral','referenceWallGradientCentroid', ...
    'referenceWallPeak','referenceAnchorWallSlope', ...
    'referenceAnchorWallWindowL2','referenceAnchorWallWindowValue', ...
    'referenceOuterWallDensityWindowValue', ...
    'referenceAnchorWallTemplate', ...
    'referenceAnchorWallTemplateProjection', ...
    'referenceAnchorBulkGradientL2','referencePeakWidths', ...
    'referenceVerticalPeakWidths','referencePeakPoints', ...
    'referenceLevelPoints','referenceVerticalPeakPoints', ...
    'referenceVerticalLevelPoints','adaptiveTargetPeakPoints', ...
    'adaptiveTargetLevelPoints','adaptiveTargetVerticalPeakPoints', ...
    'adaptiveTargetVerticalLevelPoints','safetyPeakPoints', ...
    'safetyLevelPoints','safetyVerticalPeakPoints', ...
    'safetyVerticalLevelPoints','peakTrackingHalfWidth'};
present = derived(isfield(settings,derived));
if ~isempty(present)
    settings = rmfield(settings,present);
end
end
