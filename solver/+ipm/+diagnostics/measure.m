function d = measure(rho, flow, ops, mass0)
%IPM.DIAGNOSTICS.MEASURE Accuracy, boundary-contamination, and blow-up indicators.

rhoX = rho*ops.Dx';
rhoY = ops.Dy*rho;
divU = flow.u1*ipm.mesh.oddDx(ops)' + ops.Dy*flow.u2;
speed = hypot(flow.u1,flow.u2);
gradRho = hypot(rhoX,rhoY);
weights = ops.integrationWeights;
mass = sum(rho.*weights,'all');
[d.wallPeak,wallPeakIndex] = max(flow.source(1,:));
d.wallPeakX = ops.x(wallPeakIndex);
if ops.rescaling.enabled && isfinite(flow.trackedPeakX)
    frameX = flow.trackedPeakX;
else
    frameX = d.wallPeakX;
end
frameVelocityU1 = interp1(ops.x,flow.u1(1,:),frameX,'linear');

outer = false(ops.ny,ops.nx);
bandX = max(2,ceil(0.05*ops.nx));
bandY = max(2,ceil(0.05*ops.ny));
outer(:,end-bandX+1:end) = true;
if ~ipm.mesh.isQuadrant(ops)
    outer(:,1:bandX) = true;
end
outer(end-bandY+1:end,:) = true;

d.mass = mass;
d.massDrift = (mass-mass0)/max(abs(mass0),eps);
d.l2 = sqrt(sum(rho.^2.*weights,'all'));
d.rhoMin = min(rho,[],'all');
d.rhoMax = max(rho,[],'all');
d.rhoXInf = flow.rhoXInf;
d.rhoYInf = max(abs(rhoY),[],'all');
d.omegaInf = flow.rhoXInf;
positiveHalfPlane = ops.x > 0;
if any(positiveHalfPlane)
    positiveSource = flow.source(:,positiveHalfPlane);
    d.positiveHalfPlaneOmegaMinimum = min(positiveSource,[],'all');
    positiveSourcePeak = max(positiveSource,[],'all');
    d.positiveHalfPlaneOmegaNegativeRatio = ...
        max(-d.positiveHalfPlaneOmegaMinimum,0)/max(positiveSourcePeak,eps);
else
    d.positiveHalfPlaneOmegaMinimum = NaN;
    d.positiveHalfPlaneOmegaNegativeRatio = NaN;
end
d.gradInf = max(gradRho,[],'all');
d.uInf = max(speed,[],'all');
d.divInf = max(abs(divU),[],'all');
d.wallNormalInf = max(abs(flow.u2(1,:)),[],'all');
[~,originIndex] = min(abs(ops.x));
d.originVelocity = speed(1,originIndex);
if strcmp(ops.symmetryMode,'double_odd_omega') && ~ipm.mesh.isQuadrant(ops)
    d.rhoEvenDefect = max(abs(rho-fliplr(rho)),[],'all');
    d.omegaOddDefect = max(abs(flow.source+fliplr(flow.source)),[],'all');
else
    d.rhoEvenDefect = NaN;
    if ipm.mesh.isQuadrant(ops)
        d.omegaOddDefect = max(abs(flow.source(:,1)));
    else
        d.omegaOddDefect = NaN;
    end
end
d.poissonResidual = flow.poissonResidual;
d.farBoundaryVelocity = max(speed(outer),[],'all');
d.farBoundaryVelocityRatio = d.farBoundaryVelocity/max(d.uInf,eps);
d.frameVelocityU1 = frameVelocityU1;
d.farBoundarySource = max(abs(flow.source(outer)),[],'all');
d.farBoundarySourceRatio = d.farBoundarySource/max(d.omegaInf,eps);
[d.wallDensityMax,wallDensityIndex] = max(rho(1,:));
d.wallDensityMaxX = ops.x(wallDensityIndex);
if d.wallPeak > 0
    peakMetrics = ipm.diagnostics.peakResolution(flow.source(1,:),ops.x,0.5);
    d.wallPeakWidth = peakMetrics.widths;
else
    d.wallPeakWidth = NaN;
end
oscillationMetrics = wall_oscillation_metrics(flow.source(1,:),ops.x);
d.wallPeakLocalMaxima = oscillationMetrics.localMaxima;
d.wallPeakTVRatio = oscillationMetrics.totalVariationRatio;
d.positiveWallNegativeRatio = oscillationMetrics.negativeRatio;
if ops.rescaling.enabled
    d.pinDensity = rho(1,ops.rescaling.pinIndex);
    d.pinOmega = flow.source(1,ops.rescaling.pinIndex);
    d.pinOffset = d.wallPeakX-ops.rescaling.pinX;
else
    d.pinDensity = NaN;
    d.pinOmega = NaN;
    d.pinOffset = NaN;
end

gradientEnergyDensity = gradRho.^2;
d.gradientEnergyX = sum(rhoX.^2.*weights,'all');
d.gradientEnergyY = sum(rhoY.^2.*weights,'all');
d.gradientEnergy = sum(gradientEnergyDensity.*weights,'all');
gradientFourthMoment = sum(gradientEnergyDensity.^2.*weights,'all');
d.gradientEffectiveLength = d.gradientEnergy / ...
    sqrt(max(gradientFourthMoment,eps));
peakScale = d.wallPeakWidth;
if ~isfinite(peakScale) || peakScale <= 0
    peakScale = 8*max(ops.dx,ops.dy);
end
if ops.rescaling.enabled
    energyCenter = ops.rescaling.pinX;
else
    energyCenter = d.wallPeakX;
end
distance = hypot(ops.X-energyCenter,ops.Y);
innerEnergy = sum(gradientEnergyDensity(distance <= 0.5*peakScale) .* ...
    weights(distance <= 0.5*peakScale),'all');
nearMask = distance > 0.5*peakScale & distance <= 2*peakScale;
nearEnergy = sum(gradientEnergyDensity(nearMask).*weights(nearMask),'all');
d.gradientInnerFraction = innerEnergy/max(d.gradientEnergy,eps);
d.gradientNearFraction = nearEnergy/max(d.gradientEnergy,eps);
d.gradientOuterFraction = max(0,1-d.gradientInnerFraction- ...
    d.gradientNearFraction);
end

function metrics = wall_oscillation_metrics(wallOmega,x)
positiveSide = x >= 0;
if ~any(positiveSide)
    metrics = struct('localMaxima',0,'totalVariationRatio',1, ...
        'negativeRatio',0);
    return;
end
positiveOmega = wallOmega(positiveSide);
positiveX = x(positiveSide);
feature = max(positiveOmega,0);
peakMetrics = ipm.diagnostics.peakResolution(feature,positiveX,0.1);
metrics = struct('localMaxima',0,'totalVariationRatio',1, ...
    'negativeRatio',max(-min(positiveOmega),0)/max(peakMetrics.peak,eps));
if peakMetrics.peak <= 0 || any(~isfinite(peakMetrics.bounds))
    return;
end
inside = positiveX >= peakMetrics.bounds(1) & ...
    positiveX <= peakMetrics.bounds(2);
trace = feature(inside);
if numel(trace) < 3
    return;
end
increments = diff(trace);
metrics.localMaxima = nnz(increments(1:end-1) > 0 & ...
    increments(2:end) <= 0);
idealVariation = 2*peakMetrics.peak-trace(1)-trace(end);
metrics.totalVariationRatio = sum(abs(increments))/max(idealVariation,eps);
end
