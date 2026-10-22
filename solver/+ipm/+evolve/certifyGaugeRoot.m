function audit = certifyGaugeRoot( ...
        phaseFunction,cR,fRoot,rootPhase,curvature, ...
        roundingGuard,maxSpeedUncertainty)
%IPM.EVOLVE.CERTIFYGAUGEROOT Audit a WENO-limited discrete gauge root.
%   The final matrix-product guard does not cover roundoff in the upstream
%   nonlinear WENO reconstruction. Probe its actual local phase response
%   on both sides of a failed fzero root. Acceptance requires a resolved
%   sign change, endpoint margin above the measured local variation, and
%   a speed uncertainty below the caller's existing 1e-8 budget. For a negative physical slope,
%   callers multiply the function and residuals by -1 before passing them.

audit = struct('accepted',false,'phaseFloor',NaN, ...
    'halfWidth',NaN,'leftPhase',NaN,'rightPhase',NaN, ...
    'slope',NaN,'speedBound',Inf,'minimumSlope',NaN, ...
    'maximumSlope',NaN,'samples',NaN(1,9));
if ~isa(phaseFunction,'function_handle') || ...
        any(~isfinite([cR,fRoot,rootPhase,curvature, ...
        roundingGuard,maxSpeedUncertainty])) || ...
        curvature >= 0 || roundingGuard < 0 || ...
        maxSpeedUncertainty <= 0
    return;
end

expectedSlope = -curvature;
rootResidual = max(abs(fRoot),abs(rootPhase));
halfWidth = min(0.8*maxSpeedUncertainty, ...
    max(0.2*maxSpeedUncertainty, ...
    4*(rootResidual+roundingGuard)/expectedSlope));
offsets = [-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1];
for attempt = 1:5
    samples = zeros(size(offsets));
    for j = 1:numel(offsets)
        if offsets(j) == 0
            samples(j) = rootPhase;
        else
            samples(j) = phaseFunction(cR+halfWidth*offsets(j));
        end
    end
    if any(~isfinite(samples))
        return;
    end
    slope = (samples(end)-samples(1))/(2*halfWidth);
    chord = samples(1)+(offsets+1)* ...
        (samples(end)-samples(1))/2;
    phaseFloor = roundingGuard+max(abs(samples-chord));
    coarseSlopes = diff(samples(1:2:end))/(0.5*halfWidth);
    minimumSlope = min(coarseSlopes);
    maximumSlope = max(coarseSlopes);
    speedBound = halfWidth+phaseFloor/ ...
        max(minimumSlope,realmin);
    audit.phaseFloor = phaseFloor;
    audit.halfWidth = halfWidth;
    audit.leftPhase = samples(1);
    audit.rightPhase = samples(end);
    audit.slope = slope;
    audit.speedBound = speedBound;
    audit.minimumSlope = minimumSlope;
    audit.maximumSlope = maximumSlope;
    audit.samples = samples;
    audit.accepted = samples(1) < -2*phaseFloor && ...
        samples(end) > 2*phaseFloor && ...
        minimumSlope > 0.25*expectedSlope && ...
        maximumSlope < 4*expectedSlope && ...
        slope > 0.5*expectedSlope && ...
        slope < 2*expectedSlope && ...
        rootResidual <= phaseFloor && ...
        speedBound <= maxSpeedUncertainty;
    if audit.accepted || halfWidth >= 0.8*maxSpeedUncertainty
        return;
    end
    halfWidth = min(1.5*halfWidth,0.8*maxSpeedUncertainty);
end
end
