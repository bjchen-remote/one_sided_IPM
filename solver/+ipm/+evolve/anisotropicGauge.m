function [rates,baseRhs,details,transportU1,transportU2,extras] = ...
    anisotropicGauge(rho,physicalRhs,flow,ops)
%IPM.EVOLVE.ANISOTROPICGAUGE Select maintained fixed/peak-translation rates.

r = ops.rescaling;
trackingOps = ops;
trackingOps.rescaling.enabled = false;
details = ipm.diagnostics.trackFeatures(rho,physicalRhs,flow,trackingOps);

cX = r.fixedCX;
cY = r.fixedCY;
cOmega = r.fixedCOmega;
cR = r.fixedCR;
details.c_lNominal = cX;
details.lengthScaleGain = 1;
details.widthCorrection = 0;
details.widthRateScale = abs(cX);
details.widthControlMode = 0;

switch r.anisotropicGaugeMode
    case 'fixed'
        targetVelocity = NaN;
    case 'peak_translation'
        if ~isfinite(details.trackedPeakX)
            error('ipm:DegeneratePeakTranslationGauge', ...
                'The half-plane peak-translation gauge lost its tracked peak.');
        end
        peakU1 = interp1(ops.x,flow.u1(1,:), ...
            details.trackedPeakX,'linear');
        cR = -peakU1-cX*details.trackedPeakX;
        targetVelocity = 0;
    otherwise
        error('ipm:BadAnisotropicGaugeMode', ...
            ['Anisotropic gauge mode must be ''fixed'' or ' ...
            '''peak_translation''.']);
end

[~,baseRhs,transportU1,transportU2,conservativeSource] = ...
    ipm.evolve.assembleRhs( ...
        rho,flow.u1,flow.u2,cX,cY,cOmega,cR,ops);
u1X = flow.u1*ipm.mesh.oddDx(ops)';
details.strain = u1X(1,r.originIndex);
details.strainError = details.strain-r.strainTarget;
details.strainCondition = abs(details.strain)/max(abs(u1X),[],'all');
if isfinite(details.trackedPeakX)
    peakU1 = interp1(ops.x,flow.u1(1,:), ...
        details.trackedPeakX,'linear');
    rawVelocity = peakU1+cX*details.trackedPeakX+cR;
    details.peakTransportVelocity = rawVelocity;
    details.trackedPeakOffset = details.trackedPeakX-r.pinX;
    residual = rawVelocity;
    details.travelingWaveResidual = residual;
else
    rawVelocity = NaN;
    residual = NaN;
end
details.pinTransportVelocity = interp1( ...
    ops.x,transportU1(1,:),r.pinX,'linear');
details.pinGaugeResidual = residual;
details.lengthGaugeResidual = residual;
details.lengthGaugeStrain = NaN;

rates = struct('cx',cX,'cy',cY,'comega',cOmega,'cr',cR, ...
    'conservativeSource',conservativeSource);
extras = struct('mode',r.anisotropicGaugeMode, ...
    'peakTranslationResidual',residual, ...
    'peakTranslationRawVelocity',rawVelocity, ...
    'peakTranslationTargetVelocity',targetVelocity);
end
