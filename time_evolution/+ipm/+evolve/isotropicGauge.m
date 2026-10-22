function [c_l,c_omega,c_r,baseRhs,details,transportU1,transportU2] = ...
    isotropicGauge(rho,physicalRhs,flow,ops,details)
%IPM.EVOLVE.ISOTROPICGAUGE Compute exact isotropic gauge rates and residuals.

if ~ops.rescaling.enabled
    c_l = 0;
    c_omega = 0;
    c_r = 0;
    baseRhs = physicalRhs;
    transportU1 = flow.u1;
    transportU2 = flow.u2;
    return;
end

r = ops.rescaling;
if strcmp(r.scalingContract,'wall_integral_centroid_practical_v1')
    [c_l,c_omega,c_r,baseRhs,details,transportU1,transportU2] = ...
        practical_integral_centroid(rho,flow,ops,details);
    return;
end
trackedPeakX = details.trackedPeakX;
trackedPeakU1 = interp1(ops.x,flow.u1(1,:),trackedPeakX,'linear');
needsStrain = any(strcmp(r.lengthGauge, ...
    {'local_strain','wall_gradient_curvature','wall_gradient_integral'})) || ...
    any(strcmp(r.cOmegaGauge,{'strain_point','anchor_wall_slope', ...
        'anchor_wall_window_l2','anchor_wall_window_l4', ...
        'anchor_wall_template_projection','anchor_bulk_gradient_l2'}));
if needsStrain
    u1X = flow.u1*ipm.mesh.oddDx(ops)';
else
    u1X = [];
end
peakLocationRhoXX = NaN;
peakLocationRhoXXScale = NaN;
peakLocationPhaseDefect = NaN;
peakLocationResponseCondition = NaN;
peakLocationFullPhaseResidual = NaN;
peakLocationNormalizedFullPhaseResidual = NaN;
phaseNodes = [];
phaseWeights = [];
phaseContext = [];
switch r.lengthGauge
    case 'omega_peak_location'
        rhoXX = (rho*ops.Dx')*ipm.mesh.oddDx(ops)';
        rhoXXX = rhoXX*ops.Dx';
        physicalRhsXX = (physicalRhs*ops.Dx')*ipm.mesh.oddDx(ops)';
        peakRhoXXX = interp1( ...
            ops.x,rhoXXX(1,:),trackedPeakX,'linear');
        peakRhoXX = interp1( ...
            ops.x,rhoXX(1,:),trackedPeakX,'linear');
        peakPhysicalRhsXX = interp1( ...
            ops.x,physicalRhsXX(1,:),trackedPeakX,'linear');
        peakLocationResponse = trackedPeakX*peakRhoXXX;
        responseScale = max(abs(ops.x.*rhoXXX(1,:)));
        if abs(peakLocationResponse) <= 1e-8*max(responseScale,eps)
            error('ipm:DegenerateLengthGauge', ...
                ['The wall-omega peak-location gauge has negligible ' ...
                'curvature response.']);
        end
        c_lNominal = peakPhysicalRhsXX/peakLocationResponse;
        peakLocationRhoXX = peakRhoXX;
        peakLocationRhoXXScale = max(abs(rhoXX(1,:)));
        peakLocationPhaseDefect = abs(peakRhoXX)/ ...
            max(peakLocationRhoXXScale,realmin);
        peakLocationResponseCondition = abs(peakLocationResponse)/ ...
            max(responseScale,realmin);
        lengthGaugeStrain = peakLocationResponse;
        lengthGaugeResidual = peakPhysicalRhsXX- ...
            c_lNominal*peakLocationResponse;
    case 'wall_omega_width'
        currentWidth = details.wallWidthCurrent;
        lengthGaugeStrain = currentWidth;
        if lengthGaugeStrain <= 100*eps(max(abs(ops.x)))
            error('ipm:DegenerateLengthGauge', ...
                ['The wall-omega half-height width is too small for ' ...
                'normalization.']);
        end
        c_lNominal = -details.wallWidthRate/lengthGaugeStrain;
        lengthGaugeResidual = details.wallWidthRate + ...
            c_lNominal*lengthGaugeStrain;
    case {'wall_gradient_curvature','wall_gradient_integral'}
        lengthGaugeStrain = NaN;
        c_lNominal = -interp1( ...
            ops.x,u1X(1,:),trackedPeakX,'linear');
        lengthGaugeResidual = NaN;
    case 'symmetry_peak'
        originX = ops.x(r.originIndex);
        originU1 = flow.u1(1,r.originIndex);
        peakSeparation = trackedPeakX-originX;
        if peakSeparation <= 100*eps(max(abs(ops.x)))
            error('ipm:DegenerateLengthGauge', ...
                ['The tracked positive omega peak is too close to the ' ...
                'symmetry origin.']);
        end
        lengthGaugeStrain = (trackedPeakU1-originU1)/peakSeparation;
        c_lNominal = -lengthGaugeStrain;
        lengthGaugeResidual = trackedPeakU1-originU1 + ...
            c_lNominal*peakSeparation;
    case 'transport_anchor'
        anchorX = r.transportAnchorX;
        anchorU1 = interp1(ops.x,flow.u1(1,:),anchorX,'linear');
        c_lNominal = -anchorU1/anchorX;
        lengthGaugeStrain = anchorX;
        lengthGaugeResidual = anchorU1+c_lNominal*anchorX;
    case 'wall_density_width'
        lowerX = details.wallDensityLowerX;
        upperX = details.wallDensityUpperX;
        lengthGaugeWidth = details.wallDensityWidth;
        lowerU1 = interp1(ops.x,flow.u1(1,:),lowerX,'linear');
        upperU1 = interp1(ops.x,flow.u1(1,:),upperX,'linear');
        lengthGaugeStrain = (upperU1-lowerU1)/lengthGaugeWidth;
        c_lNominal = -lengthGaugeStrain;
        lengthGaugeResidual = upperU1-lowerU1 + ...
            c_lNominal*lengthGaugeWidth;
    case 'local_strain'
        lengthGaugeStrain = interp1( ...
            ops.x,u1X(1,:),trackedPeakX,'linear');
        c_lNominal = -lengthGaugeStrain;
        lengthGaugeResidual = lengthGaugeStrain+c_lNominal;
    otherwise
        error('ipm:BadLengthGaugeRuntime', ...
            'Unknown runtime length gauge: %s.',r.lengthGauge);
end

c_l = c_lNominal;
coupledWidthPeakGauge = strcmp(ops.symmetryMode,'half_plane') && ...
    strcmp(r.lengthGauge,'wall_omega_width') && ...
    isfield(r,'travelingWaveGauge') && ...
    strcmp(r.travelingWaveGauge,'gradient_peak');
coupledCurvaturePeakGauge = strcmp(r.lengthGauge, ...
    'wall_gradient_curvature');
coupledIntegralPeakGauge = strcmp(r.lengthGauge,'wall_gradient_integral');
coupledRestoringPeakGauge = coupledCurvaturePeakGauge || ...
    coupledIntegralPeakGauge;
widthCorrection = 0;
widthControlMode = 0;
resolutionFactor = details.resolutionFactor;
safetyFactor = details.safetyFactor;
peakExpansionRatio = details.peakExpansionRatio;
if strcmp(ops.symmetryMode,'double_odd_omega')
    c_r = 0;
elseif isfield(r,'travelingWaveGauge') && ...
        strcmp(r.travelingWaveGauge,'gradient_peak')
    % Pin the actual three-node quadratic vertex used by peakLocation.
    % Its derivative functional annihilates R_X, so the amplitude rate
    % cancels exactly from the discrete phase equation.  Evaluate the
    % remaining phase with the same nonlinear transport used by the RHS.
    if coupledCurvaturePeakGauge
        [c_l,c_r,phaseNodes,phaseWeights,secondWeights, ...
            phaseContext,peakP,peakK] = ...
            ipm.evolve.peakCurvatureGauge( ...
            rho,flow,ops,details,c_lNominal);
        phaseContext.curvature = ...
            phaseContext.smoothCurvature;
        phaseContext.secondWeights = secondWeights;
        phaseContext.peakP = peakP;
        phaseContext.phaseTarget = -r.peakPositionRestoringRate* ...
            phaseContext.peakX*phaseContext.peakK;
        lengthGaugeStrain = phaseContext.curvature;
    elseif coupledIntegralPeakGauge
        [c_l,c_r,phaseNodes,phaseWeights,~, ...
            phaseContext,peakP,~] = ipm.evolve.wallIntegralGauge( ...
            rho,flow,ops,details,c_lNominal);
        phaseContext.peakP = peakP;
        phaseContext.phaseTarget = -r.peakPositionRestoringRate* ...
            phaseContext.peakX*phaseContext.peakK;
        lengthGaugeStrain = phaseContext.integralValue;
    elseif coupledWidthPeakGauge
        [c_l,c_r,phaseNodes,phaseWeights,phaseContext] = ...
            pin_width_and_gradient_peak( ...
            rho,flow,ops,details,c_l,trackedPeakU1);
    else
        [c_r,phaseNodes,phaseWeights,phaseContext] = pin_gradient_peak( ...
            rho,flow,ops,details,c_l,trackedPeakU1);
    end
    % Retain the continuum phase defect as a separate resolution audit.
    if ~coupledRestoringPeakGauge
        rhoXX = (rho*ops.Dx')*ipm.mesh.oddDx(ops)';
        rhoXXX = rhoXX*ops.Dx';
        physicalRhsXX = (physicalRhs*ops.Dx')*ipm.mesh.oddDx(ops)';
        peakRhoXXX = interp1(ops.x,rhoXXX(1,:),trackedPeakX,'linear');
        peakRhoXX = interp1(ops.x,rhoXX(1,:),trackedPeakX,'linear');
        peakPhysicalRhsXX = interp1( ...
            ops.x,physicalRhsXX(1,:),trackedPeakX,'linear');
        curvatureScale = max(abs(rhoXXX(1,:)));
        if ~isfinite(peakRhoXXX) || ...
                abs(peakRhoXXX) <= 1e-8*max(curvatureScale,realmin)
            error('ipm:DegenerateTravelingWaveGauge', ...
                'The tracked wall-gradient peak has negligible curvature.');
        end
        peakLocationRhoXX = peakRhoXX;
        peakLocationRhoXXScale = max(abs(rhoXX(1,:)));
        peakLocationPhaseDefect = abs(peakRhoXX)/ ...
            max(peakLocationRhoXXScale,realmin);
        peakLocationResponseCondition = abs(peakRhoXXX)/ ...
            max(curvatureScale,realmin);
    end
else
    c_r = -trackedPeakU1-c_l*trackedPeakX;
end

rhoX = flow.source;
wallOmegaPeak = details.trackedWallPeak;
wallOmegaPeakIndex = details.wallOmegaPeakIndex;
if isstruct(phaseContext) && isfield(phaseContext,'generatorBaseRhs')
    % The candidate's actual discrete RHS is the fixed signed-basis sum.
    % Net velocities remain continuum diagnostics and do not define its CFL.
    baseRhs=phaseContext.generatorBaseRhs;
    transportU1=flow.u1+c_l*ops.X+c_r;
    transportU2=flow.u2+c_l*ops.Y;
    details.generatorCflRate=phaseContext.generatorCflRate;
    details.generatorReferenceVelocity=phaseContext.generatorReferenceVelocity;
    details.generatorTranslationRate=phaseContext.generatorTranslationRate;
    details.generatorRootCondition=phaseContext.generatorRootCondition;
else
    [~,baseRhs,transportU1,transportU2] = ipm.evolve.assembleRhs( ...
        rho,flow.u1,flow.u2,c_l,c_l,0,c_r,ops,'open');
end
if strcmp(r.lengthGauge,'transport_anchor')
    transportAnchorX = r.transportAnchorX;
    transportAnchorVelocity = interp1( ...
        ops.x,transportU1(1,:),transportAnchorX,'linear');
else
    transportAnchorX = NaN;
    transportAnchorVelocity = NaN;
end

if needsStrain
    strain = u1X(1,r.originIndex);
    strainScale = max(abs(u1X),[],'all');
    strainCondition = abs(strain)/max(strainScale,realmin);
else
    strain = NaN;
    strainScale = NaN;
    strainCondition = NaN;
end
omegaGaugeAnchorX = NaN;
omegaGaugeAnchorValue = NaN;
omegaGaugeAnchorCondition = NaN;
omegaGaugeAnchorRateResidual = NaN;
omegaGaugeReferenceAnchorValue = NaN;
omegaGaugeWindowCenter = NaN;
omegaGaugeWindowRadius = NaN;
omegaGaugeWindowValue = NaN;
omegaGaugeWindowCondition = NaN;
omegaGaugeWindowEnergy = NaN;
omegaGaugeWindowRateResidual = NaN;
omegaGaugeReferenceWindowValue = NaN;
omegaGaugeWindowMomentOrder = NaN;
omegaGaugeWindowMoment = NaN;
omegaGaugeWindowSupportPoints = NaN;
omegaGaugeWindowRateScale = NaN;
omegaGaugeTemplateCenter = NaN;
omegaGaugeTemplateRadius = NaN;
omegaGaugeTemplateSupportPoints = NaN;
omegaGaugeTemplateValue = NaN;
omegaGaugeTemplateCondition = NaN;
omegaGaugeTemplateRateResidual = NaN;
omegaGaugeTemplateRateScale = NaN;
omegaGaugeReferenceTemplateValue = NaN;
omegaGaugeBulkCenterX = NaN;
omegaGaugeBulkCenterY = NaN;
omegaGaugeBulkRadius = NaN;
omegaGaugeBulkSupportPointsX = NaN;
omegaGaugeBulkSupportPointsY = NaN;
omegaGaugeBulkValue = NaN;
omegaGaugeBulkCondition = NaN;
omegaGaugeBulkEnergy = NaN;
omegaGaugeBulkRateResidual = NaN;
omegaGaugeBulkRateScale = NaN;
omegaGaugeReferenceBulkValue = NaN;
omegaGaugeEnergyX = NaN;
omegaGaugeEnergyY = NaN;
omegaGaugeForcingX = NaN;
omegaGaugeForcingY = NaN;
omegaGaugeForcing = NaN;
omegaGaugeAbsoluteForcing = NaN;
omegaGaugeForcingCancellationRatio = NaN;
omegaGaugeQuadraticPeakX = NaN;
omegaGaugeQuadraticPeakValue = NaN;
omegaGaugeQuadraticForcing = NaN;
omegaGaugeQuadraticCurvature = NaN;
omegaGaugeQuadraticStationarityResidual = NaN;
omegaGaugeQuadraticNormalizedStationarityResidual = NaN;
omegaGaugeQuadraticCurvatureCondition = NaN;
omegaGaugeQuadraticWeightCondition = NaN;
omegaGaugeQuadraticActiveNodeMargin = NaN;
omegaGaugeQuadraticActiveSelectionMargin = NaN;
omegaGaugeQuadraticStencilCenterX = NaN;
omegaGaugeQuadraticVertexInteriorFraction = NaN;
omegaGaugeQuadraticRelativeSubgridLift = NaN;
omegaGaugeQuadraticRateResidual = NaN;
omegaGaugeQuadraticRateScale = NaN;
switch r.cOmegaGauge
    case 'wall_omega_peak'
        baseRhsX = baseRhs*ops.Dx';
        if wallOmegaPeak < r.cOmegaGaugeFloor
            error('ipm:DegenerateNormalization', ...
                ['The wall-omega peak c_omega gauge denominator became ' ...
                'too small.']);
        end
        gaugeForcing = baseRhsX(1,wallOmegaPeakIndex);
        gaugeEnergy = wallOmegaPeak;
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaPeakRateResidual = gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeResidual = omegaPeakRateResidual / ...
            max([abs(gaugeForcing),abs(c_omega*gaugeEnergy),realmin]);
        strainRateResidual = NaN;
    case {'wall_omega_quadratic_peak','wall_omega_tracked_spline_peak'}
        trackedSpline=strcmp(r.cOmegaGauge,'wall_omega_tracked_spline_peak');
        if trackedSpline
            peak=details.trackedPeakFunctional;
        else
            peak = ipm.evolve.quadraticPeakFunctional( ...
                rhoX(1,:),ops.x,wallOmegaPeakIndex);
        end
        if ~trackedSpline && any(abs(ops.x(peak.indices)-r.pinX) > ...
                r.peakTrackingHalfWidth)
            error('ipm:DegenerateNormalization', ...
                ['The quadratic wall-omega peak stencil crosses the ' ...
                'registered tracking-window boundary.']);
        end
        sourceScale = max(abs(rhoX),[],'all');
        peakCondition = peak.value/max(sourceScale,realmin);
        if peak.value < r.cOmegaGaugeFloor || ...
                peakCondition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The quadratic wall-omega peak became too small ' ...
                'relative to the full R_X field.']);
        end
        if isstruct(phaseContext) && isfield(phaseContext,'generatorBaseRhs')
            [forcingHigh,forcingLow]=ipm.field.compensatedMultiply( ...
                ops.Dx(peak.indices,:),baseRhs(1,:)');
            [forcingHigh,forcingLow]=ipm.field.compensatedMultiply( ...
                sparse(peak.weights),forcingHigh,forcingLow);
            gaugeForcing=forcingHigh+forcingLow;
        else
            stencilForcing = baseRhs(1,:)*ops.Dx(peak.indices,:)';
            gaugeForcing = sum(peak.weights.*stencilForcing);
        end
        gaugeEnergy = peak.value;
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The quadratic wall-omega peak forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeQuadraticRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeQuadraticRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),realmin]);
        omegaGaugeResidual = omegaGaugeQuadraticRateResidual / ...
            omegaGaugeQuadraticRateScale;
        omegaPeakRateResidual = omegaGaugeQuadraticRateResidual;
        omegaGaugeQuadraticPeakX = peak.x;
        omegaGaugeQuadraticPeakValue = peak.value;
        omegaGaugeQuadraticForcing = gaugeForcing;
        omegaGaugeQuadraticCurvature = peak.curvature;
        omegaGaugeQuadraticStationarityResidual = ...
            peak.stationarityResidual;
        omegaGaugeQuadraticNormalizedStationarityResidual = ...
            peak.normalizedStationarityResidual;
        omegaGaugeQuadraticCurvatureCondition = ...
            peak.curvatureCondition;
        omegaGaugeQuadraticWeightCondition = peak.weightCondition;
        omegaGaugeQuadraticActiveNodeMargin = peak.activeNodeMargin;
        trackingMask = abs(ops.x-r.pinX) <= r.peakTrackingHalfWidth;
        competitorMask = trackingMask;
        competitorMask(wallOmegaPeakIndex) = false;
        selectionCompetitors = rhoX(1,competitorMask);
        if isempty(selectionCompetitors)
            error('ipm:DegenerateNormalization', ...
                'The quadratic wall-omega tracking window has no peer node.');
        end
        omegaGaugeQuadraticActiveSelectionMargin = ( ...
            rhoX(1,wallOmegaPeakIndex)-max(selectionCompetitors)) / ...
            peak.value;
        selectionRoundoffTolerance = 512*eps(1)*max( ...
            max(abs(rhoX(1,trackingMask))),realmin)/peak.value;
        if ~trackedSpline && omegaGaugeQuadraticActiveSelectionMargin < ...
                -selectionRoundoffTolerance
            error('ipm:DegenerateNormalization', ...
                ['The quadratic wall-omega stencil center is not the ' ...
                'active maximum in the registered tracking window.']);
        end
        omegaGaugeQuadraticStencilCenterX = ...
            ops.x(wallOmegaPeakIndex);
        omegaGaugeQuadraticVertexInteriorFraction = ...
            peak.vertexInteriorFraction;
        omegaGaugeQuadraticRelativeSubgridLift = ...
            peak.relativeSubgridLift;
        details.omegaGaugePeakX = peak.x;
        details.omegaGaugePeakCondition = peakCondition;
        strainRateResidual = NaN;
    case 'gradient_energy'
        integrationWeights = ops.integrationWeights;
        rhoY = ops.Dy*rho;
        baseRhsX = baseRhs*ops.Dx';
        baseRhsY = ops.Dy*baseRhs;
        omegaGaugeEnergyX = sum( ...
            rhoX.^2.*integrationWeights,'all');
        omegaGaugeEnergyY = sum( ...
            rhoY.^2.*integrationWeights,'all');
        omegaGaugeForcingX = sum( ...
            rhoX.*baseRhsX.*integrationWeights,'all');
        omegaGaugeForcingY = sum( ...
            rhoY.*baseRhsY.*integrationWeights,'all');
        gaugeEnergy = sum((rhoX.^2+rhoY.^2).*integrationWeights,'all');
        gaugeForcing = sum((rhoX.*baseRhsX+rhoY.*baseRhsY) .* ...
            integrationWeights,'all');
        omegaGaugeAbsoluteForcing = sum(abs( ...
            rhoX.*baseRhsX+rhoY.*baseRhsY).*integrationWeights,'all');
        finiteGaugeValues = [omegaGaugeEnergyX,omegaGaugeEnergyY, ...
            gaugeEnergy,omegaGaugeForcingX,omegaGaugeForcingY, ...
            gaugeForcing,omegaGaugeAbsoluteForcing];
        if any(~isfinite(finiteGaugeValues))
            error('ipm:DegenerateNormalization', ...
                ['The full-gradient c_omega gauge energy or forcing ' ...
                'diagnostic is not finite.']);
        end
        if gaugeEnergy < r.cOmegaGaugeFloor
            error('ipm:DegenerateNormalization', ...
                ['The full-gradient c_omega gauge energy became too ' ...
                'small.']);
        end
        omegaGaugeForcing = gaugeForcing;
        omegaGaugeForcingCancellationRatio = abs(gaugeForcing) / ...
            max(omegaGaugeAbsoluteForcing,realmin);
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeResidual = (gaugeForcing+c_omega*gaugeEnergy) / ...
            max([abs(gaugeForcing),abs(c_omega*gaugeEnergy),realmin]);
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'strain_point'
        if strainCondition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The point-strain c_omega gauge is singular: ' ...
                'd_x u_1(0,0) is too small relative to the strain field.']);
        end
        [baseU1,~] = ipm.field.velocity(baseRhs,ops);
        baseStrain = baseU1*ipm.mesh.oddDx(ops)';
        gaugeForcing = baseStrain(1,r.originIndex);
        gaugeEnergy = strain;
        c_omega = -gaugeForcing/strain;
        strainRateResidual = gaugeForcing+c_omega*strain;
        omegaGaugeResidual = strainRateResidual / ...
            max([abs(gaugeForcing),abs(c_omega*strain),realmin]);
        omegaPeakRateResidual = NaN;
    case {'anchor_wall_window_l2','anchor_wall_window_l4'}
        baseRhsX = baseRhs*ops.Dx';
        if strcmp(r.cOmegaGauge,'anchor_wall_window_l4')
            momentOrder = 4;
        else
            momentOrder = 2;
        end
        window = ipm.evolve.wallWindowMoment( ...
            rhoX(1,:),ops,r,momentOrder);
        if ~isfield(r,'referenceAnchorWallWindowValue') || ...
                ~isnumeric(r.referenceAnchorWallWindowValue) || ...
                ~isreal(r.referenceAnchorWallWindowValue) || ...
                ~isscalar(r.referenceAnchorWallWindowValue) || ...
                ~isfinite(r.referenceAnchorWallWindowValue) || ...
                r.referenceAnchorWallWindowValue <= 0
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall window gauge lacks a positive ' ...
                'frozen reference value.']);
        end
        omegaGaugeWindowCenter = window.center;
        omegaGaugeWindowRadius = window.radius;
        omegaGaugeWindowValue = window.value;
        omegaGaugeWindowCondition = window.condition;
        omegaGaugeWindowEnergy = window.energy;
        omegaGaugeWindowMomentOrder = window.order;
        omegaGaugeWindowMoment = window.moment;
        omegaGaugeWindowSupportPoints = window.supportPoints;
        omegaGaugeReferenceWindowValue = ...
            r.referenceAnchorWallWindowValue;
        gaugeEnergy = window.moment;
        gaugeForcing = sum(window.normalizedWeights .* ...
            rhoX(1,:).^(momentOrder-1).*baseRhsX(1,:));
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The anchor-wall moment forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeWindowRateResidual = gaugeForcing + ...
            c_omega*gaugeEnergy;
        omegaGaugeWindowRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy), ...
            strainScale*gaugeEnergy]);
        omegaGaugeResidual = omegaGaugeWindowRateResidual / ...
            max(omegaGaugeWindowRateScale,realmin);
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'outer_wall_density_window_l2'
        window = ipm.evolve.outerWallDensityWindow( ...
            rho(1,:),ops,r);
        if ~isfield(r,'referenceOuterWallDensityWindowValue') || ...
                ~isnumeric(r.referenceOuterWallDensityWindowValue) || ...
                ~isreal(r.referenceOuterWallDensityWindowValue) || ...
                ~isscalar(r.referenceOuterWallDensityWindowValue) || ...
                ~isfinite(r.referenceOuterWallDensityWindowValue) || ...
                r.referenceOuterWallDensityWindowValue <= 0
            error('ipm:DegenerateNormalization', ...
                ['The outer wall-density gauge lacks a positive ' ...
                'frozen reference value.']);
        end
        gaugeEnergy = window.moment;
        gaugeForcing = sum(window.normalizedWeights .* ...
            rho(1,:).*baseRhs(1,:));
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The outer wall-density forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeWindowRateResidual = gaugeForcing + ...
            c_omega*gaugeEnergy;
        omegaGaugeWindowRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),realmin]);
        omegaGaugeResidual = omegaGaugeWindowRateResidual / ...
            omegaGaugeWindowRateScale;
        omegaGaugeWindowCenter = window.center;
        omegaGaugeWindowRadius = window.radius;
        omegaGaugeWindowValue = window.value;
        omegaGaugeWindowCondition = window.condition;
        omegaGaugeWindowEnergy = window.energy;
        omegaGaugeWindowMomentOrder = window.order;
        omegaGaugeWindowMoment = window.moment;
        omegaGaugeWindowSupportPoints = window.supportPoints;
        omegaGaugeReferenceWindowValue = ...
            r.referenceOuterWallDensityWindowValue;
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'anchor_wall_template_projection'
        baseRhsX = baseRhs*ops.Dx';
        geometry = ipm.evolve.wallWindowGeometry(ops,r);
        if ~isfield(r,'referenceAnchorWallTemplate') || ...
                ~isnumeric(r.referenceAnchorWallTemplate) || ...
                ~isreal(r.referenceAnchorWallTemplate) || ...
                ~isequal(size(r.referenceAnchorWallTemplate),size(ops.x)) || ...
                any(~isfinite(r.referenceAnchorWallTemplate),'all') || ...
                ~isfield(r,'referenceAnchorWallTemplateProjection') || ...
                ~isnumeric(r.referenceAnchorWallTemplateProjection) || ...
                ~isreal(r.referenceAnchorWallTemplateProjection) || ...
                ~isscalar(r.referenceAnchorWallTemplateProjection) || ...
                ~isfinite(r.referenceAnchorWallTemplateProjection) || ...
                r.referenceAnchorWallTemplateProjection <= 0
            error('ipm:DegenerateNormalization', ...
                'The anchor-wall template gauge lacks a valid frozen template.');
        end
        template = r.referenceAnchorWallTemplate;
        weights = geometry.normalizedWeights;
        templateNorm = sqrt(sum(weights.*template.^2));
        currentNorm = sqrt(sum(weights.*rhoX(1,:).^2));
        projection = sum(weights.*rhoX(1,:).*template);
        condition = abs(projection)/max(templateNorm*currentNorm,realmin);
        if ~all(isfinite([templateNorm,currentNorm,projection,condition])) || ...
                abs(projection) < r.cOmegaGaugeFloor || ...
                projection*r.referenceAnchorWallTemplateProjection <= 0 || ...
                condition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall template projection became singular or ' ...
                'lost the sign of its frozen reference.']);
        end
        gaugeEnergy = projection;
        gaugeForcing = sum(weights.*baseRhsX(1,:).*template);
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The anchor-wall template forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeTemplateRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeTemplateRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),strainScale*abs(gaugeEnergy)]);
        omegaGaugeResidual = omegaGaugeTemplateRateResidual / ...
            max(omegaGaugeTemplateRateScale,realmin);
        omegaGaugeTemplateCenter = geometry.center;
        omegaGaugeTemplateRadius = geometry.radius;
        omegaGaugeTemplateSupportPoints = geometry.supportPoints;
        omegaGaugeTemplateValue = projection;
        omegaGaugeTemplateCondition = condition;
        omegaGaugeReferenceTemplateValue = ...
            r.referenceAnchorWallTemplateProjection;
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'anchor_bulk_gradient_l2'
        rhoY = ops.Dy*rho;
        baseRhsX = baseRhs*ops.Dx';
        baseRhsY = ops.Dy*baseRhs;
        bulk = ipm.evolve.bulkGradientWindow(rhoX,rhoY,ops,r);
        if ~isfield(r,'referenceAnchorBulkGradientL2') || ...
                ~isnumeric(r.referenceAnchorBulkGradientL2) || ...
                ~isreal(r.referenceAnchorBulkGradientL2) || ...
                ~isscalar(r.referenceAnchorBulkGradientL2) || ...
                ~isfinite(r.referenceAnchorBulkGradientL2) || ...
                r.referenceAnchorBulkGradientL2 <= 0
            error('ipm:DegenerateNormalization', ...
                ['The anchor-bulk gradient L2 gauge lacks a positive ' ...
                'frozen reference value.']);
        end
        gaugeEnergy = bulk.energy;
        gaugeForcing = sum(bulk.normalizedWeights .* ...
            (rhoX.*baseRhsX+rhoY.*baseRhsY),'all');
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                'The anchor-bulk gradient forcing is not finite.');
        end
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeBulkRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        omegaGaugeBulkRateScale = max([abs(gaugeForcing), ...
            abs(c_omega*gaugeEnergy),strainScale*gaugeEnergy]);
        omegaGaugeResidual = omegaGaugeBulkRateResidual / ...
            max(omegaGaugeBulkRateScale,realmin);
        omegaGaugeBulkCenterX = bulk.centerX;
        omegaGaugeBulkCenterY = bulk.centerY;
        omegaGaugeBulkRadius = bulk.radius;
        omegaGaugeBulkSupportPointsX = bulk.supportPointsX;
        omegaGaugeBulkSupportPointsY = bulk.supportPointsY;
        omegaGaugeBulkValue = bulk.value;
        omegaGaugeBulkCondition = bulk.condition;
        omegaGaugeBulkEnergy = bulk.energy;
        omegaGaugeReferenceBulkValue = ...
            r.referenceAnchorBulkGradientL2;
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'anchor_wall_slope'
        baseRhsX = baseRhs*ops.Dx';
        omegaGaugeAnchorX = r.transportAnchorX;
        % Evaluate every member of the anchor constraint with one common
        % fixed-coordinate linear interpolation functional.
        anchorSamples = interp1(ops.x,[rhoX(1,:);baseRhsX(1,:)]', ...
            omegaGaugeAnchorX,'linear');
        omegaGaugeAnchorValue = anchorSamples(1);
        gaugeForcing = anchorSamples(2);
        wallSlopeScale = max(abs(rhoX(1,:)));
        omegaGaugeAnchorCondition = abs(omegaGaugeAnchorValue) / ...
            max(wallSlopeScale,realmin);
        if abs(omegaGaugeAnchorValue) < r.cOmegaGaugeFloor || ...
                omegaGaugeAnchorCondition < ...
                r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall c_omega gauge slope became too small ' ...
                'relative to the wall R_X profile.']);
        end
        if ~isfield(r,'referenceAnchorWallSlope') || ...
                omegaGaugeAnchorValue*r.referenceAnchorWallSlope <= 0
            error('ipm:DegenerateNormalization', ...
                ['The anchor-wall c_omega gauge lost the sign of its ' ...
                'frozen reference slope.']);
        end
        omegaGaugeReferenceAnchorValue = r.referenceAnchorWallSlope;
        gaugeEnergy = omegaGaugeAnchorValue;
        c_omega = -gaugeForcing/gaugeEnergy;
        omegaGaugeAnchorRateResidual = ...
            gaugeForcing+c_omega*gaugeEnergy;
        anchorSlopeRateScale = max([ ...
            abs(gaugeForcing),abs(c_omega*gaugeEnergy), ...
            strainScale*abs(gaugeEnergy)]);
        omegaGaugeResidual = omegaGaugeAnchorRateResidual / ...
            max(anchorSlopeRateScale,realmin);
        omegaPeakRateResidual = NaN;
        strainRateResidual = NaN;
    case 'none'
        c_omega = 0;
        gaugeEnergy = sum(rho.^2.*ops.integrationWeights,'all');
        omegaGaugeResidual = 0;
        % C_omega=1 is prescribed by this identity/material-amplitude
        % gauge, so its selected zero rate is not empirical convergence
        % evidence.  Independently evaluate the local quadratic wall-peak
        % functional.  Here F/P is the actual peak log-rate and -F/P is
        % the counterfactual c_omega that would freeze the peak.  Neither
        % quantity feeds back into the selected zero rate.
        omegaPeakRateResidual = NaN;
        try
        peak = ipm.evolve.quadraticPeakFunctional( ...
            rhoX(1,:),ops.x,wallOmegaPeakIndex);
        if any(abs(ops.x(peak.indices)-r.pinX) > ...
                r.peakTrackingHalfWidth)
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic stencil crosses ' ...
                'the registered tracking-window boundary.']);
        end
        sourceScale = max(abs(rhoX),[],'all');
        peakCondition = peak.value/max(sourceScale,realmin);
        if peak.value < r.cOmegaGaugeFloor || ...
                peakCondition < r.pointStrainConditionFloor
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic peak became too ' ...
                'small relative to the full R_X field.']);
        end
        stencilForcing = baseRhs(1,:)*ops.Dx(peak.indices,:)';
        gaugeForcing = sum(peak.weights.*stencilForcing);
        if ~isfinite(gaugeForcing)
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic forcing is not ' ...
                'finite.']);
        end
        omegaGaugeQuadraticRateResidual = gaugeForcing;
        omegaGaugeQuadraticRateScale = max(abs(gaugeForcing),realmin);
        omegaGaugeQuadraticPeakX = peak.x;
        omegaGaugeQuadraticPeakValue = peak.value;
        omegaGaugeQuadraticForcing = gaugeForcing;
        omegaGaugeQuadraticCurvature = peak.curvature;
        omegaGaugeQuadraticStationarityResidual = ...
            peak.stationarityResidual;
        omegaGaugeQuadraticNormalizedStationarityResidual = ...
            peak.normalizedStationarityResidual;
        omegaGaugeQuadraticCurvatureCondition = peak.curvatureCondition;
        omegaGaugeQuadraticWeightCondition = peak.weightCondition;
        omegaGaugeQuadraticActiveNodeMargin = peak.activeNodeMargin;
        trackingMask = abs(ops.x-r.pinX) <= r.peakTrackingHalfWidth;
        competitorMask = trackingMask;
        competitorMask(wallOmegaPeakIndex) = false;
        selectionCompetitors = rhoX(1,competitorMask);
        if isempty(selectionCompetitors)
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic window has no ' ...
                'peer node.']);
        end
        omegaGaugeQuadraticActiveSelectionMargin = ( ...
            rhoX(1,wallOmegaPeakIndex)-max(selectionCompetitors)) / ...
            peak.value;
        selectionRoundoffTolerance = 512*eps(1)*max( ...
            max(abs(rhoX(1,trackingMask))),realmin)/peak.value;
        if omegaGaugeQuadraticActiveSelectionMargin < ...
                -selectionRoundoffTolerance
            error('ipm:DegenerateNormalization', ...
                ['The identity-gauge quadratic diagnostic center is not ' ...
                'the active maximum in the registered window.']);
        end
        omegaGaugeQuadraticStencilCenterX = ...
            ops.x(wallOmegaPeakIndex);
        omegaGaugeQuadraticVertexInteriorFraction = ...
            peak.vertexInteriorFraction;
        omegaGaugeQuadraticRelativeSubgridLift = peak.relativeSubgridLift;
        details.omegaGaugePeakX = peak.x;
        details.omegaGaugePeakCondition = peakCondition;
        omegaPeakRateResidual = gaugeForcing;
        catch exception
            % A diagnostic peak may disappear in a generic identity-gauge
            % calculation.  That is recorded by absent/NaN telemetry and
            % must not turn an otherwise valid prescribed gauge into a
            % hidden normalization failure.
            if ~strcmp(exception.identifier,'ipm:DegenerateNormalization')
                rethrow(exception)
            end
        end
        strainRateResidual = NaN;
    case 'anchor_wall_strain'
        error('ipm:RetiredCOmegaGauge', ...
            ['anchor_wall_strain is not an exact discrete normalization ' ...
            'and is retired from the current runtime contract.']);
    otherwise
        error('ipm:BadCOmegaGaugeRuntime', ...
            'Unknown runtime c_omega gauge: %s.',r.cOmegaGauge);
end

if coupledWidthPeakGauge
    % The relative-height width is homogeneous in R_X: c_omega*R must
    % cancel from its directional derivative. Audit the actual complete
    % RHS with the same tracking mask and crossings used by trackFeatures.
    completeRhs = baseRhs+c_omega*rho;
    audited = ipm.diagnostics.trackFeatures( ...
        rho,completeRhs,flow,ops);
    lengthGaugeResidual = audited.wallWidthRate;
    [~,levelIndex] = max(r.adaptiveLevels);
    width = width_geometry(flow.source(1,:),ops.x, ...
        details.wallOmegaPeakIndex,r.adaptiveLevels(levelIndex), ...
        r.pinX,r.peakTrackingHalfWidth);
    widthMagnitude = width_functional_magnitude( ...
        width,abs(baseRhs(1,:))+abs(c_omega*rho(1,:)),ops.Dx);
    gamma512 = (512*eps)/(1-512*eps);
    widthAuditTolerance = min( ...
        1e-8*details.wallWidthCurrent,max( ...
            1e-11*max([1,abs(details.wallWidthRate), ...
            abs(c_l*details.wallWidthCurrent)]), ...
            gamma512*widthMagnitude));
    if ~isempty(phaseContext) && ...
            isfield(phaseContext,'widthRootAudit')
        % A local response certificate may exceed the final-product
        % gamma bound. Carry its measured floor into this post-audit only
        % after verifying the same wall operator and amplitude cancellation
        % against the complete two-dimensional RHS.
        columns = phaseContext.columns;
        localWall = phaseContext.widthRootWallRhs;
        wallDifferenceBound = width_functional_magnitude( ...
            width,localWall-baseRhs(1,columns), ...
            ops.Dx(:,columns));
        baseSamples = baseRhs(1,:)* ...
            ops.Dx(width.rateNodes,:)';
        sampledRate = zeros(1,ops.nx);
        sampledRate(width.rateNodes) = baseSamples;
        baseWidthRate = width_rate_from_samples(width,sampledRate);
        if wallDifferenceBound > gamma512*widthMagnitude || ...
                abs(baseWidthRate-phaseContext.widthRootRate) > ...
                wallDifferenceBound+gamma512*widthMagnitude || ...
                abs(lengthGaugeResidual-baseWidthRate) > ...
                gamma512*widthMagnitude
            error('ipm:WidthGaugeLocalFullMismatch', ...
                ['The certified local width root is inconsistent with ' ...
                'the complete RHS or amplitude cancellation.']);
        end
        widthAuditTolerance = min( ...
            1e-8*details.wallWidthCurrent,max( ...
            widthAuditTolerance, ...
            phaseContext.widthRootAudit.phaseFloor + ...
            gamma512*widthMagnitude));
    end
    if abs(lengthGaugeResidual) > widthAuditTolerance
        error('ipm:WidthGaugeFullRhsMismatch', ...
            ['The complete RHS has wall-width rate %.3e, exceeding ' ...
            'the conditioned audit tolerance %.3e.'], ...
            lengthGaugeResidual,widthAuditTolerance);
    end
end

if coupledCurvaturePeakGauge
    % Check both locked wall-gradient observables on the complete RHS.
    fullWallRate = (baseRhs(1,:)+c_omega*rho(1,:))*ops.Dx';
    peak = ipm.evolve.quadraticPeakFunctional( ...
        rhoX(1,:),ops.x,details.wallOmegaPeakIndex);
    peakRate = peak.weights*fullWallRate(phaseNodes)';
    phaseRate = phaseWeights*fullWallRate(phaseNodes)';
    peakPositionRate = phaseRate/phaseContext.peakK;
    lengthGaugeResidual = phaseContext.curvatureWeights* ...
        fullWallRate(phaseContext.curvatureNodes)' + ...
        phaseContext.curvatureThird*peakPositionRate;
    curvatureScale = max(1,abs(phaseContext.curvature));
    if abs(peakRate) > 1e-8*max(1,peak.value) || ...
            abs(lengthGaugeResidual) > 1e-8*curvatureScale
        error('ipm:PeakCurvatureGaugeFullRhsMismatch', ...
            ['The complete RHS does not preserve the peak and ' ...
            'curvature (Pdot %.3e, Kdot %.3e, tolerances %.3e, %.3e).'], ...
            peakRate,lengthGaugeResidual, ...
            1e-8*max(1,peak.value),1e-8*curvatureScale);
    end
    details.peakCurvature = phaseContext.curvature;
    details.innerCurvatureWidth = sqrt( ...
        -2*peak.value/phaseContext.curvature);
    details.peakCurvatureRate = lengthGaugeResidual;
    details.peakHeightRate = peakRate;
    details.peakPositionTargetRate = ...
        -r.peakPositionRestoringRate*peak.x;
end

if coupledIntegralPeakGauge
    % Audit the integral and peak on the actual COMPLETE transport RHS.
    fullWallRhs = baseRhs(1,:)+c_omega*rho(1,:);
    fullWallRate = fullWallRhs*ops.Dx(phaseNodes,:)';
    peak=phaseContext.peakFunctional;
    peakRate = peak.weights*fullWallRate';
    peakPositionRate = (phaseWeights*fullWallRate')/phaseContext.peakK;
    integralNodes = phaseContext.integralNodes;
    integralRate = fullWallRhs*ops.Dx(integralNodes,:)';
    lengthGaugeResidual = ...
        phaseContext.integralWeights(integralNodes)*integralRate';
    if abs(lengthGaugeResidual) > ...
            1e-8*max(1,abs(phaseContext.integralValue)) || ...
            abs(peakRate) > 1e-8*max(1,peak.value) || ...
            abs(peakPositionRate+r.peakPositionRestoringRate*peak.x) > 1e-8
        error('ipm:WallIntegralGaugeFullRhsMismatch', ...
            ['The complete RHS fails the integral/peak/position gauges ' ...
            '(Idot %.3e, Pdot %.3e, position residual %.3e).'], ...
            lengthGaugeResidual,peakRate, ...
            peakPositionRate+r.peakPositionRestoringRate*peak.x);
    end
    if isfield(phaseContext,'generatorBaseRhs')
        % Audit the rates actually selected by the amplitude case against
        % this candidate's complete full-array RHS at the original stricter
        % 1e-9 root gate, in addition to the unchanged audit above.
        [rateHigh,rateLow]=ipm.field.compensatedMultiply( ...
            ops.Dx(phaseNodes,:),fullWallRhs');
        [pHigh,pLow]=ipm.field.compensatedMultiply( ...
            sparse(peak.weights),rateHigh,rateLow);
        [aHigh,aLow]=ipm.field.compensatedMultiply( ...
            sparse(phaseWeights),rateHigh,rateLow);
        [rateHigh,rateLow]=ipm.field.compensatedMultiply( ...
            ops.Dx(integralNodes,:),fullWallRhs');
        [iHigh,iLow]=ipm.field.compensatedMultiply( ...
            sparse(phaseContext.integralWeights(integralNodes)),rateHigh,rateLow);
        generatorResidual=[(iHigh+iLow)/phaseContext.integralValue; ...
            (pHigh+pLow)/peak.value; ...
            (aHigh+aLow)/phaseContext.peakK+r.peakPositionRestoringRate*peak.x];
        if any(~isfinite(generatorResidual)) || max(abs(generatorResidual))>1e-9
            error('ipm:UnresolvedWallIntegralGauge', ...
                'The complete fixed-generator RHS fails the unchanged root gate: [%g %g %g].', ...
                generatorResidual(1),generatorResidual(2),generatorResidual(3));
        end
        details.generatorRootResidual=max(abs(generatorResidual));
        details.generatorIntegralResidual=generatorResidual(1);
        details.generatorPeakResidual=generatorResidual(2);
        details.generatorPhaseResidual=generatorResidual(3);
    end
    details.wallIntegral = phaseContext.integralValue;
    details.wallIntegralRadius = r.wallIntegralRadius;
    details.wallIntegralSupportPoints = phaseContext.integralSupportPoints;
    details.wallIntegralCondition = phaseContext.integralCondition;
    details.wallIntegralRate = lengthGaugeResidual;
    details.wallIntegralRelativeDrift = ...
        (phaseContext.integralValue-r.referenceWallIntegral)/ ...
        max(abs(r.referenceWallIntegral),realmin);
    details.peakHeightRate = peakRate;
    details.peakPositionTargetRate = -r.peakPositionRestoringRate*peak.x;
    details.peakPositionRate = peakPositionRate;
end

if strcmp(r.lengthGauge,'omega_peak_location') || ...
        (isfield(r,'travelingWaveGauge') && ...
        strcmp(r.travelingWaveGauge,'gradient_peak') && ...
        ~coupledRestoringPeakGauge)
    if ~strcmp(r.lengthGauge,'omega_peak_location')
        peakLocationResponse = trackedPeakX*peakRhoXXX;
    end
    fullLengthContribution = c_l*( ...
        2*peakLocationRhoXX+peakLocationResponse);
    fullTranslationContribution = c_r*peakRhoXXX;
    fullAmplitudeContribution = c_omega*peakLocationRhoXX;
    peakLocationFullPhaseResidual = peakPhysicalRhsXX - ...
        fullLengthContribution-fullTranslationContribution + ...
        fullAmplitudeContribution;
    fullPhaseScale = max([abs(peakPhysicalRhsXX), ...
        abs(fullLengthContribution),abs(fullTranslationContribution), ...
        abs(fullAmplitudeContribution),realmin]);
    peakLocationNormalizedFullPhaseResidual = ...
        peakLocationFullPhaseResidual/fullPhaseScale;
end

details.strain = strain;
details.strainCondition = strainCondition;
details.c_lNominal = c_lNominal;
details.lengthScaleGain = 1;
details.widthCorrection = widthCorrection;
details.widthRateScale = abs(c_lNominal);
details.widthControlMode = widthControlMode;
details.peakExpansionRatio = peakExpansionRatio;
details.resolutionFactor = resolutionFactor;
details.safetyFactor = safetyFactor;
details.strainError = strain-r.strainTarget;
details.strainRateResidual = strainRateResidual;
details.omegaGaugeEnergy = gaugeEnergy;
details.omegaGaugeResidual = omegaGaugeResidual;
details.omegaGaugeEnergyX = omegaGaugeEnergyX;
details.omegaGaugeEnergyY = omegaGaugeEnergyY;
details.omegaGaugeForcingX = omegaGaugeForcingX;
details.omegaGaugeForcingY = omegaGaugeForcingY;
details.omegaGaugeForcing = omegaGaugeForcing;
details.omegaGaugeAbsoluteForcing = omegaGaugeAbsoluteForcing;
details.omegaGaugeForcingCancellationRatio = ...
    omegaGaugeForcingCancellationRatio;
details.omegaGaugeQuadraticPeakX = omegaGaugeQuadraticPeakX;
details.omegaGaugeQuadraticPeakValue = omegaGaugeQuadraticPeakValue;
details.omegaGaugeQuadraticForcing = omegaGaugeQuadraticForcing;
details.omegaGaugeQuadraticCurvature = omegaGaugeQuadraticCurvature;
details.omegaGaugeQuadraticStationarityResidual = ...
    omegaGaugeQuadraticStationarityResidual;
details.omegaGaugeQuadraticNormalizedStationarityResidual = ...
    omegaGaugeQuadraticNormalizedStationarityResidual;
details.omegaGaugeQuadraticCurvatureCondition = ...
    omegaGaugeQuadraticCurvatureCondition;
details.omegaGaugeQuadraticWeightCondition = ...
    omegaGaugeQuadraticWeightCondition;
details.omegaGaugeQuadraticActiveNodeMargin = ...
    omegaGaugeQuadraticActiveNodeMargin;
details.omegaGaugeQuadraticActiveSelectionMargin = ...
    omegaGaugeQuadraticActiveSelectionMargin;
details.omegaGaugeQuadraticStencilCenterX = ...
    omegaGaugeQuadraticStencilCenterX;
details.omegaGaugeQuadraticVertexInteriorFraction = ...
    omegaGaugeQuadraticVertexInteriorFraction;
details.omegaGaugeQuadraticRelativeSubgridLift = ...
    omegaGaugeQuadraticRelativeSubgridLift;
details.omegaGaugeQuadraticRateResidual = ...
    omegaGaugeQuadraticRateResidual;
details.omegaGaugeQuadraticRateScale = omegaGaugeQuadraticRateScale;
details.omegaPeakRateResidual = omegaPeakRateResidual;
details.omegaGaugeAnchorX = omegaGaugeAnchorX;
details.omegaGaugeAnchorValue = omegaGaugeAnchorValue;
details.omegaGaugeAnchorCondition = omegaGaugeAnchorCondition;
details.omegaGaugeAnchorRateResidual = omegaGaugeAnchorRateResidual;
details.omegaGaugeReferenceAnchorValue = ...
    omegaGaugeReferenceAnchorValue;
details.omegaGaugeWindowCenter = omegaGaugeWindowCenter;
details.omegaGaugeWindowRadius = omegaGaugeWindowRadius;
details.omegaGaugeWindowValue = omegaGaugeWindowValue;
details.omegaGaugeWindowCondition = omegaGaugeWindowCondition;
details.omegaGaugeWindowEnergy = omegaGaugeWindowEnergy;
details.omegaGaugeWindowRateResidual = omegaGaugeWindowRateResidual;
details.omegaGaugeReferenceWindowValue = ...
    omegaGaugeReferenceWindowValue;
details.omegaGaugeWindowMomentOrder = omegaGaugeWindowMomentOrder;
details.omegaGaugeWindowMoment = omegaGaugeWindowMoment;
details.omegaGaugeWindowSupportPoints = omegaGaugeWindowSupportPoints;
details.omegaGaugeWindowRateScale = omegaGaugeWindowRateScale;
details.omegaGaugeTemplateCenter = omegaGaugeTemplateCenter;
details.omegaGaugeTemplateRadius = omegaGaugeTemplateRadius;
details.omegaGaugeTemplateSupportPoints = ...
    omegaGaugeTemplateSupportPoints;
details.omegaGaugeTemplateValue = omegaGaugeTemplateValue;
details.omegaGaugeTemplateCondition = omegaGaugeTemplateCondition;
details.omegaGaugeTemplateRateResidual = ...
    omegaGaugeTemplateRateResidual;
details.omegaGaugeTemplateRateScale = omegaGaugeTemplateRateScale;
details.omegaGaugeReferenceTemplateValue = ...
    omegaGaugeReferenceTemplateValue;
details.omegaGaugeBulkCenterX = omegaGaugeBulkCenterX;
details.omegaGaugeBulkCenterY = omegaGaugeBulkCenterY;
details.omegaGaugeBulkRadius = omegaGaugeBulkRadius;
details.omegaGaugeBulkSupportPointsX = omegaGaugeBulkSupportPointsX;
details.omegaGaugeBulkSupportPointsY = omegaGaugeBulkSupportPointsY;
details.omegaGaugeBulkValue = omegaGaugeBulkValue;
details.omegaGaugeBulkCondition = omegaGaugeBulkCondition;
details.omegaGaugeBulkEnergy = omegaGaugeBulkEnergy;
details.omegaGaugeBulkRateResidual = omegaGaugeBulkRateResidual;
details.omegaGaugeBulkRateScale = omegaGaugeBulkRateScale;
details.omegaGaugeReferenceBulkValue = omegaGaugeReferenceBulkValue;
details.lengthGaugeStrain = lengthGaugeStrain;
details.lengthGaugeResidual = lengthGaugeResidual;
details.peakLocationRhoXX = peakLocationRhoXX;
details.peakLocationRhoXXScale = peakLocationRhoXXScale;
details.peakLocationPhaseDefect = peakLocationPhaseDefect;
details.peakLocationResponseCondition = peakLocationResponseCondition;
details.peakLocationFullPhaseResidual = ...
    peakLocationFullPhaseResidual;
details.peakLocationNormalizedFullPhaseResidual = ...
    peakLocationNormalizedFullPhaseResidual;
details.pinGaugeResidual = details.lengthGaugeResidual;
details.transportAnchorX = transportAnchorX;
details.transportAnchorVelocity = transportAnchorVelocity;
if strcmp(ops.symmetryMode,'double_odd_omega')
    details.pinTransportVelocity = transportU1(1,r.originIndex);
else
    details.pinTransportVelocity = ...
        interp1(ops.x,transportU1(1,:),r.pinX,'linear');
end
details.peakTransportVelocity = interp1( ...
    ops.x,transportU1(1,:),trackedPeakX,'linear');
if strcmp(ops.symmetryMode,'double_odd_omega')
    details.travelingWaveResidual = transportU1(1,r.originIndex);
elseif isfield(r,'travelingWaveGauge') && ...
        strcmp(r.travelingWaveGauge,'gradient_peak')
    wallRhsX = (baseRhs+c_omega*rho)*ops.Dx';
    fullSample = wallRhsX(1,phaseNodes);
    details.travelingWaveResidual = phaseWeights*fullSample';
    if coupledRestoringPeakGauge
        details.travelingWaveResidual = ...
            details.travelingWaveResidual-phaseContext.phaseTarget;
    end
    if ~isempty(phaseContext)
        [localPhase,localScale,~,localWallRhs] = discrete_peak_phase( ...
            c_r,rho,flow,ops,c_l,phaseNodes,phaseWeights,phaseContext);
        baseSample = baseRhs(1,:)*ops.Dx(phaseNodes,:)';
        fullBasePhase = phaseWeights*baseSample';
        fullWallRhs = baseRhs(1,phaseContext.columns);
        absoluteFunctional = abs(phaseWeights)* ...
            abs(ops.Dx(phaseNodes,phaseContext.columns));
        rhsDifferenceBound = absoluteFunctional* ...
            abs(localWallRhs-fullWallRhs)';
        gamma64 = 64*eps/(1-64*eps);
        productRoundingBound = gamma64*( ...
            absoluteFunctional*(abs(localWallRhs)+abs(fullWallRhs))' + ...
            sum(abs(phaseWeights).* ...
            (abs(localWallRhs* ...
            ops.Dx(phaseNodes,phaseContext.columns)')+ ...
            abs(baseSample))));
        auditScale = max([1,localScale, ...
            sum(abs(phaseWeights).*abs(baseSample))]);
        if abs(localPhase-fullBasePhase) > ...
                1e-11*auditScale
            error('ipm:TravelingWaveLocalPhaseMismatch', ...
                ['The local wall phase differs from the full transport ' ...
                'phase by %.3e (scale %.3e, RHS difference bound %.3e, ' ...
                'product rounding bound %.3e, max wall RHS diff %.3e).'], ...
                abs(localPhase-fullBasePhase),auditScale, ...
                rhsDifferenceBound,productRoundingBound, ...
                max(abs(localWallRhs-fullWallRhs)));
        end
    end
else
    details.travelingWaveResidual = details.peakTransportVelocity;
end
end

function [cL,cR,nodes,weights,phaseContext] = ...
        pin_width_and_gradient_peak( ...
        rho,flow,ops,details,cLGuess,peakU1)
% Enforce both active discrete invariants on the same nonlinear wall RHS.
% The width and quadratic-peak functionals depend on only a few Dx rows.
% WENO still sees the full wall row; its vertical pass is rebuilt for every
% trial cL because u2+cL*Y changes with the trial rate.
r = ops.rescaling;
[~,levelIndex] = max(r.adaptiveLevels);
width = width_geometry(flow.source(1,:),ops.x, ...
    details.wallOmegaPeakIndex,r.adaptiveLevels(levelIndex), ...
    r.pinX,r.peakTrackingHalfWidth);
widthColumns = find(any(ops.Dx(width.rateNodes,:) ~= 0,1));
widthScale = max([1,abs(details.wallWidthRate), ...
    abs(cLGuess*details.wallWidthCurrent)]);
tolerance = 1e-11*widthScale;

[f0,cR,nodes,weights,phaseContext] = width_phase_at( ...
    cLGuess,rho,flow,ops,details,peakU1,width,widthColumns);
if abs(f0) <= tolerance
    cL = cLGuess;
    return;
end
probeStep = max(0.02,0.05*max(1,abs(cLGuess)));
[fProbe,probeCR,probeNodes,probeWeights,probeContext] = ...
    width_phase_at(cLGuess+probeStep,rho,flow,ops,details, ...
    peakU1,width,widthColumns);
if abs(fProbe) <= tolerance
    cL = cLGuess+probeStep;
    cR = probeCR;
    nodes = probeNodes;
    weights = probeWeights;
    phaseContext = probeContext;
    return;
end
slope = (fProbe-f0)/probeStep;
if ~isfinite(slope) || abs(slope) <= ...
        1e-10*max(details.wallWidthCurrent,realmin)
    error('ipm:DegenerateLengthGauge', ...
        'The discrete wall-width rate has negligible c_l response.');
end
direction = -sign(f0/slope);
step = min(max(probeStep,abs(f0/slope)), ...
    4*max(1,abs(cLGuess)));
bracketed = false;
for attempt = 1:8
    candidate = cLGuess+direction*step;
    [candidateRate,candidateCR,candidateNodes, ...
        candidateWeights,candidateContext] = width_phase_at( ...
        candidate,rho,flow,ops,details,peakU1,width,widthColumns);
    if abs(candidateRate) <= tolerance
        cL = candidate;
        cR = candidateCR;
        nodes = candidateNodes;
        weights = candidateWeights;
        phaseContext = candidateContext;
        return;
    end
    if f0*candidateRate <= 0
        bracketed = true;
        break;
    end
    step = 2*step;
end
if ~bracketed
    error('ipm:UnbracketedLengthGauge', ...
        'The discrete wall-width rate root was not bracketed.');
end
interval = sort([cLGuess,candidate]);
options = optimset('Display','off','TolX',1e-12, ...
    'MaxFunEvals',24,'MaxIter',24);
[cL,fRoot,exitFlag] = fzero( ...
    @(rate)width_phase_at( ...
    rate,rho,flow,ops,details,peakU1,width,widthColumns), ...
    interval,options);
[rootRate,cR,nodes,weights,phaseContext,rootWallRhs] = width_phase_at( ...
    cL,rho,flow,ops,details,peakU1,width,widthColumns);
if exitFlag <= 0 || ~isfinite(cL) || ~isfinite(rootRate)
    error('ipm:UnresolvedLengthGauge', ...
        'The discrete wall-width rate root did not converge.');
end
if max(abs(fRoot),abs(rootRate)) > tolerance
    % The selected-level crossing combines five Dx rows and a nested
    % nonlinear phase solve. Bound floating cancellation by the absolute
    % wall functional, then require a narrow sign-changing cL bracket.
    if isempty(phaseContext)
        sampledColumns = 1:ops.nx;
    else
        sampledColumns = phaseContext.columns;
    end
    magnitude = width_functional_magnitude( ...
        width,rootWallRhs,ops.Dx(:,sampledColumns));
    gamma512 = (512*eps)/(1-512*eps);
    roundingGuard = gamma512*magnitude;
    maxSpeedUncertainty = 1e-8*max(1,abs(cL));
    bracketHalfWidth = min(maxSpeedUncertainty, ...
        max(16*eps*max(1,abs(cL)), ...
        4*roundingGuard/max(abs(slope),realmin)));
    certified = false;
    for guardAttempt = 1:8
        [leftRate,~,~,~,~,leftWallRhs] = width_phase_at( ...
            cL-bracketHalfWidth,rho,flow,ops,details, ...
            peakU1,width,widthColumns);
        [rightRate,~,~,~,~,rightWallRhs] = width_phase_at( ...
            cL+bracketHalfWidth,rho,flow,ops,details, ...
            peakU1,width,widthColumns);
        endpointMagnitude = width_functional_magnitude( ...
            width,max(abs(leftWallRhs),abs(rightWallRhs)), ...
            ops.Dx(:,sampledColumns));
        roundingGuard = gamma512*max(magnitude,endpointMagnitude);
        speedBound = bracketHalfWidth + ...
            roundingGuard/max(abs(slope),realmin);
        certified = isfinite(roundingGuard) && ...
            max(abs(fRoot),abs(rootRate)) <= roundingGuard && ...
            leftRate*rightRate < 0 && ...
            min(abs(leftRate),abs(rightRate)) > roundingGuard && ...
            speedBound <= maxSpeedUncertainty;
        if certified || bracketHalfWidth >= maxSpeedUncertainty
            break;
        end
        bracketHalfWidth = min(2*bracketHalfWidth,maxSpeedUncertainty);
    end
    if ~certified
        % A nested WENO phase solve can leave a reproducible width-rate
        % floor just above the final crossing-product bound. Audit that
        % floor only after the original guard fails. Every cL sample
        % rebuilds its own vertical WENO contribution and resolves cR.
        slopeSign = sign(slope);
        widthFunction = @(rate)slopeSign*width_phase_at( ...
            rate,rho,flow,ops,details,peakU1,width,widthColumns);
        widthAudit = ipm.evolve.certifyGaugeRoot( ...
            widthFunction,cL,slopeSign*fRoot, ...
            slopeSign*rootRate,-abs(slope),roundingGuard, ...
            maxSpeedUncertainty);
        if ~widthAudit.accepted
            error('ipm:UnresolvedLengthGauge', ...
                ['The discrete wall-width rate root failed its ' ...
                'measured-noise bracket: fRoot %.3e, recomputed %.3e, ' ...
                'tolerance %.3e, guard %.3e, local floor %.3e, ' ...
                'half-width %.3e, signed endpoint rates ' ...
                '[%.3e, %.3e], local slope %.3e, speed bound %.3e, ' ...
                'c_l %.16g.'], ...
                fRoot,rootRate,tolerance,roundingGuard, ...
                widthAudit.phaseFloor,widthAudit.halfWidth, ...
                widthAudit.leftPhase,widthAudit.rightPhase, ...
                widthAudit.slope,widthAudit.speedBound,cL);
        end
        if ~isempty(phaseContext)
            phaseContext.widthRootAudit = widthAudit;
            phaseContext.widthRootRate = rootRate;
            phaseContext.widthRootWallRhs = rootWallRhs;
        end
    end
end
end

function [widthRate,cR,nodes,weights,phaseContext,wallRhs] = ...
        width_phase_at( ...
        cL,rho,flow,ops,details,peakU1,width,widthColumns)
[cR,nodes,weights,phaseContext] = pin_gradient_peak( ...
    rho,flow,ops,details,cL,peakU1,widthColumns);
if isempty(phaseContext)
    [~,baseRhs] = ipm.evolve.assembleRhs( ...
        rho,flow.u1,flow.u2,cL,cL,0,cR,ops,'open');
    wallRhs = baseRhs(1,:);
    sampledRate = wallRhs*ops.Dx(width.rateNodes,:)';
else
    [wallRhs,phaseContext.vertical] = ...
        ipm.field.weno5WallTransport( ...
        rho,flow.u1,flow.u2,cL,cR,ops, ...
        phaseContext.columns,phaseContext.vertical);
    sampledRate = wallRhs* ...
        ops.Dx(width.rateNodes,phaseContext.columns)';
end
featureRate = zeros(1,ops.nx);
featureRate(width.rateNodes) = sampledRate;
widthRate = width_rate_from_samples(width,featureRate);
end

function width = width_geometry(source,x,peakIndex,level,pinX,halfWidth)
feature = max(source(:)',0);
feature(abs(x-pinX) > halfWidth) = 0;
peak = feature(peakIndex);
threshold = level*peak;
leftNode = peakIndex;
while leftNode > 1 && feature(leftNode-1) >= threshold
    leftNode = leftNode-1;
end
rightNode = peakIndex;
while rightNode < numel(feature) && feature(rightNode+1) >= threshold
    rightNode = rightNode+1;
end
if leftNode == 1 || rightNode == numel(feature)
    error('ipm:DegenerateLengthGauge', ...
        'The tracked wall-omega selected-level component is not interior.');
end
[width.leftSlope,width.leftFraction] = width_crossing_data( ...
    leftNode-1,leftNode,feature,x,threshold);
[width.rightSlope,width.rightFraction] = width_crossing_data( ...
    rightNode,rightNode+1,feature,x,threshold);
width.leftNodes = [leftNode-1,leftNode];
width.rightNodes = [rightNode,rightNode+1];
width.peakIndex = peakIndex;
width.level = level;
width.rateNodes = unique([width.leftNodes,width.rightNodes,peakIndex]);
end

function [slope,fraction] = width_crossing_data( ...
        left,right,feature,x,threshold)
valueDifference = feature(right)-feature(left);
if abs(valueDifference) <= 100*eps(max(abs(feature)))
    error('ipm:DegenerateLengthGauge', ...
        'A wall-omega selected-level crossing has zero discrete slope.');
end
fraction = (threshold-feature(left))/valueDifference;
slope = valueDifference/(x(right)-x(left));
end

function widthRate = width_rate_from_samples(width,featureRate)
left = width.leftNodes;
right = width.rightNodes;
leftRate = featureRate(left(1))+width.leftFraction* ...
    (featureRate(left(2))-featureRate(left(1)));
rightRate = featureRate(right(1))+width.rightFraction* ...
    (featureRate(right(2))-featureRate(right(1)));
peakRate = featureRate(width.peakIndex);
leftVelocity = -(leftRate-width.level*peakRate)/width.leftSlope;
rightVelocity = -(rightRate-width.level*peakRate)/width.rightSlope;
widthRate = rightVelocity-leftVelocity;
end

function magnitude = width_functional_magnitude(width,wallRhs,Dx)
left = width.leftNodes;
right = width.rightNodes;
leftWeights = [1-width.leftFraction,width.leftFraction]/ ...
    abs(width.leftSlope);
rightWeights = [1-width.rightFraction,width.rightFraction]/ ...
    abs(width.rightSlope);
absoluteRows = leftWeights*abs(Dx(left,:)) + ...
    rightWeights*abs(Dx(right,:)) + ...
    width.level*(1/abs(width.leftSlope)+ ...
    1/abs(width.rightSlope))*abs(Dx(width.peakIndex,:));
magnitude = absoluteRows*abs(wallRhs)';
end

function [cR,nodes,weights,phaseContext] = pin_gradient_peak( ...
        rho,flow,ops,details,cL,peakU1,extraColumns)
% Solve the active quadratic vertex phase with the actual semi-discrete RHS.
center = details.wallOmegaPeakIndex;
if center <= 1 || center >= ops.nx
    error('ipm:DegenerateTravelingWaveGauge', ...
        'The tracked wall-gradient peak must have two interior neighbors.');
end
nodes = center-1:center+1;
a = details.trackedPeakX;
x = ops.x(nodes);
weights = zeros(1,3);
secondWeights = zeros(1,3);
for j = 1:3
    others = setdiff(1:3,j);
    denominator = (x(j)-x(others(1)))*(x(j)-x(others(2)));
    weights(j) = (2*a-x(others(1))-x(others(2)))/denominator;
    secondWeights(j) = 2/denominator;
end
q = flow.source(1,nodes);
curvature = secondWeights*q';
curvatureScale = sum(abs(secondWeights).*abs(q));
stationarity = weights*q';
stationarityScale = sum(abs(weights).*abs(q));
if ~isfinite(curvature) || curvature >= 0 || ...
        abs(curvature) <= 1e-8*max(curvatureScale,realmin) || ...
        abs(stationarity) > 1e-8*max(stationarityScale,realmin)
    error('ipm:DegenerateTravelingWaveGauge', ...
        'The tracked quadratic wall-gradient peak is not well conditioned.');
end

initial = -peakU1-cL*a;
phaseContext = [];
if strcmp(ops.transportScheme,'weno5_fd') && ...
        strcmp(ops.transportBoundaryMode,'open') && ...
        ~ipm.mesh.isQuadrant(ops)
    % Dx at three peak nodes samples at most nine wall RHS columns. Keep
    % every horizontal wall node in WENO's global LF scales, but cache
    % the cR-independent vertical derivative on these sampled columns.
    columns = find(any(ops.Dx(nodes,:) ~= 0,1));
    if nargin >= 7 && ~isempty(extraColumns)
        columns = union(columns,extraColumns);
    end
    phaseContext = struct('columns',columns,'vertical',[]);
end
[f0,scale0,phaseContext] = discrete_peak_phase( ...
    initial,rho,flow,ops,cL,nodes,weights,phaseContext);
tolerance = 1e-11;
if abs(f0) <= tolerance*max(scale0,1)
    cR = initial;
    return;
end
probeStep = max(0.02,0.05*max(1,abs(initial)));
[fProbe,probeScale] = discrete_peak_phase( ...
    initial+probeStep,rho,flow,ops,cL,nodes,weights,phaseContext);
if abs(fProbe) <= tolerance*max(probeScale,1)
    cR = initial+probeStep;
    return;
end
slope = (fProbe-f0)/probeStep;
if ~isfinite(slope) || abs(slope) <= ...
        1e-10*max(abs(curvature),realmin)
    error('ipm:DegenerateTravelingWaveGauge', ...
        'The discrete wall-gradient phase has negligible translation response.');
end
direction = -sign(f0/slope);
step = max(probeStep,abs(f0/slope));
step = min(step,4*max(1,abs(initial)));
bracketed = false;
for attempt = 1:8
    candidate = initial+direction*step;
    [fCandidate,candidateScale] = discrete_peak_phase( ...
        candidate,rho,flow,ops,cL,nodes,weights,phaseContext);
    if abs(fCandidate) <= tolerance*max(candidateScale,1)
        cR = candidate;
        return;
    end
    if f0*fCandidate <= 0
        bracketed = true;
        break;
    end
    step = 2*step;
end
if ~bracketed
    error('ipm:UnbracketedTravelingWaveGauge', ...
        'The discrete wall-gradient phase root was not bracketed.');
end
interval = sort([initial,candidate]);
options = optimset('Display','off','TolX',1e-12, ...
    'MaxFunEvals',24,'MaxIter',24);
[cR,fRoot,exitFlag] = fzero( ...
    @(rate)discrete_peak_phase( ...
    rate,rho,flow,ops,cL,nodes,weights,phaseContext), ...
    interval,options);
[rootPhase,rootScale,~,rootWallRhs] = discrete_peak_phase( ...
    cR,rho,flow,ops,cL,nodes,weights,phaseContext);
if exitFlag <= 0 || ~isfinite(cR) || ~isfinite(rootPhase)
    error('ipm:UnresolvedTravelingWaveGauge', ...
        'The discrete wall-gradient phase root did not converge.');
end
if max(abs(fRoot),abs(rootPhase)) > tolerance*max(rootScale,1)
    if isempty(phaseContext)
        error('ipm:UnresolvedTravelingWaveGauge', ...
            'The discrete wall-gradient phase root did not converge.');
    end
    % At fine spacing, the two discrete x derivatives cancel much more
    % strongly than the three final phase samples reveal. Use the absolute
    % magnitude of the full wall-RHS functional in a conservative 128-op
    % floating-point guard. Accept only a narrow, sign-certified bracket;
    % the resulting uncertainty in cR is bounded independently of fzero's
    % attainable residual. The ordinary tolerance path above is unchanged.
    absoluteFunctional = abs(weights)* ...
        abs(ops.Dx(nodes,phaseContext.columns));
    magnitude = absoluteFunctional*abs(rootWallRhs)';
    gamma128 = (128*eps)/(1-128*eps);
    roundingGuard = gamma128*magnitude;
    maxSpeedUncertainty = 1e-8*max(1,abs(cR));
    bracketHalfWidth = min(maxSpeedUncertainty, ...
        max(16*eps*max(1,abs(cR)), ...
        4*roundingGuard/max(abs(curvature),realmin)));
    certified = false;
    for guardAttempt = 1:8
        [leftPhase,~,~,leftWallRhs] = discrete_peak_phase( ...
            cR-bracketHalfWidth,rho,flow,ops,cL, ...
            nodes,weights,phaseContext);
        [rightPhase,~,~,rightWallRhs] = discrete_peak_phase( ...
            cR+bracketHalfWidth,rho,flow,ops,cL, ...
            nodes,weights,phaseContext);
        endpointMagnitude = absoluteFunctional* ...
            max(abs(leftWallRhs),abs(rightWallRhs))';
        roundingGuard = gamma128*max(magnitude,endpointMagnitude);
        speedBound = bracketHalfWidth + ...
            roundingGuard/max(abs(curvature),realmin);
        certified = isfinite(roundingGuard) && ...
            max(abs(fRoot),abs(rootPhase)) <= roundingGuard && ...
            leftPhase*rightPhase < 0 && ...
            min(abs(leftPhase),abs(rightPhase)) > roundingGuard && ...
            speedBound <= maxSpeedUncertainty;
        if certified || bracketHalfWidth >= maxSpeedUncertainty
            break;
        end
        bracketHalfWidth = min(2*bracketHalfWidth,maxSpeedUncertainty);
    end
    if ~certified
        % The wall WENO reconstruction itself can have a larger local
        % floating-point floor than the final Dx/vertex product bounded
        % above. Measure that floor from symmetric phase probes, and only
        % accept a sign-changing, well-conditioned bracket whose complete
        % speed uncertainty remains below the original 1e-8 budget.
        phaseFunction = @(rate)discrete_peak_phase( ...
            rate,rho,flow,ops,cL,nodes,weights,phaseContext);
        phaseAudit = ipm.evolve.certifyGaugeRoot( ...
            phaseFunction,cR,fRoot,rootPhase,curvature, ...
            roundingGuard,maxSpeedUncertainty);
        if ~phaseAudit.accepted
            error('ipm:UnresolvedTravelingWaveGauge', ...
                ['The discrete wall-gradient phase root failed its ' ...
                'measured-noise bracket: fRoot %.3e, rootScale %.3e, ' ...
                'guard %.3e, local floor %.3e, half-width %.3e, ' ...
                'endpoint phases [%.3e, %.3e], local slope %.3e, ' ...
                'speed bound %.3e, curvature %.3e.'], ...
                fRoot,rootScale,roundingGuard,phaseAudit.phaseFloor, ...
                phaseAudit.halfWidth,phaseAudit.leftPhase, ...
                phaseAudit.rightPhase,phaseAudit.slope, ...
                phaseAudit.speedBound,curvature);
        end
    end
end
end

function [phase,scale,phaseContext,wallRhs] = discrete_peak_phase( ...
        cR,rho,flow,ops,cL,nodes,weights,phaseContext)
wallRhs = [];
if isempty(phaseContext)
    [~,baseRhs] = ipm.evolve.assembleRhs( ...
        rho,flow.u1,flow.u2,cL,cL,0,cR,ops,'open');
    wallRhsX = baseRhs(1,:)*ops.Dx';
    sample = wallRhsX(nodes);
elseif isfield(phaseContext,'generatorBaseRhs')
    % The candidate has no optimized net-velocity wall operator. Compare
    % its selected wall subset to the same complete signed-basis assembly.
    columns=phaseContext.columns;
    wallRhs=phaseContext.generatorBaseRhs(1,columns);
    sample=wallRhs*ops.Dx(nodes,columns)';
else
    columns = phaseContext.columns;
    [wallRhs,phaseContext.vertical] = ...
        ipm.field.weno5WallTransport( ...
        rho,flow.u1,flow.u2,cL,cR,ops,columns, ...
        phaseContext.vertical);
    sample = wallRhs*ops.Dx(nodes,columns)';
end
phase = weights*sample';
scale = sum(abs(weights).*abs(sample));
end


function [cL,cOmega,cR,baseRhs,details,v1,v2]= ...
        practical_integral_centroid(rho,flow,ops,details)
% The practical mode preserves I and global spline P, and restores the
% q^2 centroid C to its saved initial C0. It never divides by peak curvature
% and does not impose a peak branch, basin, stage-time, or TV acceptance.
r=ops.rescaling;
peak=details.trackedPeakFunctional;
integral=ipm.evolve.wallGradientIntegral(flow.source(1,:),ops,r);
centroid=ipm.evolve.wallGradientCentroid(flow.source(1,:),ops);
if ~isfinite(peak.value) || peak.value<=0 || centroid.degenerateEnergy
    error('ipm:DegenerateWallIntegralGauge','The nonzero integral requires positive amplitude and centroid energy.');
end
context=struct('centroidNodes',centroid.indices,'centroidValue',centroid.value, ...
    'centroidDisplacement',centroid.value-centroid.referenceValue);
u1X=flow.u1*ipm.mesh.oddDx(ops)';
nominal=-interp1(ops.x,u1X(1,:),centroid.value,'linear');
[cL,cR,context]=ipm.evolve.generatorWallIntegralGauge( ...
    rho,flow,ops,peak,centroid.weights,integral,context,nominal);
cOmega=context.generatorAmplitudeRate;baseRhs=context.generatorBaseRhs;
v1=flow.u1+cL*ops.X+cR;v2=flow.u2+cL*ops.Y;
fullWallRhs=baseRhs(1,:)+cOmega*rho(1,:);
qRate=fullWallRhs*ops.Dx';
iRate=integral.weights*qRate';
pRate=peak.weights*qRate(peak.indices)';
cRate=centroid.weights*qRate(centroid.indices)';
cResidual=cRate+r.peakPositionRestoringRate*context.centroidDisplacement;
normalized=[iRate/integral.value,pRate/peak.value,cResidual];
if any(~isfinite(normalized)) || max(abs(normalized))>1e-5
    error('ipm:UnresolvedWallIntegralGauge', ...
        'Practical complete RHS I/P/C residuals exceed 1e-5: [%g %g %g].',normalized);
end
details.c_lNominal=nominal;details.lengthScaleGain=1;
details.widthCorrection=0;details.widthControlMode=0;details.widthRateScale=0;
details.lengthGaugeStrain=integral.value;details.lengthGaugeResidual=iRate;
details.omegaGaugeEnergy=peak.value;details.omegaGaugeForcing=-cOmega*peak.value;
details.omegaGaugeAbsoluteForcing=abs(details.omegaGaugeForcing);
details.omegaGaugeForcingCancellationRatio=1;
details.omegaGaugeResidual=pRate/max(1,abs(peak.value));
details.omegaPeakRateResidual=pRate;details.omegaGaugePeakX=peak.x;
details.omegaGaugePeakCondition=peak.value/max(max(abs(flow.source),[],'all'),realmin);
details.globalWallPeak=peak.value;details.globalWallPeakX=peak.x;
details.globalWallPeakTie=double(peak.tieNonunique);
details.wallCentroid=centroid.value;details.wallCentroidReference=centroid.referenceValue;
details.wallCentroidEnergy=centroid.energy;details.wallCentroidRate=cRate;
details.wallCentroidRateResidual=cResidual;
details.wallIntegral=integral.value;details.wallIntegralRadius=integral.radius;
details.wallIntegralSupportPoints=integral.supportPoints;details.wallIntegralCondition=integral.condition;
details.wallIntegralRate=iRate;
details.wallIntegralRelativeDrift=(integral.value-r.referenceWallIntegral)/max(abs(r.referenceWallIntegral),realmin);
details.peakHeightRate=pRate;
% The former peak-position columns are unused; explicit C fields carry phase.
details.peakPositionTargetRate=NaN;details.peakPositionRate=NaN;
details.pinGaugeResidual=cResidual;details.travelingWaveResidual=cResidual;
details.pinTransportVelocity=interp1(ops.x,v1(1,:),centroid.value,'linear');
details.peakTransportVelocity=interp1(ops.x,v1(1,:),peak.x,'linear');
details.transportAnchorX=centroid.value;
details.transportAnchorVelocity=details.pinTransportVelocity;
names={'generatorCflRate','generatorRootCondition','generatorReferenceVelocity', ...
    'generatorTranslationRate','generatorRootTolerance'};
for k=1:numel(names),details.(names{k})=context.(names{k});end
details.generatorRootResidual=max(abs(context.generatorRootResidual));
details.generatorIntegralResidual=normalized(1);details.generatorPeakResidual=normalized(2);
details.generatorPhaseResidual=normalized(3);
end
