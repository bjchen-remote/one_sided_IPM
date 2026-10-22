function metadata = initializeMetadata(config,ops)
%IPM.OUTPUT.INITIALIZEMETADATA Create one identity and reserve output names.

scaling = config.scaling;
output = config.output;

timestamp = datetime('now','TimeZone','local');
createdTimestamp = timestamp;
createdTimestamp.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX';
caseTimestamp = timestamp;
caseTimestamp.Format = 'yyyyMMdd''T''HHmmssSSS';
[~,token] = fileparts(tempname);
caseId = sprintf('%s_%s_%s_%dx%d_%s',char(caseTimestamp), ...
    scaling.rescalingMode,ops.dynamicScaleGeometry,ops.nx,ops.ny,token);

resultFile = '';
if output.saveResults
    resultFile = ipm.output.uniquePath(output.resultFile,caseId);
end
videoFile = '';
if output.writeVideo
    videoFile = ipm.output.uniquePath(output.videoFile,caseId);
end
metadata = struct('caseId',caseId, ...
    'createdAt',string(createdTimestamp),'solver','IPM', ...
    'caseMetadata',output.caseMetadata, ...
    'gaugeContract',ipm.evolve.gaugeContract(config,ops), ...
    'resultFile',string(resultFile),'videoFile',string(videoFile));
if isfield(ops,'timeIntegrationContract')
    metadata.timeIntegrationContract=ops.timeIntegrationContract;
    metadata.transportCflPolicy='complete_rate_recorded; implicit_stage_convergence_and_time_accuracy_control';
    metadata.implicitTimeMethodExperimental=true;
end
if strcmp(scaling.cOmegaGauge,'wall_omega_tracked_spline_peak')
    metadata.peakReconstruction='tracked_natural_cubic_spline';
    metadata.peakObservablePreserved=false;
    metadata.productionSpatialOperatorPreserved= ...
        ~any(strcmp(scaling.scalingContract, ...
        {'wall_integral_generators_v2','wall_integral_centered_generators_v3'}));
    metadata.legacyQuadraticFieldsAreAliases=true;
end
if strcmp(scaling.scalingContract,'wall_integral_generators_v2')
    metadata.transportDiscretization='fixed_signed_generator_weno5_fd';
    metadata.productionSpatialOperatorPreserved=false;
end
if strcmp(scaling.scalingContract,'wall_integral_centered_generators_v3')
    metadata.transportDiscretization='centered_signed_generator_weno5_fd';
    metadata.productionSpatialOperatorPreserved=false;
    metadata.generatorReferenceVelocity='Uref=interp1(x,U1(1,:),a,''linear'')';
    metadata.generatorTranslationRate='delta=c_r+Uref';
    metadata.generatorCflMask='forbidden';
end
if strcmp(scaling.scalingContract,'wall_integral_centroid_practical_v1')
    metadata.peakReconstruction='global_full_domain_C2_natural_cubic_spline_maximum';
    metadata.peakObservablePreserved=false;
    metadata.productionSpatialOperatorPreserved=false;
    metadata.legacyQuadraticFieldsAreAliases=false;
    metadata.legacyQuadraticFields='unused_in_practical_mode';
    metadata.legacyPeakPositionColumns='unused_and_not_recorded';
    metadata.maximumDerivative='selected_maximizer_evaluation_at_ties_not_unique_maximum_derivative';
    metadata.positionObservable='positive_trapezoidal_wall_R_X_squared_centroid_on_fixed_[-L,L]';
    metadata.centroidReference='C0=initial_centroid_persisted_in_runtime_references_and_checkpoint';
    metadata.exactPeakPinning=false;
    metadata.transportDiscretization='centered_signed_generator_weno5_fd';
    metadata.generatorReferenceVelocity='Uref=interp1(x,U1(1,:),C,''linear'')';
    metadata.generatorTranslationRate='delta=c_r+Uref';
    metadata.generatorCflMask='forbidden';
    metadata.gaugeRootGate=1e-6;
    metadata.gaugeFullAuditGate=1e-5;
    metadata.qualityAcceptance='finite_state_full_generator_CFL_valid_mesh';
    metadata.diagnosticOnlyQualityMetrics={'wall_TV','wall_peak_count','density_range', ...
        'positiveWallNegativeRatio','feature_point_count'};
    metadata.gaugeTransportAssembly=scaling.gaugeTransportAssembly;
    if ismember(scaling.gaugeTransportAssembly,{'net_velocity','net_velocity_density_centered', ...
            'net_velocity_density_centered_advective'})
        metadata.transportDiscretization='centered_net_velocity_weno5_fd';
        if ismember(scaling.gaugeTransportAssembly, ...
                {'net_velocity_density_centered','net_velocity_density_centered_advective'})
            metadata.transportDiscretization='density_centered_net_velocity_weno5_fd';
        end
        if strcmp(scaling.gaugeTransportAssembly,'net_velocity_density_centered_advective')
            metadata.transportDiscretization='density_centered_net_velocity_advective_weno5_fd';
            metadata.transportGradient='one_sided_density_WENO5_with_unit_dummy_velocities_then_physical_J_division';
            metadata.transportPhysicalInflow='same_actual_net_velocity_endpoint_masks_for_both_density_derivatives';
        end
        metadata.gaugeSolve='finite_difference_Newton_complete_wall_then_full_2D_RHS_audit';
        metadata.legacyGeneratorCflRateColumn='actual_net_J_dXi_transport_rate';
        metadata.qualityAcceptance='finite_state_full_net_CFL_valid_mesh';
    end
end
end
