function validateV2(result)
%IPM.OUTPUT.VALIDATEV2 Enforce the grouped result contract invariants.

if ~isstruct(result) || ~isscalar(result)
    error('ipm:ResultScalar','A solver result must be a scalar structure.');
end
if ~isfield(result,'schemaVersion') || ...
        ~isnumeric(result.schemaVersion) || ~isreal(result.schemaVersion) || ...
        ~isscalar(result.schemaVersion) || ~isfinite(result.schemaVersion) || ...
        result.schemaVersion ~= 2
    error('ipm:ResultSchemaVersion', ...
        ['Only results with explicit schemaVersion=2 are supported; ' ...
        'legacy, unversioned, and other result schemas are rejected.']);
end

required = {'metadata','state','grid','scale','physical','history', ...
    'quality','fit','config','elliptic','snapshots'};
require_fields(result,required,'version-2 result');
for index = 1:numel(required)
    name = required{index};
    if ~isstruct(result.(name)) || ~isscalar(result.(name))
        error('ipm:ResultGroup', ...
            'Version-2 result group "%s" must be a scalar structure.',name);
    end
end
validate_metadata(result.metadata);
validate_groups(result);
validate_config(result.config);

x = result.grid.x;
y = result.grid.y;
validateattributes(x,{'numeric'},{'vector','real','finite','increasing'});
validateattributes(y,{'numeric'},{'vector','real','finite','increasing'});
expectedSize = [numel(y),numel(x)];
validate_state(result.state,expectedSize);
validate_physical(result.physical,expectedSize);
validate_grid(result.grid,x,y);
validate_scale(result.scale);

assert_terminal_time(result.history.common.t, ...
    result.state.normalizedTime,'normalized');
assert_terminal_time(result.history.common.canonicalTau, ...
    result.state.canonicalTime,'canonical');
assert_terminal_time(result.history.common.physicalTime, ...
    result.state.physicalTime,'physical');
validate_current_canonical_contract(result);
if isfield(result.history,'wallCore')
    ipm.output.validateWallCoreHistory( ...
        result.history.wallCore,result.history.common);
end
validate_snapshots(result.snapshots,result.history.common);
validate_autonomous_mesh_controller(result);
end

function validate_autonomous_mesh_controller(result)
hasPolicy = isfield(result.config.remesh,'autonomousMesh');
if ~hasPolicy && ~isfield(result.metadata,'autonomousMesh')
    return;
end
current = struct();
if hasPolicy && result.config.remesh.autonomousMesh.enabled
    mesh = result.history.mesh;
    require_fields(mesh,{'coreGridPoints','verticalCoreGridPoints', ...
        'safetyFactor'},'autonomous result.history.mesh');
    if isempty(mesh.coreGridPoints) || isempty(mesh.verticalCoreGridPoints) || ...
            isempty(mesh.safetyFactor)
        error('ipm:AutonomousMeshController', ...
            'An enabled autonomous result requires terminal mesh observations.');
    end
    current = struct('step',result.state.steps, ...
        'canonicalTime',result.state.canonicalTime, ...
        'physicalTime',result.state.physicalTime, ...
        'normalizedTime',result.state.normalizedTime, ...
        'remeshCount',result.grid.remeshCount, ...
        'coreCells',[mesh.coreGridPoints(end),mesh.verticalCoreGridPoints(end)], ...
        'safety',mesh.safetyFactor(end),'x',result.grid.x,'y',result.grid.y, ...
        'history',result.history);
    if any(result.config.remesh.autonomousMesh.version == [2,3,4,5])
        require_fields(result.metadata,{'autonomousMesh'},'version-two metadata');
        require_fields(result.metadata.autonomousMesh,{'currentLevelId'},'version-two controller');
        current.levelId=result.metadata.autonomousMesh.currentLevelId;
        current.nodeCount=[numel(result.grid.x),numel(result.grid.y)];
        current.snapshots=struct('rho',{result.snapshots.rho}, ...
            'x',{result.snapshots.x},'y',{result.snapshots.y});
    end
end
ipm.remesh.validateController(result.metadata,result.config,current);
end

function validate_metadata(metadata)
require_fields(metadata,{'caseId','createdAt','solver'},'result.metadata');
validate_text(metadata.caseId,'result.metadata.caseId',false);
validate_text(metadata.createdAt,'result.metadata.createdAt',true);
validate_text(metadata.solver,'result.metadata.solver',false);
end

function validate_config(config)
schema = ipm.config.schema();
domains = schema.domainNames;
required = [{'schemaVersion','frozen'},domains];
if ~isstruct(config) || ~isscalar(config)
    error('ipm:ResultConfig', ...
        'result.config must be a grouped scalar configuration.');
end
missing = setdiff(required,fieldnames(config),'stable');
if ~isempty(missing)
    error('ipm:ResultConfig', ...
        'Grouped result.config is missing required field "%s".',missing{1});
end
unexpected = setdiff(fieldnames(config),required,'stable');
if ~isempty(unexpected)
    error('ipm:ResultConfig', ...
        'Unexpected top-level result.config field "%s".',unexpected{1});
end
if ~isnumeric(config.schemaVersion) || ~isreal(config.schemaVersion) || ...
        ~isscalar(config.schemaVersion) || ~isfinite(config.schemaVersion) || ...
        ~any(config.schemaVersion == [1,2,3,4]) || ...
        ~islogical(config.frozen) || ...
        ~isscalar(config.frozen) || ~config.frozen
    error('ipm:ResultConfig', ...
        ['result.config must be a frozen grouped schema-version-1, ' ...
        'schema-version-2, schema-version-3, or schema-version-4 ' ...
        'configuration.']);
end
for index = 1:numel(domains)
    name = domains{index};
    if ~isstruct(config.(name)) || ~isscalar(config.(name))
        error('ipm:ResultConfig', ...
            'result.config.%s must be a scalar structure.',name);
    end
    [declaredFields,requiredFields] = ...
        config_domain_fields(schema,name,config.schemaVersion);
    require_fields(config.(name),requiredFields, ...
        ['result.config.',name]);
    unexpected = setdiff(fieldnames(config.(name)), ...
        declaredFields,'stable');
    if ~isempty(unexpected)
        error('ipm:ResultConfig', ...
            'Unexpected result.config.%s field "%s".', ...
            name,unexpected{1});
    end
end
hasCustomX = isfield(config.grid,'customX');
hasCustomY = isfield(config.grid,'customY');
if hasCustomX ~= hasCustomY
    error('ipm:ResultConfig', ...
        ['Optional result.config.grid fields customX and customY must ' ...
        'occur together.']);
end
validate_config_values(config,schema);
end

function validate_config_values(config,schema)
% Check persisted normalized selectors without running defaults, derivation,
% or configuration resolution a second time during result finalization.
if config.schemaVersion <= 2
    legacyRateLimit = config.scaling.maxDynamicRate;
    if ~isnumeric(legacyRateLimit) || ~isreal(legacyRateLimit) || ...
            ~isscalar(legacyRateLimit) || ~isfinite(legacyRateLimit) || ...
            legacyRateLimit <= 0
        error('ipm:ResultConfig', ...
            ['Legacy result.config.scaling.maxDynamicRate must be a ' ...
            'positive finite scalar.']);
    end
    if ischar(config.scaling.cOmegaGauge) && ...
            strcmp(config.scaling.cOmegaGauge,'anchor_wall_window_l2')
        error('ipm:ResultConfig', ...
            ['The anchor_wall_window_l2 gauge was introduced with ' ...
            'configuration schema version 3.']);
    end
end
if config.schemaVersion <= 3
    validate_legacy_feedback_fields(config.scaling);
    futureGauges = {'anchor_wall_template_projection', ...
        'anchor_bulk_gradient_l2','anchor_wall_window_l4', ...
        'wall_omega_quadratic_peak','wall_omega_tracked_spline_peak', ...
        'wall_omega_global_spline_maximum', ...
        'outer_wall_density_window_l2'};
    if any(strcmp(config.scaling.cOmegaGauge,futureGauges))
        error('ipm:ResultConfig', ...
            ['The selected c_omega gauge was introduced with ' ...
            'configuration schema version 4.']);
    end
    if strcmp(config.scaling.travelingWaveGauge,'wall_gradient_centroid')
        error('ipm:ResultConfig', ...
            'The wall-gradient centroid position was introduced with configuration schema version 4.');
    end
else
    validate_scaling_contract(config.scaling);
end
if config.schemaVersion == 4 && ...
        any(strcmp(config.scaling.scalingContract, ...
        {'wall_integral_generators_v2','wall_integral_centered_generators_v3'})) && ...
        (~strcmp(config.transport.spatialDiscretization,'high_order') || ...
        ~strcmp(config.physics.symmetryMode,'half_plane') || ...
        ~strcmp(config.scaling.rescalingMode,'dynamic') || ...
        ~strcmp(config.scaling.dynamicScaleGeometry,'isotropic') || ...
        ~strcmp(config.scaling.lengthGauge,'wall_gradient_integral') || ...
        ~strcmp(config.transport.transportScheme,'weno5_fd') || ...
        ~strcmp(config.transport.transportBoundaryMode,'open') || ...
        ~strcmp(config.scaling.travelingWaveGauge,'gradient_peak') || ...
        ~any(strcmp(config.scaling.cOmegaGauge, ...
        {'wall_omega_quadratic_peak','wall_omega_tracked_spline_peak'})))
    error('ipm:ResultConfig', ...
        ['The fixed-generator wall-integral contract requires high_order, half_plane ' ...
        'isotropic dynamic scaling, open WENO5-FD, wall_gradient_integral, ' ...
        'gradient_peak, and a declared quadratic or tracked-spline peak.']);
end
if config.schemaVersion==4 && ...
        strcmp(config.scaling.scalingContract,'wall_integral_centroid_practical_v1') && ...
        (~strcmp(config.transport.spatialDiscretization,'high_order') || ...
        ~strcmp(config.physics.symmetryMode,'half_plane') || ...
        ~strcmp(config.scaling.rescalingMode,'dynamic') || ...
        ~strcmp(config.scaling.dynamicScaleGeometry,'isotropic') || ...
        ~strcmp(config.scaling.lengthGauge,'wall_gradient_integral') || ...
        ~strcmp(config.transport.transportScheme,'weno5_fd') || ...
        ~strcmp(config.transport.transportBoundaryMode,'open') || ...
        ~strcmp(config.scaling.travelingWaveGauge,'wall_gradient_centroid') || ...
        ~strcmp(config.scaling.cOmegaGauge,'wall_omega_global_spline_maximum'))
    error('ipm:ResultConfig','Invalid practical integral/global-maximum/centroid tuple.');
end
if config.schemaVersion == 4 && ...
        (strcmp(config.scaling.cOmegaGauge,'wall_omega_global_spline_maximum') || ...
        strcmp(config.scaling.travelingWaveGauge,'wall_gradient_centroid')) && ...
        ~strcmp(config.scaling.scalingContract,'wall_integral_centroid_practical_v1')
    error('ipm:ResultConfig', ...
        'Global spline amplitude and wall centroid position require the practical contract.');
end
optionNames = fieldnames(schema.options);
for index = 1:numel(optionNames)
    name = optionNames{index};
    descriptor = schema.options.(name);
    domain = descriptor.domain;
    if ~isfield(config.(domain),name)
        continue;
    end
    value = config.(domain).(name);
    if strcmp(descriptor.kind,'enum')
        if ~ischar(value) || ~isrow(value) || ...
                ~any(strcmp(value,descriptor.allowed))
            error('ipm:ResultConfig', ...
                'result.config.%s.%s has an invalid selector value.', ...
                domain,name);
        end
    elseif strcmp(descriptor.kind,'boolean') && ...
            (~islogical(value) || ~isscalar(value))
        error('ipm:ResultConfig', ...
            'result.config.%s.%s must be a normalized logical scalar.', ...
            domain,name);
    elseif strcmp(descriptor.kind,'checkpoint')
        validate_checkpoint_policy(value);
    elseif strcmp(descriptor.kind,'autonomous_mesh')
        validate_autonomous_mesh_policy(value,config,schema);
    elseif strcmp(descriptor.kind,'initial_mesh_observation')
        choices = struct();
        for d = 1:numel(schema.domainNames)
            group = config.(schema.domainNames{d});names = fieldnames(group);
            for j = 1:numel(names),choices.(names{j}) = group.(names{j});end
        end
        normalized = ipm.config.initialMeshObservationPolicy(value,choices);
        if ~isequaln(normalized,value) || ~same_policy_types(normalized,value)
            error('ipm:ResultConfig','The initial observation policy must be fully canonical.');
        end
    end
end

if config.schemaVersion == 4 && strcmp(config.scaling.cOmegaGauge, ...
        'anchor_wall_template_projection') && config.remesh.adaptiveRemesh
    error('ipm:ResultConfig', ...
        ['The schema-version-4 anchor-wall template projection ' ...
        'gauge requires adaptiveRemesh=false.']);
end
if config.schemaVersion == 4 && ...
        strcmp(config.scaling.cOmegaGauge,'anchor_wall_strain')
    error('ipm:ResultConfig', ...
        ['The approximate historical anchor_wall_strain selector is not ' ...
        'available under the exact schema-version-4 contract.']);
end

transportScheme = config.transport.transportScheme;
if config.schemaVersion == 1
    % Schema 1 can record MUSCL or the former nonuniform WENO choice, but
    % its absent fields imply legacy spatial operators, SSPRK3, and PCHIP.
    spatialDiscretization = 'legacy_second_order';
    timeIntegrator = 'ssprk3';
    remeshTransferScheme = 'pchip';
else
    spatialDiscretization = config.transport.spatialDiscretization;
    timeIntegrator = config.time.timeIntegrator;
    remeshTransferScheme = config.remesh.remeshTransferScheme;
end

[tupleIdentifier,tupleMessage] = ipm.config.tupleViolation( ...
    transportScheme,spatialDiscretization,timeIntegrator, ...
    remeshTransferScheme, ...
    config.scaling.dynamicScaleGeometry, ...
    config.elliptic.anisotropicPoissonSolver);
if ~isempty(tupleIdentifier)
    error('ipm:ResultConfig', ...
        'Persisted numerical choices are invalid: %s',tupleMessage);
end
end

function validate_legacy_feedback_fields(scaling)
if ~islogical(scaling.adaptiveLengthScaling) || ...
        ~isscalar(scaling.adaptiveLengthScaling)
    error('ipm:ResultConfig', ...
        ['Legacy result.config.scaling.adaptiveLengthScaling must be a ' ...
        'normalized logical scalar.']);
end
validate_legacy_feedback_scalar( ...
    scaling.widthExpansionStrength,0,false,'widthExpansionStrength');
validate_legacy_feedback_scalar( ...
    scaling.widthContractionOnset,1,true,'widthContractionOnset');
validate_legacy_feedback_scalar( ...
    scaling.widthContractionStrength,0,true,'widthContractionStrength');
validate_legacy_feedback_scalar( ...
    scaling.maxWidthRateCorrection,0,true,'maxWidthRateCorrection');
end

function validate_legacy_feedback_scalar(value,bound,strict,name)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value);
if strict
    valid = valid && value > bound;
else
    valid = valid && value >= bound;
end
if ~valid
    error('ipm:ResultConfig', ...
        'Legacy result.config.scaling.%s is invalid.',name);
end
end

function validate_scaling_contract(scaling)
if ~ischar(scaling.scalingContract) || ...
        ~isrow(scaling.scalingContract) || ...
        ~any(strcmp(scaling.scalingContract, ...
        {'exact_gauge_no_feedback_v1','peak_curvature_restoring_v1', ...
        'wall_integral_restoring_v1','wall_integral_generators_v2', ...
        'wall_integral_centered_generators_v3','wall_integral_centroid_practical_v1'}))
    error('ipm:ResultConfig', ...
        'Schema-version-4 scalingContract is invalid.');
end
if strcmp(scaling.scalingContract,'peak_curvature_restoring_v1') ~= ...
        strcmp(scaling.lengthGauge,'wall_gradient_curvature')
    error('ipm:ResultConfig', ...
        'The peak-curvature restoring contract and length gauge must agree.');
end
if strcmp(scaling.scalingContract,'peak_curvature_restoring_v1') && ...
        (~isfield(scaling,'peakPositionRestoringRate') || ...
        ~isfield(scaling,'peakCurvatureRadius'))
    error('ipm:ResultConfig', ...
        'Peak-curvature results must persist both gauge parameters.');
end
wallIntegralContract = any(strcmp(scaling.scalingContract, ...
    {'wall_integral_restoring_v1','wall_integral_generators_v2', ...
    'wall_integral_centered_generators_v3','wall_integral_centroid_practical_v1'}));
if wallIntegralContract ~= ...
        strcmp(scaling.lengthGauge,'wall_gradient_integral')
    error('ipm:ResultConfig', ...
        'The wall-integral contract and length gauge must agree.');
end
if wallIntegralContract && ...
        (~isfield(scaling,'peakPositionRestoringRate') || ...
        ~isfield(scaling,'wallIntegralRadius'))
    error('ipm:ResultConfig', ...
        'Wall-integral results must persist both gauge parameters.');
end
names = {'lengthScaleGain','widthGaugeGain','omegaGaugeGain', ...
    'travelingWaveGain'};
expected = [1,0,0,0];
for index = 1:numel(names)
    value = scaling.(names{index});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || value ~= expected(index)
        error('ipm:ResultConfig', ...
            ['Schema-version-4 result.config.scaling.%s must equal ' ...
            'the no-feedback sentinel %.0f.'],names{index},expected(index));
    end
end
end

function validate_checkpoint_policy(policy)
required = {'enabled','file','every','atExit'};
if ~isstruct(policy) || ~isscalar(policy) || ...
        ~isempty(setxor(fieldnames(policy),required)) || ...
        ~islogical(policy.enabled) || ~isscalar(policy.enabled) || ...
        ~islogical(policy.atExit) || ~isscalar(policy.atExit) || ...
        ~(ischar(policy.file) && isrow(policy.file) && ~isempty(policy.file)) || ...
        ~isnumeric(policy.every) || ~isreal(policy.every) || ...
        ~isscalar(policy.every) || isnan(policy.every) || policy.every <= 0
    error('ipm:ResultConfig', ...
        'result.config.output.checkpoint has an invalid normalized policy.');
end
end

function validate_autonomous_mesh_policy(value,config,schema)
% Reuse the policy contract against only the choices actually persisted.
% Normalization is compared, never installed into the historical config.
choices = struct();
for index = 1:numel(schema.domainNames)
    group = config.(schema.domainNames{index});
    names = fieldnames(group);
    for fieldIndex = 1:numel(names)
        name = names{fieldIndex};
        choices.(name) = group.(name);
    end
end
try
    normalized = ipm.config.autonomousMeshPolicy(value,choices);
catch exception
    error('ipm:ResultConfig', ...
        'Persisted autonomousMesh policy is invalid: %s',exception.message);
end
if ~isequaln(normalized,value) || ~same_policy_types(normalized,value)
    error('ipm:ResultConfig', ...
        ['result.config.remesh.autonomousMesh must contain the complete ' ...
        'normalized policy exactly as persisted.']);
end
end

function same = same_policy_types(normalized,value)
% isequaln permits logical/numeric and char/string equality; persisted
% normalized policies must also retain each field's canonical type.
same = strcmp(class(normalized),class(value)) && ...
    isequal(size(normalized),size(value));
if ~same || ~isstruct(normalized)
    return;
end
names = fieldnames(normalized);
for index = 1:numel(names)
    name = names{index};
    if ~same_policy_types(normalized.(name),value.(name))
        same = false;
        return;
    end
end
end

function [declaredFields,requiredFields] = ...
        config_domain_fields(schema,domain,configVersion)
optionNames = fieldnames(schema.options);
inDomain = false(size(optionNames));
hasDefault = false(size(optionNames));
for index = 1:numel(optionNames)
    descriptor = schema.options.(optionNames{index});
    inDomain(index) = strcmp(descriptor.domain,domain);
    hasDefault(index) = descriptor.hasDefault;
end
if configVersion <= 2
    % Configuration schema 3 introduced this window-gauge control. Do not
    % require or accept it when validating an unchanged historical config.
    addedInVersionThree = {'omegaGaugeWindowRadius'};
    inDomain(ismember(optionNames,addedInVersionThree)) = false;
end
if configVersion <= 3
    % Configuration schema 4 made the exact no-feedback contract explicit
    % and retired every state-dependent length-rate feedback control.
    addedInVersionFour = {'scalingContract','autonomousMesh', ...
        'initialMeshObservationFallback', ...
        'peakPositionRestoringRate','peakCurvatureRadius','wallIntegralRadius', ...
        'wenoAlphaMode'};
    inDomain(ismember(optionNames,addedInVersionFour)) = false;
end
% This experimental selector is optional and belongs to schema-4 only.
if configVersion <= 3
    inDomain(ismember(optionNames,{'wenoInflowDensityDegree'})) = false;
end
% Historical schema-4 results predate this optional restoring selector.
% New runs persist it, while old files remain readable.
hasDefault(ismember(optionNames,{'peakPositionRestoringRate'})) = false;
hasDefault(ismember(optionNames,{'peakCurvatureRadius'})) = false;
hasDefault(ismember(optionNames,{'wallIntegralRadius'})) = false;
hasDefault(ismember(optionNames,{'wenoAlphaMode'})) = false;
hasDefault(ismember(optionNames,{'wenoInflowDensityDegree'})) = false;
if configVersion == 1
    % Configuration schema 1 predates the explicitly selectable high-order
    % numerical path. Keep existing version-1 result files readable without
    % pretending that their absent choices were persisted.
    addedInVersionTwo = {'timeIntegrator','spatialDiscretization', ...
        'remeshTransferScheme'};
    inDomain(ismember(optionNames,addedInVersionTwo)) = false;
end
declaredFields = optionNames(inDomain);
requiredFields = optionNames(inDomain & hasDefault);
if configVersion <= 3 && strcmp(domain,'scaling')
    retiredInVersionFour = {'adaptiveLengthScaling'; ...
        'widthExpansionStrength';'widthContractionOnset'; ...
        'widthContractionStrength';'maxWidthRateCorrection'};
    declaredFields = [declaredFields;retiredInVersionFour];
    requiredFields = [requiredFields;retiredInVersionFour];
end
if configVersion <= 2 && strcmp(domain,'scaling')
    % Config schemas 1 and 2 persisted the former time-reparameterization
    % limiter. It remains part of their immutable historical record even
    % though later schemas removed it from the evolution equations.
    declaredFields = [declaredFields;{'maxDynamicRate'}];
    requiredFields = [requiredFields;{'maxDynamicRate'}];
end
end

function validate_groups(result)
require_fields(result.state, ...
    {'rho','psi','omega','velocity','normalizedTime','canonicalTime', ...
    'physicalTime','steps','stopReason'},'result.state');
if ~isstruct(result.state.velocity) || ~isscalar(result.state.velocity)
    error('ipm:ResultVelocity','result.state.velocity must be scalar.');
end
require_fields(result.state.velocity,{'x','y'},'result.state.velocity');
require_fields(result.grid, ...
    {'x','y','physicalX','physicalY','comovingPhysicalX', ...
    'trackedPeakX','physicalTrackedPeakX','remeshCount'},'result.grid');
require_fields(result.scale, ...
    {'geometry','Cx','Cy','Comega','Lx','Ly','aspect','Xshift', ...
    'rates','canonicalRates','conservativeSource', ...
    'canonicalConservativeSource'},'result.scale');
if ~isstruct(result.scale.rates) || ~isscalar(result.scale.rates) || ...
        ~isstruct(result.scale.canonicalRates) || ...
        ~isscalar(result.scale.canonicalRates)
    error('ipm:ResultScaleRates', ...
        'Scale rates must be scalar structures.');
end
require_fields(result.scale.rates,{'cx','cy','comega','cr'}, ...
    'result.scale.rates');
require_fields(result.scale.canonicalRates,{'cx','cy','comega','cr'}, ...
    'result.scale.canonicalRates');
require_fields(result.physical, ...
    {'psi','rho','omega','rhoX','rhoY','u1','u2'},'result.physical');
validate_history(result.history);
require_fields(result.elliptic,{'solveInfo'},'result.elliptic');
require_fields(result.snapshots, ...
    {'rho','normalizedTime','physicalTime','canonicalTime','x','y'}, ...
    'result.snapshots');
end

function validate_history(history)
groups = {'common','gauge','mesh','anisotropic'};
require_fields(history,groups,'result.history');
for index = 1:numel(groups)
    name = groups{index};
    if ~isstruct(history.(name)) || ~isscalar(history.(name))
        error('ipm:ResultHistoryGroup', ...
            'result.history.%s must be a scalar structure.',name);
    end
end
if isfield(history.anisotropic,'mode')
    error('ipm:ResultHistoryField', ...
        ['The early-v2 result.history.anisotropic.mode field is not ' ...
        'supported; use anisotropicGaugeMode in schemaVersion=2 results.']);
end
reject_duplicate_history_fields(history,groups);
require_fields(history.common,{'t','canonicalTau','physicalTime'}, ...
    'result.history.common');
commonCount = aligned_history_length(history.common,'common');
if commonCount < 1
    error('ipm:ResultHistoryEmpty', ...
        'result.history.common must contain at least one recorded state.');
end
for index = 2:numel(groups)
    name = groups{index};
    count = aligned_history_length(history.(name),name);
    if count ~= 0 && count ~= commonCount
        error('ipm:ResultHistoryLength', ...
            ['Nonempty result.history.%s must have the same record count ' ...
            'as result.history.common.'],name);
    end
end
end

function reject_duplicate_history_fields(history,groups)
seen = {};
for groupIndex = 1:numel(groups)
    names = fieldnames(history.(groups{groupIndex}));
    duplicateIndex = find(ismember(names,seen),1);
    if ~isempty(duplicateIndex)
        name = names{duplicateIndex};
        error('ipm:HistoryDuplicateField', ...
            'History field "%s" occurs in more than one group.',name);
    end
    seen = [seen;names]; %#ok<AGROW>
end
end

function validate_state(state,expectedSize)
fieldNames = {'rho','psi','omega'};
for index = 1:numel(fieldNames)
    validate_field(state.(fieldNames{index}),expectedSize, ...
        ['result.state.',fieldNames{index}]);
end
for name = {'x','y'}
    validate_field(state.velocity.(name{1}),expectedSize, ...
        ['result.state.velocity.',name{1}]);
end
for name = {'normalizedTime','canonicalTime','physicalTime'}
    validate_scalar(state.(name{1}),['result.state.',name{1}]);
end
validateattributes(state.steps,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'},mfilename, ...
    'result.state.steps');
validate_text(state.stopReason,'result.state.stopReason',false);
end

function validate_physical(physical,expectedSize)
for name = {'psi','rho','omega','rhoX','rhoY','u1','u2'}
    validate_field(physical.(name{1}),expectedSize, ...
        ['result.physical.',name{1}]);
end
end

function validate_field(value,expectedSize,name)
if ~isnumeric(value) || ~isreal(value) || ...
        any(~isfinite(value),'all') || ~isequal(size(value),expectedSize)
    error('ipm:ResultField', ...
        '%s must be finite, real, numeric, and match the state grid.',name);
end
end

function validate_grid(grid,x,y)
validate_coordinate(grid.physicalX,numel(x),'result.grid.physicalX');
validate_coordinate(grid.comovingPhysicalX,numel(x), ...
    'result.grid.comovingPhysicalX');
validate_coordinate(grid.physicalY,numel(y),'result.grid.physicalY');
validate_scalar(grid.trackedPeakX,'result.grid.trackedPeakX');
validate_scalar(grid.physicalTrackedPeakX, ...
    'result.grid.physicalTrackedPeakX');
validateattributes(grid.remeshCount,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'},mfilename, ...
    'result.grid.remeshCount');
end

function validate_scale(scale)
validate_text(scale.geometry,'result.scale.geometry',false);
geometry = char(string(scale.geometry));
if ~any(strcmp(geometry,{'isotropic','anisotropic'}))
    error('ipm:ResultScaleGeometry', ...
        'result.scale.geometry must be isotropic or anisotropic.');
end
for name = {'Cx','Cy','Comega','Lx','Ly','aspect'}
    validateattributes(scale.(name{1}),{'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename, ...
        ['result.scale.',name{1}]);
end
validate_scalar(scale.Xshift,'result.scale.Xshift');
for name = {'cx','cy','comega','cr'}
    validate_scalar(scale.rates.(name{1}), ...
        ['result.scale.rates.',name{1}]);
    validate_scalar(scale.canonicalRates.(name{1}), ...
        ['result.scale.canonicalRates.',name{1}]);
end
validate_scalar(scale.conservativeSource, ...
    'result.scale.conservativeSource');
validate_scalar(scale.canonicalConservativeSource, ...
    'result.scale.canonicalConservativeSource');
assert_close(scale.Lx,1/scale.Cx,'Lx and Cx');
assert_close(scale.Ly,1/scale.Cy,'Ly and Cy');
assert_close(scale.aspect,scale.Cy/scale.Cx,'aspect and Cy/Cx');
assert_close(scale.conservativeSource, ...
    scale.rates.cx+scale.rates.cy+scale.rates.comega, ...
    'conservative source and rates');
assert_close(scale.canonicalConservativeSource, ...
    scale.canonicalRates.cx+scale.canonicalRates.cy+ ...
    scale.canonicalRates.comega,'canonical conservative source and rates');
end

function validate_current_canonical_contract(result)
% Configuration schemas 3 and 4 use the direct canonical clock. Keep
% that removal machine-checkable without retroactively rejecting preserved
% schema-1/2 results produced by the historical evolution contract.
if result.config.schemaVersion < 3
    return;
end

common = result.history.common;
directNames = {'c_l','c_x','c_y','c_omega','c_r', ...
    'conservativeSource'};
canonicalNames = {'canonicalCL','canonicalCX','canonicalCY', ...
    'canonicalCOmega','canonicalCR','canonicalConservativeSource'};
require_fields(common,[{'timeSpeed'},directNames,canonicalNames], ...
    'direct-clock result.history.common');

assert_contract_close(common.t,common.canonicalTau, ...
    'normalized and canonical history clocks', ...
    'ipm:ResultCanonicalClock');
assert_contract_close(result.state.normalizedTime, ...
    result.state.canonicalTime,'normalized and canonical terminal clocks', ...
    'ipm:ResultCanonicalClock');
assert_contract_close(common.timeSpeed,ones(size(common.timeSpeed)), ...
    'unit canonical time speed','ipm:ResultCanonicalClock');

for index = 1:numel(directNames)
    assert_contract_close(common.(directNames{index}), ...
        common.(canonicalNames{index}), ...
        [directNames{index},' and ',canonicalNames{index}], ...
        'ipm:ResultCanonicalRate');
end
terminalDirect = {'cx','cy','comega','cr'};
for index = 1:numel(terminalDirect)
    name = terminalDirect{index};
    assert_contract_close(result.scale.rates.(name), ...
        result.scale.canonicalRates.(name), ...
        ['terminal ',name,' rate and canonical rate'], ...
        'ipm:ResultCanonicalRate');
end
assert_contract_close(result.scale.conservativeSource, ...
    result.scale.canonicalConservativeSource, ...
    'terminal conservative source and canonical source', ...
    'ipm:ResultCanonicalRate');
end

function validate_snapshots(snapshots,history)
count = numel(snapshots.rho);
for name = {'rho','x','y'}
    if ~iscell(snapshots.(name{1})) || ...
            numel(snapshots.(name{1})) ~= count
        error('ipm:ResultSnapshotContract', ...
            'result.snapshots.%s must be a cell array of snapshot count.', ...
            name{1});
    end
end
for name = {'normalizedTime','physicalTime','canonicalTime'}
    value = snapshots.(name{1});
    validEmpty = count == 0 && isnumeric(value) && isempty(value);
    validSeries = count > 0 && isnumeric(value) && isreal(value) && ...
        isvector(value) && numel(value) == count && ...
        all(isfinite(value),'all');
    if ~(validEmpty || validSeries)
        error('ipm:ResultSnapshotTime', ...
            'result.snapshots.%s must contain one finite time per snapshot.', ...
            name{1});
    end
end
historyCount = numel(history.t);
if count ~= 0 && count ~= historyCount
    error('ipm:ResultSnapshotHistoryCount', ...
        'Stored snapshots must be paired with every recorded history row.');
end
for index = 1:count
    snapshotX = snapshots.x{index};
    snapshotY = snapshots.y{index};
    validateattributes(snapshotX,{'numeric'}, ...
        {'vector','real','finite','increasing'});
    validateattributes(snapshotY,{'numeric'}, ...
        {'vector','real','finite','increasing'});
    snapshotRho = snapshots.rho{index};
    if ~isnumeric(snapshotRho) || ~isreal(snapshotRho) || ...
            any(~isfinite(snapshotRho),'all') || ...
            ~isequal(size(snapshotRho),[numel(snapshotY),numel(snapshotX)])
        error('ipm:ResultSnapshotSize', ...
            'Snapshot %d must be finite and match its x/y grid.',index);
    end
end
if count > 0
    assert_time_series(snapshots.normalizedTime,history.t,'normalized');
    assert_time_series(snapshots.canonicalTime, ...
        history.canonicalTau,'canonical');
    assert_time_series(snapshots.physicalTime, ...
        history.physicalTime,'physical');
end
end

function count = aligned_history_length(group,name)
fields = fieldnames(group);
if isempty(fields)
    count = 0;
    return;
end
count = numel(group.(fields{1}));
for index = 1:numel(fields)
    value = group.(fields{index});
    if ~isvector(value) || numel(value) ~= count
        error('ipm:ResultHistoryLength', ...
            'Fields in result.history.%s must be aligned vectors.',name);
    end
end
end

function validate_coordinate(value,count,name)
validateattributes(value,{'numeric'}, ...
    {'vector','real','finite','increasing','numel',count},mfilename,name);
end

function validate_scalar(value,name)
validateattributes(value,{'numeric'}, ...
    {'scalar','real','finite'},mfilename,name);
end

function validate_text(value,name,allowEmpty)
valid = (ischar(value) && (isrow(value) || (allowEmpty && isempty(value)))) || ...
    (isstring(value) && isscalar(value) && ~ismissing(value));
if valid && ~allowEmpty
    valid = strlength(string(value)) > 0;
end
if ~valid
    error('ipm:ResultText','%s must be scalar text.',name);
end
end

function assert_terminal_time(values,terminalValue,name)
if ~isnumeric(values) || ~isreal(values) || ...
        any(~isfinite(values)) || isempty(values)
    error('ipm:ResultHistoryTime', ...
        'The %s history must be a nonempty finite numeric vector.',name);
end
assert_close(values(end),terminalValue,[name,' terminal time']);
end

function assert_time_series(actual,expected,name)
actual = actual(:);
expected = expected(:);
if numel(actual) ~= numel(expected)
    error('ipm:ResultSnapshotTimeAlignment', ...
        'Snapshot and history %s times have different lengths.',name);
end
scale = max([ones(size(actual)),abs(actual),abs(expected)],[],2);
if any(abs(actual-expected) > 100*eps(scale))
    error('ipm:ResultSnapshotTimeAlignment', ...
        'Snapshot and history %s times are not paired.',name);
end
end

function assert_contract_close(actual,expected,name,identifier)
if ~isnumeric(actual) || ~isnumeric(expected) || ...
        ~isreal(actual) || ~isreal(expected) || ...
        ~isvector(actual) || ~isvector(expected) || ...
        numel(actual) ~= numel(expected) || isempty(actual)
    error(identifier, ...
        'Direct-clock fields for %s must be aligned numeric vectors.', ...
        name);
end
actual = actual(:);
expected = expected(:);
if any(~isfinite(actual)) || any(~isfinite(expected))
    error(identifier, ...
        'Direct-clock fields for %s must be finite.',name);
end
scale = max([ones(size(actual)),abs(actual),abs(expected)],[],2);
if any(abs(actual-expected) > 100*eps(scale))
    error(identifier, ...
        'Direct-clock fields for %s violate the canonical identity.', ...
        name);
end
end

function assert_close(actual,expected,name)
scale = max([1,abs(actual),abs(expected)]);
if abs(actual-expected) > 100*eps(scale)
    error('ipm:ResultConsistency', ...
        'Result fields for %s are inconsistent.',name);
end
end

function require_fields(value,names,context)
for index = 1:numel(names)
    if ~isfield(value,names{index})
        error('ipm:ResultMissingField', ...
            '%s is missing required field "%s".',context,names{index});
    end
end
end
