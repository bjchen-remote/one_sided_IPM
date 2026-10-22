function config = resolve(userOpts)
%IPM.CONFIG.RESOLVE Resolve flat overrides into one grouped configuration.

if nargin < 1
    userOpts = struct();
end
if ~isstruct(userOpts) || ~isscalar(userOpts)
    error('ipm:BadOptions','userOpts must be a scalar structure.');
end

schema = ipm.config.schema();
optionNames = fieldnames(schema.options);
values = default_values(schema,optionNames);

overrideNames = fieldnames(userOpts);
unknownNames = setdiff(overrideNames,optionNames);
if ~isempty(unknownNames)
    error('ipm:UnknownOption','Unknown IPM option(s): %s.',strjoin(unknownNames,', '));
end
for index = 1:numel(overrideNames)
    name = overrideNames{index};
    values.(name) = userOpts.(name);
end

% Normalize every derivation control before any dependent branch uses it.
for index = 1:numel(optionNames)
    name = optionNames{index};
    descriptor = schema.options.(name);
    if isfield(values,name) && strcmp(descriptor.kind,'boolean')
        values.(name) = normalize_boolean(values.(name),name,descriptor.errorId);
    end
end
for index = 1:numel(optionNames)
    name = optionNames{index};
    descriptor = schema.options.(name);
    if isfield(values,name) && strcmp(descriptor.kind,'enum')
        values.(name) = normalize_enum(values.(name),name,descriptor.allowed,descriptor.errorId);
    end
end
initialDescriptor = schema.options.initialCondition;
values.initialCondition = normalize_initial_condition(values.initialCondition,initialDescriptor);
values = apply_custom_domain(values,schema);

validate_derivation_inputs(values,schema);

if strcmp(values.dynamicScaleGeometry,'anisotropic') && ...
        ~isfield(userOpts,'adaptiveRemesh')
    values.adaptiveRemesh = false;
end
if ~isfield(userOpts,'lengthGauge') && ...
        strcmp(values.symmetryMode,'half_plane')
    values.lengthGauge = 'local_strain';
end
if isfield(userOpts,'physicalFinalTime') && ...
        ~isfield(userOpts,'finalTime') && isfinite(values.physicalFinalTime)
    values.finalTime = Inf;
end
defaultPaths = ipm.config.outputPaths(values);
if ~isfield(userOpts,'resultFile')
    values.resultFile = defaultPaths.resultFile;
end
if ~isfield(userOpts,'videoFile')
    values.videoFile = defaultPaths.videoFile;
end
if isfield(userOpts,'gridStretchAutomatic')
    automaticGridStretch = values.gridStretchAutomatic;
elseif isfield(userOpts,'gridStretch')
    automaticGridStretch = isempty(userOpts.gridStretch);
else
    automaticGridStretch = true;
end
if automaticGridStretch
    if strcmp(values.gridMode,'uniform')
        values.gridStretch = [0,0];
    else
        targetSpacing = values.targetCenterSpacing(:)';
        if isscalar(targetSpacing)
            targetSpacing = [targetSpacing,targetSpacing];
        end
        if strcmp(values.spatialDiscretization,'sixth_order')
            requestedStretch = [ ...
                automatic_stretch(diff(values.xlim),values.nx-1, ...
                targetSpacing(1)), ...
                automatic_stretch(values.ymax,values.ny-1, ...
                targetSpacing(2))];
            policy = ipm.mesh.sixthOrderPolicy();
            stretchCap = [(values.nx-1)* ...
                log(policy.maximumAdjacentCellRatio)/2, ...
                (values.ny-1)* ...
                log(policy.maximumAdjacentCellRatio)];
            values.gridStretch = min(requestedStretch,stretchCap);
        else
            values.gridStretch = [automatic_stretch( ...
                diff(values.xlim),values.nx-1,targetSpacing(1)), ...
                automatic_stretch(values.ymax,values.ny-1, ...
                targetSpacing(2))];
        end
    end
end
values.gridStretchAutomatic = automaticGridStretch;

validate_numeric_options(values,schema,optionNames);
if numel(values.gridStretch) > 2 || numel(values.targetCenterSpacing) > 2
    error('ipm:GridStretchSize', ...
        'gridStretch and targetCenterSpacing must be scalar or two-vectors.');
end
if mod(values.ccfTailVerticalPower,2) ~= 0
    error('ipm:CcfTailVerticalPower', ...
        'ccfTailVerticalPower must be even so the wall-normal envelope is smooth.');
end
if numel(values.adaptiveLevels) ~= numel(values.targetLevelPoints)
    error('ipm:AdaptiveLevels', ...
        'adaptiveLevels and targetLevelPoints must have the same length.');
end
if any(diff(values.adaptiveLevels) <= 0)
    error('ipm:AdaptiveLevelOrder', ...
        'adaptiveLevels must be strictly increasing.');
end
if numel(values.minimumLevelPoints) ~= numel(values.adaptiveLevels)
    error('ipm:MinimumLevelPoints', ...
        'minimumLevelPoints and adaptiveLevels must have the same length.');
end
validateattributes(values.remeshTargetSafety,{'numeric'}, ...
    {'<',values.remeshSafetyTrigger});
for index = 1:numel(optionNames)
    name = optionNames{index};
    descriptor = schema.options.(name);
    if strcmp(descriptor.kind,'path')
        values.(name) = normalize_path(values.(name),name,descriptor.errorId);
    end
end
metadataDescriptor = schema.options.caseMetadata;
if ~isstruct(values.caseMetadata) || ~isscalar(values.caseMetadata)
    error(metadataDescriptor.errorId,'opts.caseMetadata must be a scalar structure.');
end
if isfield(values,'checkpoint')
    values.checkpoint = normalize_checkpoint_policy( ...
        values.checkpoint,values.resultFile,schema.options.checkpoint.errorId);
end
if isfield(values,'autonomousMesh')
    values.autonomousMesh = ipm.config.autonomousMeshPolicy(values.autonomousMesh,values);
end
if isfield(values,'initialMeshObservationFallback')
    values.initialMeshObservationFallback = ipm.config.initialMeshObservationPolicy( ...
        values.initialMeshObservationFallback,values);
end
quadrantOnly = isfield(values,'quadrantOnly') && values.quadrantOnly;
if quadrantOnly && ( ~strcmp(values.symmetryMode,'double_odd_omega') || ...
        values.xlim(1) ~= 0 || values.xlim(2) <= 0)
    error('ipm:QuadrantDomain', ...
        'quadrantOnly requires double_odd_omega and xlim=[0,H], H>0.');
end
shiftedPrimitive = ischar(values.initialCondition) && any(strcmp( ...
    values.initialCondition,{'shifted_degenerate_primitive', ...
    'shifted_degenerate_primitive_cutoff'}));
if shiftedPrimitive && ...
        (quadrantOnly || ~strcmp(values.symmetryMode,'half_plane') || ...
        values.xlim(1) >= 0 || values.xlim(2) <= 0)
    error('ipm:ShiftedPrimitiveDomain', ...
        ['A shifted degenerate primitive requires a full half-plane ' ...
        'grid with xlim spanning x=0.']);
end
if ischar(values.initialCondition) && ...
        strcmp(values.initialCondition, ...
        'shifted_degenerate_primitive_cutoff') && ...
        values.shiftedCutoffB <= ...
        (values.degeneratePower-1)^(1/values.degeneratePower)
    error('ipm:ShiftedCutoffPeak', ...
        'shiftedCutoffB must exceed the uncut wall-gradient peak coordinate.');
end
if quadrantOnly && isfield(values,'autonomousMesh') && ...
        values.autonomousMesh.enabled
    error('ipm:QuadrantAutonomousMesh', ...
        'The legacy symmetric autonomousMesh policy is not a quadrant policy.');
end
if strcmp(values.symmetryMode,'double_odd_omega') && ~quadrantOnly && ...
        (mod(values.nx,2) ~= 1 || abs(sum(values.xlim)) > ...
        100*eps(max(abs(values.xlim))))
    error('ipm:SymmetryGrid', ...
        'Double-odd omega symmetry requires odd nx and symmetric xlim.');
end
if strcmp(values.dynamicScaleGeometry,'anisotropic') && ...
        ~strcmp(values.rescalingMode,'dynamic')
    error('ipm:AnisotropicRequiresDynamic', ...
        'Anisotropic scale geometry is available only in dynamic mode.');
end
[tupleIdentifier,tupleMessage] = ipm.config.tupleViolation( ...
    values.transportScheme,values.spatialDiscretization, ...
    values.timeIntegrator,values.remeshTransferScheme, ...
    values.dynamicScaleGeometry, ...
    values.anisotropicPoissonSolver);
if ~isempty(tupleIdentifier)
    error(tupleIdentifier,'%s',tupleMessage);
end
if quadrantOnly && ~(strcmp(values.spatialDiscretization,'high_order') && ...
        strcmp(values.transportScheme,'weno5_fd') && ...
        strcmp(values.timeIntegrator,'ssprk54') && ...
        strcmp(values.remeshTransferScheme,'high_order') && ...
        strcmp(values.dynamicScaleGeometry,'isotropic'))
    error('ipm:QuadrantNumericalTuple', ...
        ['quadrantOnly currently supports the high_order/weno5_fd/' ...
        'ssprk54/high_order isotropic tuple.']);
end
if strcmp(values.dynamicScaleGeometry,'anisotropic') && ...
        strcmp(values.symmetryMode,'double_odd_omega') && ...
        values.anisotropicFixedCR ~= 0
    error('ipm:AnisotropicOddTranslation', ...
        'Double-odd omega symmetry requires anisotropicFixedCR=0.');
end
if strcmp(values.dynamicScaleGeometry,'anisotropic') && ...
        strcmp(values.anisotropicGaugeMode,'peak_translation') && ...
        ~strcmp(values.symmetryMode,'half_plane')
    error('ipm:AnisotropicPeakTranslationSymmetry', ...
        'peak_translation is a half-plane gauge; double-odd requires c_r=0.');
end
if strcmp(values.travelingWaveGauge,'gradient_peak') && ...
        ~(strcmp(values.symmetryMode,'half_plane') && ...
        strcmp(values.rescalingMode,'dynamic') && ...
        strcmp(values.dynamicScaleGeometry,'isotropic') && ...
        values.xlim(1) < 0 && values.xlim(2) > 0 && ...
        mod(values.nx,2) == 1)
    error('ipm:GradientPeakGaugeDomain', ...
        ['gradient_peak requires isotropic dynamic rescaling on an odd ' ...
        'full half-plane grid with xlim spanning x=0.']);
end
if strcmp(values.wenoAlphaMode,'local') && ...
        ~strcmp(values.transportScheme,'weno5_fd')
    error('ipm:LocalWenoAlphaMode', ...
        'wenoAlphaMode=local currently requires weno5_fd.');
end
if isfield(values,'wenoInflowDensityDegree') && ...
        (~ismember(values.wenoInflowDensityDegree,[1,2]) || ...
        ~strcmp(values.transportScheme,'weno5_fd') || ...
        ~strcmp(values.transportBoundaryMode,'open') || ...
        ~strcmp(values.symmetryMode,'half_plane'))
    error('ipm:WenoInflowDensityDegree', ...
        ['wenoInflowDensityDegree=1 or 2 requires open WENO5-FD ' ...
        'on the full half-plane.']);
end
if strcmp(values.lengthGauge,'wall_gradient_curvature')
    if ~strcmp(values.scalingContract,'peak_curvature_restoring_v1') || ...
            ~strcmp(values.symmetryMode,'half_plane') || ...
            ~strcmp(values.rescalingMode,'dynamic') || ...
            ~strcmp(values.dynamicScaleGeometry,'isotropic') || ...
            ~strcmp(values.transportScheme,'weno5_fd') || ...
            ~strcmp(values.transportBoundaryMode,'open') || ...
            ~strcmp(values.travelingWaveGauge,'gradient_peak') || ...
            ~strcmp(values.cOmegaGauge,'wall_omega_quadratic_peak')
        error('ipm:PeakCurvatureGaugeTuple', ...
            ['wall_gradient_curvature requires the restoring contract, ' ...
            'half_plane isotropic dynamic scaling, open WENO5-FD, ' ...
            'gradient_peak, and wall_omega_quadratic_peak.']);
    end
    if values.peakCurvatureRadius >= ...
            min(-values.xlim(1),values.xlim(2))
        error('ipm:PeakCurvatureGaugeDomain', ...
            'peakCurvatureRadius must fit on both sides of the origin.');
    end
elseif strcmp(values.scalingContract,'peak_curvature_restoring_v1')
    error('ipm:PeakCurvatureGaugeTuple', ...
        'The restoring contract requires wall_gradient_curvature.');
end
if any(strcmp(values.lengthGauge, ...
        {'omega_peak_location','symmetry_peak', ...
        'transport_anchor'})) && ...
        ~strcmp(values.symmetryMode,'double_odd_omega')
    error('ipm:SymmetryLengthGauge', ...
        'The selected symmetry length gauge requires double_odd_omega symmetry.');
end
if any(strcmp(values.cOmegaGauge, ...
        {'anchor_wall_slope','anchor_wall_strain', ...
        'anchor_wall_window_l2','anchor_wall_template_projection', ...
        'anchor_bulk_gradient_l2','anchor_wall_window_l4'})) && ...
        ~strcmp(values.lengthGauge,'transport_anchor')
    error('ipm:AnchorWallGaugeRequiresTransportAnchor', ...
        ['The anchor-wall c_omega gauges require ' ...
        'opts.lengthGauge=''transport_anchor''.']);
end
if strcmp(values.cOmegaGauge,'anchor_wall_strain')
    error('ipm:RetiredCOmegaGauge', ...
        ['anchor_wall_strain was only an approximate historical relation ' ...
        'and is not available under the exact schema-version-4 contract.']);
end
if any(strcmp(values.cOmegaGauge, ...
        {'anchor_wall_window_l2','anchor_wall_template_projection', ...
        'anchor_bulk_gradient_l2','anchor_wall_window_l4'}))
    anchorTolerance = 100*eps(max(1,abs(values.transportAnchorX)));
    if abs(values.transportAnchorX-1) > anchorTolerance
        error('ipm:OmegaGaugeWindowAnchor', ...
            ['The anchor-wall window c_omega gauge requires ' ...
            'opts.transportAnchorX=1.']);
    end
    windowLeft = values.transportAnchorX-values.omegaGaugeWindowRadius;
    windowRight = values.transportAnchorX+values.omegaGaugeWindowRadius;
    if windowLeft <= 0 || windowRight >= values.xlim(2)
        error('ipm:OmegaGaugeWindowDomain', ...
            ['The open support of the anchor-wall window must lie ' ...
            'strictly in x>0 and strictly inside the positive x domain.']);
    end
    if strcmp(values.cOmegaGauge,'anchor_bulk_gradient_l2') && ...
            values.omegaGaugeWindowRadius >= values.ymax
        error('ipm:OmegaGaugeWindowDomain', ...
            ['The open support of the anchor-bulk window must lie ' ...
            'strictly inside the positive y boundary.']);
    end
end
if strcmp(values.cOmegaGauge,'outer_wall_density_window_l2') && ...
        (2-values.omegaGaugeWindowRadius <= 0 || ...
        2+values.omegaGaugeWindowRadius >= values.xlim(2))
    error('ipm:OmegaGaugeWindowDomain', ...
        ['The X=2 wall-density window must lie strictly in x>0 ' ...
        'and inside the positive x boundary.']);
end
if strcmp(values.cOmegaGauge,'anchor_wall_template_projection') && ...
        values.adaptiveRemesh
    error('ipm:TemplateGaugeRequiresFixedMesh', ...
        ['The frozen wall-template projection is currently defined only ' ...
        'on a fixed mesh; set opts.adaptiveRemesh=false.']);
end
if strcmp(values.lengthGauge,'transport_anchor')
    if values.transportAnchorX >= values.xlim(2)
        error('ipm:TransportAnchorDomain', ...
            'opts.transportAnchorX must lie strictly inside the positive x domain.');
    end
end
if values.minDt >= values.maxDt
    error('ipm:BadTimeStep','opts.minDt must be smaller than opts.maxDt.');
end
if isfinite(values.maxSteps) && values.maxSteps ~= floor(values.maxSteps)
    error('ipm:BadMaximumSteps', ...
        'opts.maxSteps must be a positive integer or Inf.');
end
if isfinite(values.maxRemeshes) && values.maxRemeshes ~= floor(values.maxRemeshes)
    error('ipm:BadMaximumRemeshes', ...
        'opts.maxRemeshes must be a nonnegative integer or Inf.');
end

config = freeze_config(values,schema,optionNames);
end

function values = default_values(schema,optionNames)
values = struct();
for index = 1:numel(optionNames)
    name = optionNames{index};
    descriptor = schema.options.(name);
    if descriptor.hasDefault
        values.(name) = descriptor.default;
    end
end
end

function config = freeze_config(values,schema,optionNames)
config = struct('schemaVersion',4,'frozen',true);
for index = 1:numel(schema.domainNames)
    config.(schema.domainNames{index}) = struct();
end
for index = 1:numel(optionNames)
    name = optionNames{index};
    descriptor = schema.options.(name);
    if ~any(strcmp(descriptor.domain,schema.domainNames))
        error('ipm:UngroupedOption', ...
            'Validated option "%s" lacks a configuration group.',name);
    end
    if isfield(values,name)
        config.(descriptor.domain).(name) = values.(name);
    elseif descriptor.hasDefault
        error('ipm:IncompleteResolvedOptions', ...
            'Resolved option "%s" is missing.',name);
    end
end
unassigned = setdiff(fieldnames(values),optionNames);
if ~isempty(unassigned)
    error('ipm:UngroupedOption', ...
        'Validated option(s) lack a configuration group: %s.', ...
        strjoin(unassigned,', '));
end
end

function validate_numeric_options(values,schema,optionNames)
for index = 1:numel(optionNames)
    name = optionNames{index};
    descriptor = schema.options.(name);
    if ~strcmp(descriptor.kind,'numeric') || ~isfield(values,name)
        continue;
    end
    attributes = schema.rules.(descriptor.rule);
    validateattributes(values.(name),{'numeric'},attributes, ...
        'ipm.config.resolve',['opts.',name]);
end
end

function value = normalize_boolean(value,name,identifier)
validLogical = islogical(value) && isscalar(value);
validNumeric = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && any(value == [0,1]);
if ~(validLogical || validNumeric)
    error(identifier, ...
        'Option "%s" must be logical or numeric 0 or 1.',name);
end
value = logical(value);
end

function value = normalize_enum(value,name,allowed,identifier)
if isstring(value)
    validText = isscalar(value) && ~ismissing(value);
elseif ischar(value)
    validText = isrow(value);
else
    validText = false;
end
if ~validText
    error(identifier,'opts.%s must be one of: %s.',name,quoted_list(allowed));
end
value = lower(char(value));
if ~any(strcmp(value,allowed))
    error(identifier,'opts.%s must be one of: %s.',name,quoted_list(allowed));
end
end

function list = quoted_list(values)
quoted = cellfun(@(value)sprintf('''%s''',value),values, ...
    'UniformOutput',false);
list = strjoin(quoted,', ');
end

function value = normalize_initial_condition(value,descriptor)
if isa(value,'function_handle')
    if ~isscalar(value)
        error(descriptor.errorId, ...
            'opts.initialCondition must be a scalar function handle.');
    end
    return;
end
value = normalize_enum(value,'initialCondition',descriptor.allowed, ...
    descriptor.errorId);
end

function values = apply_custom_domain(values,schema)
hasCustomX = isfield(values,'customX');
hasCustomY = isfield(values,'customY');
if hasCustomX ~= hasCustomY || ...
        (hasCustomX && xor(isempty(values.customX),isempty(values.customY)))
    error('ipm:CustomGridPair', ...
        'customX and customY must be supplied together.');
end
if ~hasCustomX || isempty(values.customX)
    return;
end

validate_numeric_options(values,schema,{'nx','ny'});
validateattributes(values.customX,{'numeric'}, ...
    {'vector','numel',values.nx,'real','finite','increasing'}, ...
    'ipm.config.resolve','opts.customX');
validateattributes(values.customY,{'numeric'}, ...
    {'vector','numel',values.ny,'real','finite','increasing','nonnegative'}, ...
    'ipm.config.resolve','opts.customY');
x = values.customX(:)';
y = values.customY(:);
yTolerance = 100*eps(max(1,max(abs(y))));
if abs(y(1)) > yTolerance
    error('ipm:CustomGridWall', ...
        'customY must start at the wall y=0.');
end
y(1) = 0;

if isfield(values,'quadrantOnly') && values.quadrantOnly
    xTolerance = 100*eps(max(1,max(abs(x))));
    if abs(x(1)) > xTolerance || x(end) <= 0
        error('ipm:CustomGridQuadrant', ...
            'quadrantOnly customX must start at x=0 and end at x>0.');
    end
    x(1) = 0;
elseif strcmp(values.symmetryMode,'double_odd_omega')
    xTolerance = 100*eps(max(1,max(abs(x))));
    mirrorError = max(abs(x+fliplr(x)));
    if mod(values.nx,2) ~= 1 || mirrorError > xTolerance || ...
            abs(x((values.nx+1)/2)) > xTolerance
        error('ipm:CustomGridSymmetry', ...
            ['double_odd_omega requires customX to contain x=0 and be ' ...
            'pointwise mirror symmetric.']);
    end
    % Canonicalize harmless roundoff for downstream exact symmetry logic.
    x = (x-fliplr(x))/2;
    x((values.nx+1)/2) = 0;
end

values.customX = x;
values.customY = y;
values.xlim = x([1,end]);
values.ymax = y(end);
end

function validate_derivation_inputs(values,schema)
% These values feed automatic stretching or filename derivation.
validate_numeric_options(values,schema,{'nx','ny','xlim'});
if values.xlim(2) <= values.xlim(1)
    error('ipm:BadDomain','opts.xlim must be strictly increasing.');
end
validate_numeric_options(values,schema,{'ymax','targetCenterSpacing'});
if numel(values.targetCenterSpacing) > 2
    error('ipm:GridStretchSize', ...
        'targetCenterSpacing must be a scalar or two-vector.');
end
validate_numeric_options(values,schema,{'physicalFinalTime'});
end

function value = normalize_path(value,name,identifier)
if isstring(value)
    validText = isscalar(value) && ~ismissing(value);
elseif ischar(value)
    validText = isrow(value);
else
    validText = false;
end
if ~validText || isempty(char(value))
    error(identifier, ...
        'Option "%s" must be a nonempty character vector or string scalar.',name);
end
value = char(value);
end

function policy = normalize_checkpoint_policy(policy,resultFile,identifier)
if ~isstruct(policy) || ~isscalar(policy)
    error(identifier,'opts.checkpoint must be a scalar structure.');
end
defaults = struct('enabled',true,'file',default_checkpoint_file(resultFile), ...
    'every',Inf,'atExit',true);
names = fieldnames(defaults);
unexpected = setdiff(fieldnames(policy),names,'stable');
if ~isempty(unexpected)
    error(identifier,'Unknown opts.checkpoint field(s): %s.', ...
        strjoin(unexpected,', '));
end
for index = 1:numel(names)
    name = names{index};
    if ~isfield(policy,name)
        policy.(name) = defaults.(name);
    end
end
policy.enabled = normalize_boolean( ...
    policy.enabled,'checkpoint.enabled',identifier);
policy.atExit = normalize_boolean( ...
    policy.atExit,'checkpoint.atExit',identifier);
policy.file = normalize_path( ...
    policy.file,'checkpoint.file',identifier);
validEvery = isnumeric(policy.every) && isreal(policy.every) && ...
    isscalar(policy.every) && ~isnan(policy.every) && policy.every > 0;
if ~validEvery
    error(identifier, ...
        'opts.checkpoint.every must be a positive scalar or Inf.');
end
policy.every = double(policy.every);
end

function fileName = default_checkpoint_file(resultFile)
[folder,stem] = fileparts(char(resultFile));
fileName = fullfile(folder,[stem,'_checkpoint.mat']);
end

function stretch = automatic_stretch(extent,intervals,targetSpacing)
uniformSpacing = extent/intervals;
ratio = min(targetSpacing/uniformSpacing,1);
if ratio >= 1-10*eps
    stretch = 0;
    return;
end
lower = 0;
upper = 32;
for iteration = 1:80
    midpoint = (lower+upper)/2;
    if midpoint/sinh(midpoint) > ratio
        lower = midpoint;
    else
        upper = midpoint;
    end
end
stretch = (lower+upper)/2;
end
