function policy = autonomousMeshPolicy(input,choices)
%IPM.CONFIG.AUTONOMOUSMESHPOLICY Normalize a registered autonomous mesh policy.
% This function is pure configuration: no geometry, fields, or solver calls.
% Mesh counts are admission/telemetry controls, never scaling-rate feedback.
if nargin < 1
    input = struct();
end
identifier = 'ipm:BadAutonomousMeshPolicy';
if ~isstruct(input) || ~isscalar(input)
    error(identifier,'opts.autonomousMesh must be a scalar structure.');
end
version = 1;
if isfield(input,'version')
    value = input.version;
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || ~any(value == [1,2,3,4,5])
        error(identifier,'autonomousMesh.version must be 1, 2, 3, 4, or 5.');
    end
    version = double(value);
end
target = [32,32];
if isfield(input,'targetCoreCells')
    target = input.targetCoreCells;
end
if ~isnumeric(target) || ~isreal(target) || ~isvector(target) || ...
        numel(target) ~= 2 || any(~isfinite(target)) || ...
        any(target < 32 | target > 256 | mod(target,1) ~= 0)
    error(identifier, ...
        'autonomousMesh.targetCoreCells must contain two integers in [32,256].');
end
target = double(target(:)');
quality = struct('maxAdjacentCellRatio',1.08, ...
    'maxLogSpacingCurvature',.01,'minStencilRcond',1e-9, ...
    'minQuadratureWeightRatio',1e-8, ...
    'minWeightToControlWidth',.35,'maxWeightToControlWidth',1.65);
search = struct('fineCells',[32,48,64,80,96,112,128], ...
    'roundingCells',[16,24,32,40,48], ...
    'coreFineCellFractions',[.5,.65],'ySigma',[.18,.25,.35,.5], ...
    'maximumAxisCandidates',70,'maximumPairCandidates',3, ...
    'xPadding',1.10,'yPadding',1.15,'maximumWarpFraction',.25);
defaults = struct('version',version,'enabled',true,'targetCoreCells',target, ...
    'transactionMinimumCoreCells',target-1,'minimumFrontCells',20, ...
    'regridCoreTrigger',(26/32)*target,'predictedCoreBuffer',(22/32)*target, ...
    'endpointCoreFloor',(20/32)*target,'historyCoreFloor',14, ...
    'maximumSafety',.70,'maximumSinglePeakJump',.002, ...
    'maximumCumulativeAbsolutePeakJump',.02, ...
    'maximumMassRelativeDefect',5e-12,'maximumRelativeRangeViolation',2e-4, ...
    'qualityLimits',quality,'search',search,'timeUnit','native_canonical', ...
    'trendWindow',.35,'maximumReviewInterval',.2);
% Explicit, frozen density opt-in. Absence preserves every prior v5 policy
% and checkpoint exactly; no implicit upgrade occurs during restoration.
if isfield(input,'densityVersion')
    if version~=5 || ~isnumeric(input.densityVersion) || ...
            ~isscalar(input.densityVersion) || ~isreal(input.densityVersion) || ...
            ~isfinite(input.densityVersion) || ~any(input.densityVersion==[1,2])
        error(identifier,'densityVersion=1 or 2 is supported only with autonomousMesh.version=5.');
    end
    defaults.densityVersion=double(input.densityVersion);
    defaults.search.roundingCells=[8,16,24,32,40];
    if defaults.densityVersion==2
        % More same-node core reserve, still under the original whole-axis
        % quality, transfer, and finite-flow gates.
        defaults.search.xPadding=1.50;
    end
end
if any(version == [2,3,4,5])
    if ~isfield(input,'nodeFamily')
        error(identifier,'Versions 2 through 5 require an explicit nodeFamily.maximumTotalNodes.');
    end
    defaults.nodeFamily = node_family(input.nodeFamily,identifier,version);
    if version == 5 && nargin >= 2 && isfield(choices,'nx') && isfield(choices,'ny') && ...
            ~isequal(defaults.nodeFamily.initialNodeCount,double([choices.nx,choices.ny]))
        error(identifier,'Version-5 nodeFamily.initialNodeCount must match the configured initial grid.');
    end
end
if any(version == [4,5])
    if version == 5
        members=size(defaults.nodeFamily.cellFactors,1);
        defaults.axisSearchPolicy = ipm.remesh.searchEvidence('registration',members,5);
    else
        defaults.axisSearchPolicy = ipm.remesh.searchEvidence('registration');
    end
end
names = fieldnames(defaults);
reject_unknown(input,names,'autonomousMesh',identifier);
policy = defaults;
for index = 1:numel(names)
    name = names{index};
    if ~isfield(input,name)
        continue;
    end
    value = input.(name);
    label = ['autonomousMesh.',name];
    if strcmp(name,'enabled')
        valid = (islogical(value) && isscalar(value)) || ...
            (isnumeric(value) && isreal(value) && isscalar(value) && ...
            isfinite(value) && any(value == [0,1]));
        if ~valid
            error(identifier,'%s must be logical or numeric 0 or 1.',label);
        end
        policy.enabled = logical(value);
    elseif strcmp(name,'qualityLimits') || strcmp(name,'search')
        policy.(name) = fixed_group(value,defaults.(name),label,identifier);
    elseif strcmp(name,'axisSearchPolicy')
        if ~isequaln(value,defaults.axisSearchPolicy)
            error(identifier,'axisSearchPolicy must equal the registered descriptor.');
        end
        policy.axisSearchPolicy=defaults.axisSearchPolicy;
    elseif strcmp(name,'nodeFamily')
        policy.nodeFamily = defaults.nodeFamily;
    elseif strcmp(name,'timeUnit')
        if isstring(value) && isscalar(value) && ~ismissing(value)
            value = char(value);
        end
        if ~ischar(value) || ~isrow(value) || ...
                ~strcmpi(value,defaults.timeUnit)
            error(identifier,'%s must be native_canonical.',label);
        end
        policy.timeUnit = defaults.timeUnit;
    else
        policy.(name) = fixed_numeric(value,defaults.(name),label,identifier);
    end
end
if nargin >= 2 && policy.enabled
    validate_combination(choices);
end
end

function group = node_family(input,identifier,version)
label = 'autonomousMesh.nodeFamily';
if ~isstruct(input) || ~isscalar(input) || ~isfield(input,'maximumTotalNodes')
    error(identifier,'%s requires an explicit maximumTotalNodes.',label);
end
cap = input.maximumTotalNodes;
if ~isnumeric(cap) || ~isreal(cap) || ~isscalar(cap) || ...
        ~isfinite(cap) || cap <= 0 || mod(cap,1) ~= 0
    error(identifier,'%s.maximumTotalNodes must be a positive finite integer.',label);
end
group = struct('generator','selected_base_index_pchip_v1', ...
    'cellFactors',[1,1;2,1;1,2;2,2], ...
    'ordering','node_product_then_registration_index', ...
    'componentwiseNondecreasing',true, ...
    'maximumAcceptedGrowthTransitions',2,'maximumTotalNodes',double(cap));
if any(version == [3,4])
    if cap ~= 310000
        error(identifier,'Versions 3 and 4 register maximumTotalNodes=310000 explicitly.');
    end
    group.generator = 'selected_base_index_pchip_integer_v2';
    group.cellFactors = [1,1;2,1;1,2;2,2;3,1;1,3;3,2;2,3];
    group.maximumAcceptedGrowthTransitions = 3;
end
if version == 5
    if ~isfield(input,'initialNodeCount')
        error(identifier,'%s.initialNodeCount is required for version 5.',label);
    end
    initial=input.initialNodeCount;
    if ~isnumeric(initial)||~isreal(initial)||~isequal(size(initial),[1,2])|| ...
            any(~isfinite(initial))||any(initial<17)||any(mod(initial,1)~=0)|| ...
            mod(initial(1),2)~=1
        error(identifier,'%s.initialNodeCount must be [odd Nx, Ny] with both >=17.',label);
    end
    initial=double(initial);
    if cap<prod(initial)
        error(identifier,'%s.maximumTotalNodes must cover the initial grid.',label);
    end
    xFactors=factor_schedule(floor((cap/initial(2)-1)/(initial(1)-1)));
    yFactors=factor_schedule(floor((cap/initial(1)-1)/(initial(2)-1)));
    rows=zeros(numel(xFactors)*numel(yFactors),4);count=0;
    for fx=xFactors
        for fy=yFactors
            nodes=[(initial(1)-1)*fx+1,(initial(2)-1)*fy+1];
            if prod(nodes)>cap,continue;end
            count=count+1;rows(count,:)=[fx,fy,prod(nodes),count];
        end
    end
    rows=sortrows(rows(1:count,:),[3,4]);
    if count>192
        error(identifier,'%s creates more than 192 registered levels.',label);
    end
    group.generator='selected_base_index_pchip_geometric_budget_v3';
    group.cellFactors=rows(:,1:2);
    group.ordering='node_product_then_registration_index';
    group.maximumAcceptedGrowthTransitions=numel(xFactors)+numel(yFactors)-2;
    group.initialNodeCount=initial;
end
names = fieldnames(group);
reject_unknown(input,names,label,identifier);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(input,name) || strcmp(name,'maximumTotalNodes')
        continue;
    end
    value = input.(name); expected = group.(name);
    if ischar(expected)
        if isstring(value) && isscalar(value) && ~ismissing(value)
            value = char(value);
        end
        valid = ischar(value) && isrow(value) && strcmpi(value,expected);
    elseif islogical(expected)
        valid = (islogical(value) || isnumeric(value)) && isreal(value) && ...
            isscalar(value) && isfinite(value) && value == 1;
    else
        valid = isnumeric(value) && isreal(value) && ...
            isequal(size(value),size(expected)) && all(isfinite(value(:))) && ...
            isequal(double(value),expected);
    end
    if ~valid
        error(identifier,'%s.%s differs from its registered version-two value.',label,name);
    end
end
end

function factors=factor_schedule(maximum)
% Every integer through six, then bounded multiplicative growth. The
% resource cap, not a fixed family endpoint, terminates this registration.
factors=1:min(maximum,6);
while factors(end)<maximum
    next=min(maximum,max(factors(end)+1,ceil(4*factors(end)/3)));
    factors(end+1)=next; %#ok<AGROW>
end
end

function validate_combination(values)
required = {'rescalingMode','dynamicScaleGeometry','symmetryMode', ...
    'lengthGauge','cOmegaGauge','spatialDiscretization','transportScheme', ...
    'timeIntegrator','remeshTransferScheme','adaptiveRemesh','initialAnalyticRemesh', ...
    'adaptiveLevels'};
supported = isstruct(values) && isscalar(values) && all(isfield(values,required));
if supported
    supported = strcmp(values.rescalingMode,'dynamic') && ...
        strcmp(values.dynamicScaleGeometry,'isotropic') && ...
        strcmp(values.symmetryMode,'double_odd_omega') && ...
        strcmp(values.lengthGauge,'transport_anchor') && ...
        any(strcmp(values.cOmegaGauge, ...
            {'wall_omega_quadratic_peak', ...
            'outer_wall_density_window_l2'})) && ...
        strcmp(values.spatialDiscretization,'high_order') && ...
        strcmp(values.transportScheme,'weno5_fd') && ...
        strcmp(values.timeIntegrator,'ssprk54') && ...
        strcmp(values.remeshTransferScheme,'high_order') && ...
        isscalar(values.adaptiveRemesh) && isequal(values.adaptiveRemesh,true) && ...
        isscalar(values.initialAnalyticRemesh) && isequal(values.initialAnalyticRemesh,true) && ...
        isequal(values.adaptiveLevels(:)',[.1,.5,.9]);
end
if ~supported
    error('ipm:AutonomousMeshCombination', ...
        ['Enabled autonomousMesh requires dynamic/isotropic/double_odd_omega, ' ...
        'transport_anchor with a registered quadratic-peak or ' ...
        'X=2 outer wall-density window amplitude gauge, ' ...
        'high_order/weno5_fd/ssprk54/high_order transfer, ' ...
        'adaptiveRemesh=true, initialAnalyticRemesh=true, and adaptiveLevels=[.1 .5 .9].']);
end
end

function group = fixed_group(input,defaults,label,identifier)
if ~isstruct(input) || ~isscalar(input)
    error(identifier,'%s must be a scalar structure.',label);
end
names = fieldnames(defaults);
reject_unknown(input,names,label,identifier);
group = defaults;
for index = 1:numel(names)
    name = names{index};
    if isfield(input,name)
        group.(name) = fixed_numeric( ...
            input.(name),defaults.(name),[label,'.',name],identifier);
    end
end
end

function value = fixed_numeric(value,expected,label,identifier)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        numel(value) ~= numel(expected) || any(~isfinite(value))
    error(identifier,'%s must have its registered finite numeric shape.',label);
end
value = double(value(:)');
if ~isequal(value,expected)
    error(identifier, ...
        ['%s differs from the registered common value or the ' ...
        'exact targetCoreCells-derived threshold.'],label);
end
end

function reject_unknown(input,names,label,identifier)
unexpected = setdiff(fieldnames(input),names,'stable');
if ~isempty(unexpected)
    error(identifier,'Unknown %s field(s): %s.',label,strjoin(unexpected,', '));
end
end
