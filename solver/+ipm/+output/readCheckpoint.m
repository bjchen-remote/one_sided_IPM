function checkpoint = readCheckpoint(inputCheckpoint)
%IPM.OUTPUT.READCHECKPOINT Read and validate one schema-4 checkpoint.

if ischar(inputCheckpoint) || ...
        (isstring(inputCheckpoint) && isscalar(inputCheckpoint))
    fileName = char(inputCheckpoint);
    if ~isfile(fileName)
        error('ipm:CheckpointMissing', ...
            'Checkpoint file does not exist: %s',fileName);
    end
    loaded = load(fileName,'checkpoint');
    if ~isfield(loaded,'checkpoint')
        error('ipm:CheckpointFile', ...
            'Checkpoint MAT-file must contain a variable named checkpoint.');
    end
    checkpoint = loaded.checkpoint;
else
    checkpoint = inputCheckpoint;
end
validate_checkpoint(checkpoint);
end

function validate_checkpoint(checkpoint)
required = {'schemaVersion','kind','createdUtc','acceptedStep', ...
    'trustedRecord','payload','signature'};
if ~isstruct(checkpoint) || ~isscalar(checkpoint)
    error('ipm:CheckpointContract', ...
        'Checkpoint must be a scalar schema-version-4 structure.');
end
if ~isfield(checkpoint,'schemaVersion') || ...
        ~isnumeric(checkpoint.schemaVersion) || ...
        ~isreal(checkpoint.schemaVersion) || ...
        ~isscalar(checkpoint.schemaVersion) || ...
        ~isfinite(checkpoint.schemaVersion)
    error('ipm:CheckpointContract', ...
        'Checkpoint must declare numeric schemaVersion=4.');
end
if checkpoint.schemaVersion <= 3
    error('ipm:CheckpointIncompatibleSchema', ...
        ['Checkpoint schema versions 1--3 predate the strict no-feedback ' ...
        'scaling contract and cannot be resumed by the current solver. ' ...
        'They remain historical artifacts; use matching legacy code.']);
end
if checkpoint.schemaVersion ~= 4
    error('ipm:CheckpointContract', ...
        'Only accepted-step checkpoint schemaVersion=4 is supported.');
end
if ~isempty(setdiff(required,fieldnames(checkpoint),'stable')) || ...
        ~strcmp(checkpoint.kind,'ipm_accepted_step_checkpoint') || ...
        ~(islogical(checkpoint.acceptedStep) && ...
        isscalar(checkpoint.acceptedStep) && checkpoint.acceptedStep) || ...
        ~(islogical(checkpoint.trustedRecord) && ...
        isscalar(checkpoint.trustedRecord) && checkpoint.trustedRecord)
    error('ipm:CheckpointContract', ...
        'Checkpoint does not satisfy the accepted-step schema-version-4 contract.');
end
if ~(ischar(checkpoint.createdUtc) && isrow(checkpoint.createdUtc) && ...
        ~isempty(checkpoint.createdUtc))
    error('ipm:CheckpointContract','Checkpoint createdUtc must be text.');
end
payload = checkpoint.payload;
if ~isstruct(payload) || ~isscalar(payload) || ...
        ~all(isfield(payload,{'state','log','cursor'}))
    error('ipm:CheckpointContract', ...
        'Checkpoint payload must contain state, log, and cursor.');
end
validate_state_payload(payload.state);
validate_log(payload.log,payload.state);
validate_controller(payload.state,payload.log);
validate_cursor(payload.cursor,payload.state);
actualSignature = ipm.output.checkpointSignature(payload);
if ~isequaln(actualSignature,checkpoint.signature)
    error('ipm:CheckpointSignature', ...
        'Checkpoint integrity signature does not match its payload.');
end
end

function validate_controller(state,log)
hasPolicy = isfield(state.config.remesh,'autonomousMesh');
if ~hasPolicy && ~isfield(state.runMetadata,'autonomousMesh')
    return;
end
current = struct();
if hasPolicy && state.config.remesh.autonomousMesh.enabled
    history=log.history;
    current = struct('step',state.step,'canonicalTime',state.scale.canonicalTime, ...
        'physicalTime',state.scale.physicalTime,'normalizedTime',state.normalizedTime, ...
        'remeshCount',state.remeshCount, ...
        'coreCells',[history.mesh.coreGridPoints(end),history.mesh.verticalCoreGridPoints(end)], ...
        'safety',history.mesh.safetyFactor(end),'x',state.x,'y',state.y, ...
        'baseX',state.baseX,'baseY',state.baseY,'history',history);
    if any(state.config.remesh.autonomousMesh.version == [2,3,4,5])
        if ~isfield(state.runMetadata,'autonomousMesh') || ...
                ~isfield(state.runMetadata.autonomousMesh,'currentLevelId')
            error('ipm:AutonomousMeshController','Version-two checkpoints require the original current level.');
        end
        current.levelId=state.runMetadata.autonomousMesh.currentLevelId;
        current.nodeCount=[numel(state.x),numel(state.y)];
        current.snapshots=struct('rho',{log.snapshotRho},'x',{log.snapshotX},'y',{log.snapshotY});
    end
end
ipm.remesh.validateController(state.runMetadata,state.config,current);
end

function validate_state_payload(state)
required = {'config','runMetadata','rho','x','y','baseX','baseY', ...
    'rescaling','remeshCount','scale','normalizedTime','step', ...
    'mass0','rhoRange0'};
if ~isstruct(state) || ~isscalar(state) || ...
        ~all(isfield(state,required))
    error('ipm:CheckpointState', ...
        'Checkpoint state payload is incomplete.');
end
if isfield(state,'rhsCache')
    error('ipm:CheckpointDerivedCache', ...
        ['Checkpoint state must not persist the transient RHS cache; ' ...
        'restore reconstructs it from the accepted state.']);
end
validateattributes(state.x,{'numeric'}, ...
    {'vector','real','finite','increasing'},mfilename,'state.x');
validateattributes(state.y,{'numeric'}, ...
    {'vector','real','finite','increasing','nonnegative'},mfilename,'state.y');
validateattributes(state.baseX,{'numeric'}, ...
    {'vector','real','finite','increasing','numel',numel(state.x)}, ...
    mfilename,'state.baseX');
validateattributes(state.baseY,{'numeric'}, ...
    {'vector','real','finite','increasing','nonnegative','numel',numel(state.y)}, ...
    mfilename,'state.baseY');
if ~isnumeric(state.rho) || ~isreal(state.rho) || ...
        any(~isfinite(state.rho),'all') || ...
        ~isequal(size(state.rho),[numel(state.y),numel(state.x)])
    error('ipm:CheckpointState', ...
        'Checkpoint rho must be finite and match its x/y axes.');
end
if state.y(1) ~= 0 || state.baseY(1) ~= 0 || ...
        ~same_endpoints(state.x,state.baseX) || ...
        ~same_endpoints(state.y,state.baseY)
    error('ipm:CheckpointGrid', ...
        'Checkpoint current/base axes must share endpoints and the wall y=0.');
end
validateattributes(state.remeshCount,{'numeric'}, ...
    {'scalar','real','finite','integer','nonnegative'},mfilename, ...
    'state.remeshCount');
validateattributes(state.normalizedTime,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,'state.normalizedTime');
validateattributes(state.step,{'numeric'}, ...
    {'scalar','real','finite','integer','nonnegative'},mfilename,'state.step');
validateattributes(state.mass0,{'numeric'}, ...
    {'scalar','real','finite'},mfilename,'state.mass0');
validateattributes(state.rhoRange0,{'numeric'}, ...
    {'vector','numel',2,'real','finite'},mfilename,'state.rhoRange0');
if state.rhoRange0(2) < state.rhoRange0(1)
    error('ipm:CheckpointState','Checkpoint initial density range is invalid.');
end
if ~isstruct(state.runMetadata) || ~isscalar(state.runMetadata) || ...
        ~isstruct(state.rescaling) || ~isscalar(state.rescaling)
    error('ipm:CheckpointState', ...
        'Checkpoint metadata and rescaling data must be scalar structures.');
end
retiredFeedback = {'maxDynamicRate','adaptiveLengthScaling', ...
    'widthExpansionStrength', ...
    'widthContractionOnset','widthContractionStrength', ...
    'maxWidthRateCorrection'};
presentRetired = retiredFeedback(isfield(state.rescaling,retiredFeedback));
if ~isempty(presentRetired)
    error('ipm:CheckpointRescaling', ...
        ['Schema-version-4 checkpoint rescaling state contains retired ' ...
        'feedback field %s.'],presentRetired{1});
end
config = validate_config(state.config);
variableNodes=isfield(config.remesh,'autonomousMesh') && ...
    config.remesh.autonomousMesh.enabled && any(config.remesh.autonomousMesh.version == [2,3,4,5]);
if ~variableNodes && (config.grid.nx ~= numel(state.x) || config.grid.ny ~= numel(state.y))
    error('ipm:CheckpointGrid', ...
        'Checkpoint axes do not match the configured node counts.');
end
validate_scale(state.scale,config.scaling.dynamicScaleGeometry);
end

function config = validate_config(config)
if ~isstruct(config) || ~isscalar(config) || ...
        ~isfield(config,'schemaVersion') || config.schemaVersion ~= 4 || ...
        ~isfield(config,'frozen') || ~isequal(config.frozen,true)
    error('ipm:CheckpointConfig', ...
        'Checkpoint config must be a frozen schema-version-4 configuration.');
end
schema = ipm.config.schema();
flat = struct();
for domainIndex = 1:numel(schema.domainNames)
    domainName = schema.domainNames{domainIndex};
    if ~isfield(config,domainName) || ~isstruct(config.(domainName)) || ...
            ~isscalar(config.(domainName))
        error('ipm:CheckpointConfig', ...
            'Checkpoint config is missing domain %s.',domainName);
    end
    names = fieldnames(config.(domainName));
    for nameIndex = 1:numel(names)
        flat.(names{nameIndex}) = config.(domainName).(names{nameIndex});
    end
end
try
    resolved = ipm.config.resolve(flat);
catch exception
    error('ipm:CheckpointConfig', ...
        'Checkpoint config is invalid: %s',exception.message);
end
if ~isequaln(orderfields(resolved),orderfields(config))
    error('ipm:CheckpointConfig', ...
        'Checkpoint config is not a canonical resolved configuration.');
end
end

function validate_scale(scale,geometry)
if ~isstruct(scale) || ~isscalar(scale)
    error('ipm:CheckpointScale','Checkpoint scale must be a scalar structure.');
end
if strcmp(geometry,'anisotropic')
    names = {'logC_l','logC_x','logC_y','logC_omega', ...
        'physicalTime','X_shift','canonicalTime'};
else
    names = {'logC_l','logC_omega','physicalTime','X_shift','canonicalTime'};
end
if ~all(isfield(scale,names))
    error('ipm:CheckpointScale','Checkpoint scale fields are incomplete.');
end
for index = 1:numel(names)
    validateattributes(scale.(names{index}),{'numeric'}, ...
        {'scalar','real','finite'},mfilename,['state.scale.',names{index}]);
end
if scale.physicalTime < 0 || scale.canonicalTime < 0 || ...
        (strcmp(geometry,'anisotropic') && scale.logC_x ~= scale.logC_l)
    error('ipm:CheckpointScale','Checkpoint scale clocks or geometry are invalid.');
end
end

function validate_log(log,state)
required = {'history','snapshotRho','snapshotNormalizedTime', ...
    'snapshotX','snapshotY'};
if ~isstruct(log) || ~isscalar(log) || ~all(isfield(log,required)) || ...
        ~isstruct(log.history) || ~isscalar(log.history)
    error('ipm:CheckpointLog','Checkpoint log is incomplete.');
end
groups = {'common','gauge','mesh','anisotropic'};
if ~all(isfield(log.history,groups)) || ...
        ~isfield(log.history.common,'t') || ...
        ~isfield(log.history.common,'canonicalTau') || ...
        ~isfield(log.history.common,'physicalTime')
    error('ipm:CheckpointLog','Checkpoint history groups are incomplete.');
end
count = aligned_group_count(log.history.common,'common');
if count < 1
    error('ipm:CheckpointLog','Checkpoint history must not be empty.');
end
for index = 2:numel(groups)
    groupCount = aligned_group_count(log.history.(groups{index}),groups{index});
    if groupCount ~= 0 && groupCount ~= count
        error('ipm:CheckpointLog','Checkpoint history groups are misaligned.');
    end
end
t = log.history.common.t(:);
canonical = log.history.common.canonicalTau(:);
physical = log.history.common.physicalTime(:);
if any(~isfinite(t)) || any(~isfinite(canonical)) || any(~isfinite(physical)) || ...
        any(diff(t) <= 0) || any(diff(canonical) <= 0) || ...
        any(diff(physical) < 0)
    error('ipm:CheckpointLog','Checkpoint history clocks are invalid.');
end
assert_terminal_clock(t(end),state.normalizedTime,'normalized');
assert_terminal_clock(canonical(end),state.scale.canonicalTime,'canonical');
assert_terminal_clock(physical(end),state.scale.physicalTime,'physical');
validate_canonical_identity(log.history.common,state,t,canonical);
if isfield(log.history,'wallCore')
    ipm.output.validateWallCoreHistory( ...
        log.history.wallCore,log.history.common);
end

snapshotCount = numel(log.snapshotRho);
if ~iscell(log.snapshotRho) || ~iscell(log.snapshotX) || ...
        ~iscell(log.snapshotY) || numel(log.snapshotX) ~= snapshotCount || ...
        numel(log.snapshotY) ~= snapshotCount || ...
        ~isnumeric(log.snapshotNormalizedTime) || ...
        numel(log.snapshotNormalizedTime) ~= snapshotCount
    error('ipm:CheckpointLog','Checkpoint snapshot arrays are misaligned.');
end
if state.config.output.storeSnapshots
    if snapshotCount ~= count
        error('ipm:CheckpointLog', ...
            'Stored checkpoint snapshots must match every history row.');
    end
else
    if snapshotCount ~= 0
        error('ipm:CheckpointLog', ...
            'A storeSnapshots=false checkpoint must not carry snapshot arrays.');
    end
end
if snapshotCount > 0 && ...
        ~isequal(log.snapshotNormalizedTime(:),t)
    error('ipm:CheckpointLog', ...
        'Checkpoint snapshot times must match the recorded history.');
end
for index = 1:snapshotCount
    x = log.snapshotX{index};
    y = log.snapshotY{index};
    rho = log.snapshotRho{index};
    if ~isnumeric(rho) || ~isreal(rho) || any(~isfinite(rho),'all') || ...
            ~isequal(size(rho),[numel(y),numel(x)])
        error('ipm:CheckpointLog','Checkpoint snapshot %d is invalid.',index);
    end
end
trusted = ipm.output.continuousTrustedPrefix( ...
    ipm.output.trustedMask(log.history,state.config));
if ~trusted(end)
    error('ipm:CheckpointTrust', ...
        'Checkpoint terminal record is outside the continuous trusted prefix.');
end
end

function validate_canonical_identity(common,state,t,canonical)
if ~isfield(common,'timeSpeed')
    error('ipm:CheckpointCanonicalClock', ...
        'Schema-version-4 checkpoint history must record timeSpeed.');
end
assert_checkpoint_identity(t,canonical, ...
    'normalized and canonical history clocks', ...
    'ipm:CheckpointCanonicalClock');
assert_checkpoint_identity(state.normalizedTime,state.scale.canonicalTime, ...
    'normalized and canonical stored-state clocks', ...
    'ipm:CheckpointCanonicalClock');
assert_checkpoint_identity(common.timeSpeed,ones(size(common.timeSpeed)), ...
    'unit canonical time speed','ipm:CheckpointCanonicalClock');

directNames = {'c_l','c_x','c_y','c_omega','c_r', ...
    'conservativeSource'};
canonicalNames = {'canonicalCL','canonicalCX','canonicalCY', ...
    'canonicalCOmega','canonicalCR','canonicalConservativeSource'};
for index = 1:numel(directNames)
    if ~isfield(common,directNames{index}) || ...
            ~isfield(common,canonicalNames{index})
        error('ipm:CheckpointCanonicalRate', ...
            ['Schema-version-4 checkpoint history must record both ' ...
            '%s and %s.'],directNames{index},canonicalNames{index});
    end
    assert_checkpoint_identity(common.(directNames{index}), ...
        common.(canonicalNames{index}), ...
        [directNames{index},' and ',canonicalNames{index}], ...
        'ipm:CheckpointCanonicalRate');
end
end

function assert_checkpoint_identity(actual,expected,name,identifier)
if ~isnumeric(actual) || ~isnumeric(expected) || ...
        ~isreal(actual) || ~isreal(expected) || ...
        ~isvector(actual) || ~isvector(expected) || ...
        numel(actual) ~= numel(expected) || isempty(actual)
    error(identifier, ...
        'Checkpoint fields for %s must be aligned numeric vectors.',name);
end
actual = actual(:);
expected = expected(:);
if any(~isfinite(actual)) || any(~isfinite(expected))
    error(identifier,'Checkpoint fields for %s must be finite.',name);
end
scale = max([ones(size(actual)),abs(actual),abs(expected)],[],2);
if any(abs(actual-expected) > 100*eps(scale))
    error(identifier, ...
        'Checkpoint fields for %s violate the canonical identity.',name);
end
end

function count = aligned_group_count(group,name)
if ~isstruct(group) || ~isscalar(group)
    error('ipm:CheckpointLog','Checkpoint history.%s must be scalar.',name);
end
names = fieldnames(group);
if isempty(names)
    count = 0;
    return;
end
count = numel(group.(names{1}));
for index = 1:numel(names)
    values = group.(names{index});
    if ~isvector(values) || numel(values) ~= count
        error('ipm:CheckpointLog', ...
            'Checkpoint history.%s fields must be aligned vectors.',name);
    end
end
end

function validate_cursor(cursor,state)
required = {'nextOutput','nextCheckpoint','lastCheckpointStep', ...
    'lastCheckpointCanonicalTime'};
if ~isstruct(cursor) || ~isscalar(cursor) || ...
        ~isempty(setxor(fieldnames(cursor),required))
    error('ipm:CheckpointCursor','Checkpoint cursor is incomplete.');
end
validateattributes(cursor.nextOutput,{'numeric'}, ...
    {'scalar','real','positive','finite'},mfilename,'cursor.nextOutput');
validateattributes(cursor.nextCheckpoint,{'numeric'}, ...
    {'scalar','real','positive','nonnan'},mfilename,'cursor.nextCheckpoint');
validateattributes(cursor.lastCheckpointStep,{'numeric'}, ...
    {'scalar','real','integer','finite','nonnegative'},mfilename, ...
    'cursor.lastCheckpointStep');
validateattributes(cursor.lastCheckpointCanonicalTime,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename, ...
    'cursor.lastCheckpointCanonicalTime');
if cursor.lastCheckpointStep ~= state.step || ...
        cursor.lastCheckpointCanonicalTime ~= state.scale.canonicalTime
    error('ipm:CheckpointCursor', ...
        'Checkpoint cursor does not identify the stored accepted state.');
end
end

function answer = same_endpoints(first,second)
scale = max([1,abs(first(1)),abs(first(end)), ...
    abs(second(1)),abs(second(end))]);
answer = max(abs([first(1)-second(1),first(end)-second(end)])) <= ...
    64*eps(scale);
end

function assert_terminal_clock(historyValue,stateValue,name)
tolerance = 64*eps(max([1,abs(historyValue),abs(stateValue)]));
if abs(historyValue-stateValue) > tolerance
    error('ipm:CheckpointClock', ...
        'Checkpoint %s history does not end at the stored state.',name);
end
end
