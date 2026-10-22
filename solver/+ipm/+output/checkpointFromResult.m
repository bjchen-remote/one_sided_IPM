function [checkpoint,fileName] = checkpointFromResult(inputResult,requestedFile)
%IPM.OUTPUT.CHECKPOINTFROMRESULT Convert a trusted terminal v2 result whose
%   frozen configuration uses the current schema-version-4 equations.
%   CHECKPOINT = IPM.OUTPUT.CHECKPOINTFROMRESULT(RESULT) reconstructs the
%   same-grid evolution state and returns a validated in-memory checkpoint.
%   [CHECKPOINT,FILE] = ...(...,REQUESTEDFILE) atomically installs an
%   immutable checkpoint beside REQUESTEDFILE. This bridge never regrids.

if nargin < 2
    requestedFile = [];
end
result = ipm.output.validate(inputResult);
if isfield(result.config.remesh,'autonomousMesh') && ...
        result.config.remesh.autonomousMesh.enabled && any(result.config.remesh.autonomousMesh.version == [2,3,4,5])
    error('ipm:CheckpointResultVariableNodeFamily', ...
        ['Version-two autonomous runs require their native accepted-step checkpoint. ' ...
        'Result conversion does not reconstruct a registered variable-node family.']);
end
if result.config.schemaVersion ~= 4
    error('ipm:CheckpointResultConfigVersion', ...
        ['Legacy result configurations remain readable for analysis, but ' ...
        'only config schema 4 implements the strict no-feedback scaling ' ...
        'contract required by current checkpoints.']);
end
trusted = ipm.output.continuousTrustedPrefix( ...
    ipm.output.trustedMask(result));
if isempty(trusted) || ~trusted(end)
    error('ipm:CheckpointResultTrust', ...
        ['Only a result whose terminal record lies in the continuous ' ...
        'trusted prefix can be converted to a restart checkpoint.']);
end
assert_initial_record(result);

seed = rebuild_seed(result.config);
validate_initial_overlap(seed,result);
ops = rebuild_terminal_ops(seed,result);
scale = rebuild_scale(result);
[~,flow] = ipm.evolve.flow(result.state.rho,ops,scale);
validate_flow_overlap(flow,result.state);

metadata = result.metadata;
metadata.checkpointSource = "trusted_terminal_result_v2";
metadata.checkpointSourceStopReason = string(result.state.stopReason);
if ischar(inputResult) || (isstring(inputResult) && isscalar(inputResult))
    metadata.checkpointSourceResult = string(inputResult);
else
    metadata.checkpointSourceResult = "in_memory_result";
end
state = struct('config',result.config,'runMetadata',metadata, ...
    'ops',ops,'rho',result.state.rho,'flow',flow,'scale',scale, ...
    'normalizedTime',result.state.normalizedTime, ...
    'step',result.state.steps, ...
    'timeStep',terminal_timestep(result,scale), ...
    'mass0',seed.mass0, ...
    'rhoRange0',seed.rhoRange0);
if ~ipm.evolve.isFinite(state)
    error('ipm:CheckpointResultState', ...
        'The terminal result could not be reconstructed as a finite state.');
end
validate_terminal_overlap(state,result);

log = result_log(result);
cursor = result_cursor(result);
checkpoint = ipm.output.makeCheckpoint(state,log,cursor);
fileName = "";
if ~isempty(requestedFile)
    [checkpoint,fileName] = ...
        ipm.output.writeCheckpoint(checkpoint,requestedFile);
end
end

function assert_initial_record(result)
common = result.history.common;
scale = max([1,abs(common.t(1)),abs(common.canonicalTau(1)), ...
    abs(common.physicalTime(1))]);
if max(abs([common.t(1),common.canonicalTau(1), ...
        common.physicalTime(1)])) > 64*eps(scale)
    error('ipm:CheckpointResultInitialRecord', ...
        ['Terminal-result conversion requires the original time-zero ' ...
        'history record so initial invariants can be reconstructed.']);
end
end

function seed = rebuild_seed(config)
flat = flatten_config(config);
% Initialization has no output side effects, but silence its optional
% analytic-remesh report and avoid carrying a checkpoint policy into it.
flat.verbose = false;
flat.saveResults = false;
flat.makePlots = false;
flat.livePlot = false;
flat.writeVideo = false;
flat.storeSnapshots = false;
if isfield(flat,'checkpoint')
    flat = rmfield(flat,'checkpoint');
end
seed = ipm.evolve.initialize(flat);
end

function flat = flatten_config(config)
schema = ipm.config.schema();
flat = struct();
for domainIndex = 1:numel(schema.domainNames)
    domain = config.(schema.domainNames{domainIndex});
    names = fieldnames(domain);
    for nameIndex = 1:numel(names)
        flat.(names{nameIndex}) = domain.(names{nameIndex});
    end
end
end

function validate_initial_overlap(seed,result)
common = result.history.common;
assert_numeric_overlap(seed.mass0,common.physicalMass(1), ...
    'initial mass',2e-11);
assert_numeric_overlap(seed.rhoRange0(1),common.physicalRhoMin(1), ...
    'initial density minimum',2e-11);
assert_numeric_overlap(seed.rhoRange0(2),common.physicalRhoMax(1), ...
    'initial density maximum',2e-11);

probeState = seed;
probeState.config.output.verbose = false;
probeState.config.output.storeSnapshots = false;
[probeLog,~] = ipm.output.record(ipm.output.initializeLog(),probeState);
compare_history_row(probeLog.history,result.history,1,'initial');
end

function ops = rebuild_terminal_ops(seed,result)
x = result.grid.x;
y = result.grid.y;
if isequal(seed.ops.x,x) && isequal(seed.ops.y,y)
    ops = seed.ops;
else
    grid = result.config.grid;
    grid.nx = numel(x);
    grid.ny = numel(y);
    grid.customX = x;
    grid.customY = y;
    ops = ipm.mesh.build(result.config,grid);
    ops.baseX = seed.ops.baseX;
    ops.baseY = seed.ops.baseY;
    ops.rescaling = ipm.output.restoreRuntimeReferences( ...
        ops.rescaling,seed.ops.rescaling,x, ...
        strcmp(result.config.scaling.rescalingMode,'dynamic'));
end
ops.remeshCount = result.grid.remeshCount;
if ~isequal(ops.x,x) || ~isequal(ops.y,y)
    error('ipm:CheckpointResultGrid', ...
        'Terminal operator reconstruction changed the stored result grid.');
end
end

function scale = rebuild_scale(result)
geometry = char(string(result.scale.geometry));
if ~strcmp(geometry,result.config.scaling.dynamicScaleGeometry)
    error('ipm:CheckpointResultScale', ...
        'Result scale geometry disagrees with its frozen configuration.');
end
if strcmp(geometry,'anisotropic')
    scale = struct('logC_l',log(result.scale.Cx), ...
        'logC_x',log(result.scale.Cx),'logC_y',log(result.scale.Cy), ...
        'logC_omega',log(result.scale.Comega), ...
        'physicalTime',result.state.physicalTime, ...
        'X_shift',result.scale.Xshift, ...
        'canonicalTime',result.state.canonicalTime);
else
    assert_numeric_overlap(result.scale.Cx,result.scale.Cy, ...
        'isotropic Cx/Cy',128*eps);
    scale = struct('logC_l',log(result.scale.Cx), ...
        'logC_omega',log(result.scale.Comega), ...
        'physicalTime',result.state.physicalTime, ...
        'X_shift',result.scale.Xshift, ...
        'canonicalTime',result.state.canonicalTime);
end
assert_numeric_overlap(exp(scale.logC_l),result.scale.Cx, ...
    'terminal Cx',128*eps);
if strcmp(geometry,'anisotropic')
    assert_numeric_overlap(exp(scale.logC_y),result.scale.Cy, ...
        'terminal Cy',128*eps);
end
assert_numeric_overlap(exp(scale.logC_omega),result.scale.Comega, ...
    'terminal Comega',128*eps);
end

function validate_flow_overlap(flow,resultState)
assert_numeric_overlap(flow.psi,resultState.psi,'terminal psi',2e-10);
assert_numeric_overlap(flow.source,resultState.omega, ...
    'terminal omega',2e-10);
assert_numeric_overlap(flow.u1,resultState.velocity.x, ...
    'terminal velocity x',2e-10);
assert_numeric_overlap(flow.u2,resultState.velocity.y, ...
    'terminal velocity y',2e-10);
end

function validate_terminal_overlap(state,result)
probeState = state;
probeState.config.output.verbose = false;
probeState.config.output.storeSnapshots = false;
[probeLog,~] = ipm.output.record(ipm.output.initializeLog(),probeState);
terminalIndex = numel(result.history.common.t);
compare_history_row( ...
    probeLog.history,result.history,terminalIndex,'terminal');
end

function compare_history_row(actualHistory,expectedHistory,index,context)
groups = {'common','gauge','mesh','anisotropic'};
for groupIndex = 1:numel(groups)
    group = groups{groupIndex};
    actual = actualHistory.(group);
    expected = expectedHistory.(group);
    names = fieldnames(actual);
    for nameIndex = 1:numel(names)
        name = names{nameIndex};
        if ~isfield(expected,name) || numel(expected.(name)) < index
            if is_optional_reconstructed_field(group,name)
                continue;
            end
            error('ipm:CheckpointResultOverlap', ...
                '%s history is missing reconstructed field %s.%s.', ...
                context,group,name);
        end
        actualValue = actual.(name)(1);
        expectedValue = expected.(name)(index);
        if isnumeric(actualValue)
            assert_numeric_overlap(actualValue,expectedValue, ...
                sprintf('%s history %s.%s',context,group,name),2e-10);
        elseif ~isequaln(actualValue,expectedValue)
            error('ipm:CheckpointResultOverlap', ...
                '%s history field %s.%s did not reproduce.', ...
                context,group,name);
        end
    end
end
end

function timeStep = terminal_timestep(result,scale)
common = result.history.common;
mapping = { ...
    'dt','acceptedCanonicalDt'; ...
    'rate','transportRate'; ...
    'dtCfl','cflDtLimit'; ...
    'maxDtLimit','maximumDtLimit'; ...
    'canonicalLimit','canonicalEndpointDtLimit'; ...
    'physicalLimit','physicalEndpointDtLimit'; ...
    'physicalClockSpeed','physicalClockSpeed'; ...
    'realizedCfl','realizedCfl'};
historyNames = mapping(:,2);
if all(isfield(common,historyNames)) && ...
        isfield(common,'timestepActiveLimiter')
    timeStep = struct();
    last = numel(common.t);
    for index = 1:size(mapping,1)
        timeStep.(mapping{index,1}) = common.(mapping{index,2})(last);
    end
    timeStep.activeLimiter = char(string( ...
        common.timestepActiveLimiter(last)));
else
    timeStep = ipm.evolve.emptyTimestep( ...
        result.config.time,scale,'unavailable');
end
end

function optional = is_optional_reconstructed_field(group,name)
common = {'acceptedStep','acceptedCanonicalDt','transportRate','cflDtLimit', ...
    'maximumDtLimit','canonicalEndpointDtLimit', ...
    'physicalEndpointDtLimit','physicalClockSpeed','realizedCfl', ...
    'timestepActiveLimiter','amplitudeSourceStep', ...
    'conservativeSourceStep','physicalDomainXMinimum', ...
    'physicalDomainXMaximum','physicalDomainXRadius', ...
    'physicalDomainYMaximum','physicalMinimumDx','physicalMinimumDy'};
gauge = {'omegaGaugeEnergyX','omegaGaugeEnergyY', ...
    'omegaGaugeForcingX','omegaGaugeForcingY','omegaGaugeForcing', ...
    'omegaGaugeAbsoluteForcing','omegaGaugeForcingCancellationRatio', ...
    'cOmegaOverCL','cLMinusCOmegaOverCL'};
optional = (strcmp(group,'common') && any(strcmp(name,common))) || ...
    (strcmp(group,'gauge') && any(strcmp(name,gauge)));
end

function assert_numeric_overlap(actual,expected,label,relativeTolerance)
if ~isnumeric(actual) || ~isnumeric(expected) || ...
        ~isequal(size(actual),size(expected)) || ...
        any(~isfinite(actual),'all') || any(~isfinite(expected),'all')
    if isequaln(actual,expected)
        return;
    end
    error('ipm:CheckpointResultOverlap', ...
        'Reconstructed %s is not comparable to the terminal result.',label);
end
scale = max([1,max(abs(actual),[],'all'),max(abs(expected),[],'all')]);
difference = max(abs(double(actual)-double(expected)),[],'all');
if difference > relativeTolerance*scale
    error('ipm:CheckpointResultOverlap', ...
        ['Reconstructed %s differs from the terminal result by %.3e ' ...
        '(allowed %.3e).'],label,difference,relativeTolerance*scale);
end
end

function log = result_log(result)
log = struct('history',result.history, ...
    'snapshotRho',{result.snapshots.rho}, ...
    'snapshotNormalizedTime',result.snapshots.normalizedTime, ...
    'snapshotX',{result.snapshots.x},'snapshotY',{result.snapshots.y});
end

function cursor = result_cursor(result)
nextOutput = result.config.time.outputEvery;
common = result.history.common;
for index = 2:numel(common.canonicalTau)
    tau = common.canonicalTau(index);
    due = tau+10*eps(tau) >= nextOutput || ...
        tau >= result.config.time.finalTime || ...
        common.physicalTime(index) >= ...
        result.config.time.physicalFinalTime;
    if due
        nextOutput = nextOutput+result.config.time.outputEvery;
    elseif index < numel(common.canonicalTau)
        error('ipm:CheckpointResultCursor', ...
            'A nonterminal result record cannot be placed on the output clock.');
    end
end

nextCheckpoint = Inf;
if isfield(result.config.output,'checkpoint') && ...
        result.config.output.checkpoint.enabled
    every = result.config.output.checkpoint.every;
    if isfinite(every)
        nextCheckpoint = every;
        tau = result.state.canonicalTime;
        tolerance = 10*eps(max(1,abs(tau)));
        while tau+tolerance >= nextCheckpoint
            nextCheckpoint = nextCheckpoint+every;
        end
    end
end
cursor = struct('nextOutput',nextOutput, ...
    'nextCheckpoint',nextCheckpoint, ...
    'lastCheckpointStep',result.state.steps, ...
    'lastCheckpointCanonicalTime',result.state.canonicalTime);
end
