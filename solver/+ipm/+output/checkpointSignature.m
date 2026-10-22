function signature = checkpointSignature(payload)
%IPM.OUTPUT.CHECKPOINTSIGNATURE Compact integrity signature for a checkpoint.
% Keep its floating-point sums independent of the solver's thread count.
% This only affects checkpoint integrity metadata, not the PDE calculation.
previousThreads = maxNumCompThreads(1);
threadCleanup = onCleanup(@() maxNumCompThreads(previousThreads)); %#ok<NASGU>

state = payload.state;
signature = struct( ...
    'config',canonical_value(state.config), ...
    'runMetadata',canonical_value(state.runMetadata), ...
    'rho',numeric_signature(state.rho), ...
    'x',state.x,'y',state.y, ...
    'baseX',state.baseX,'baseY',state.baseY, ...
    'rescaling',canonical_value(state.rescaling), ...
    'remeshCount',state.remeshCount, ...
    'scale',canonical_value(state.scale), ...
    'normalizedTime',state.normalizedTime, ...
    'step',state.step,'mass0',state.mass0, ...
    'rhoRange0',state.rhoRange0, ...
    'history',canonical_value(payload.log.history), ...
    'snapshotSummary',snapshot_signature(payload.log), ...
    'cursor',canonical_value(payload.cursor));
end

function output = canonical_value(input)
if isa(input,'function_handle')
    output = struct('kind','function_handle','text',func2str(input));
elseif isstruct(input)
    if ~isscalar(input)
        output = arrayfun(@canonical_value,input,'UniformOutput',false);
        return;
    end
    output = struct();
    names = sort(fieldnames(input));
    for index = 1:numel(names)
        output.(names{index}) = canonical_value(input.(names{index}));
    end
elseif iscell(input)
    output = cell(size(input));
    for index = 1:numel(input)
        output{index} = canonical_value(input{index});
    end
else
    output = input;
end
end

function signature = numeric_signature(values)
flat = double(values(:));
indices = (1:numel(flat))';
firstWeights = mod(indices,104729)+1;
secondWeights = mod(65537*indices,130363)+1;
signature = struct('class',class(values),'size',size(values), ...
    'minimum',min(flat),'maximum',max(flat),'sum',sum(flat), ...
    'frobeniusNorm',norm(flat), ...
    'weightedSum1',sum(flat.*firstWeights), ...
    'weightedSum2',sum(flat.*secondWeights));
end

function signature = snapshot_signature(log)
count = numel(log.snapshotRho);
signature = struct('count',count, ...
    'normalizedTime',log.snapshotNormalizedTime);
if count == 0
    signature.rho = struct([]);
    signature.x = struct([]);
    signature.y = struct([]);
    return;
end
rho = repmat(numeric_signature(log.snapshotRho{1}),count,1);
x = repmat(numeric_signature(log.snapshotX{1}),count,1);
y = repmat(numeric_signature(log.snapshotY{1}),count,1);
for index = 1:count
    rho(index) = numeric_signature(log.snapshotRho{index});
    x(index) = numeric_signature(log.snapshotX{index});
    y(index) = numeric_signature(log.snapshotY{index});
end
signature.rho = rho;
signature.x = x;
signature.y = y;
end
