function snapshots = snapshotAt(inputResult,index)
%IPM.OUTPUT.SNAPSHOTAT Return stored solver snapshots through one contract.
%   S = IPM.OUTPUT.SNAPSHOTAT(RESULT,K) returns snapshot K from a version-2
%   result. K='all' returns every stored snapshot. With no K,
%   the latest stored snapshot is returned.  When no snapshots were stored,
%   the terminal state is exposed as a single snapshot.

result = ipm.output.validate(inputResult);
stored = result.snapshots;
count = numel(stored.rho);
if count == 0
    allSnapshots = terminal_snapshot(result);
    count = 1;
else
    allSnapshots = stored_snapshots(result);
end

if nargin < 2 || isempty(index)
    index = count;
elseif (ischar(index) || (isstring(index) && isscalar(index))) && ...
        strcmpi(string(index),'all')
    index = 1:count;
else
    validateattributes(index,{'numeric'}, ...
        {'vector','integer','positive','<=',count});
end
snapshots = allSnapshots(index);
end

function snapshots = stored_snapshots(result)
stored = result.snapshots;
count = numel(stored.rho);
common = result.history.common;
isPhysical = strcmp(result.config.scaling.rescalingMode,'physical');
Cx = history_series(common,{'C_x','C_l'},count, ...
    scale_fallback(isPhysical,count,1,result.scale.Cx));
Cy = history_series(common,{'C_y'},count, ...
    scale_fallback(isPhysical,count,1,result.scale.Cy));
Comega = history_series(common,{'C_omega'},count, ...
    scale_fallback(isPhysical,count,1,result.scale.Comega));
Xshift = history_series(common,{'X_shift'},count, ...
    scale_fallback(isPhysical,count,0,result.scale.Xshift));
if any(Cx <= 0) || any(Cy <= 0) || any(Comega <= 0)
    error('ipm:SnapshotScale', ...
        'Snapshot scale factors Cx, Cy, and Comega must be positive.');
end

template = struct('index',0,'rho',[],'physicalRho',[], ...
    'x',[],'y',[],'physicalX',[],'physicalY',[], ...
    'normalizedTime',NaN,'canonicalTime',NaN,'physicalTime',NaN);
snapshots = repmat(template,count,1);
for snapshotIndex = 1:count
    rho = stored.rho{snapshotIndex};
    x = stored.x{snapshotIndex};
    y = stored.y{snapshotIndex};
    snapshots(snapshotIndex) = struct( ...
        'index',snapshotIndex,'rho',rho, ...
        'physicalRho',rho/Comega(snapshotIndex), ...
        'x',x,'y',y, ...
        'physicalX',(x-Xshift(snapshotIndex))/Cx(snapshotIndex), ...
        'physicalY',y/Cy(snapshotIndex), ...
        'normalizedTime',stored.normalizedTime(snapshotIndex), ...
        'canonicalTime',stored.canonicalTime(snapshotIndex), ...
        'physicalTime',stored.physicalTime(snapshotIndex));
end
end

function values = scale_fallback(isPhysical,count,physicalValue,finalValue)
if isPhysical
    values = repmat(physicalValue,count,1);
elseif count == 1
    values = finalValue;
else
    values = [];
end
end

function snapshot = terminal_snapshot(result)
snapshot = struct( ...
    'index',1,'rho',result.state.rho, ...
    'physicalRho',result.physical.rho, ...
    'x',result.grid.x,'y',result.grid.y, ...
    'physicalX',result.grid.physicalX, ...
    'physicalY',result.grid.physicalY, ...
    'normalizedTime',result.state.normalizedTime, ...
    'canonicalTime',result.state.canonicalTime, ...
    'physicalTime',result.state.physicalTime);
end

function values = history_series(history,names,count,fallback)
values = [];
for index = 1:numel(names)
    if isfield(history,names{index}) && ~isempty(history.(names{index}))
        values = history.(names{index});
        break;
    end
end
if isempty(values)
    values = fallback;
end
if isempty(values)
    error('ipm:SnapshotScaleHistory', ...
        'Dynamic snapshots require a matching scale history.');
end
values = values(:);
if isscalar(values)
    values = repmat(values,count,1);
elseif numel(values) ~= count
    error('ipm:SnapshotHistoryCount', ...
        'Snapshot count does not match the stored scaling history.');
end
if any(~isfinite(values))
    error('ipm:SnapshotScale', ...
        'Snapshot scaling history must be finite.');
end
end
