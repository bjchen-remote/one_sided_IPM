function [newValues,info] = interpolate( ...
        oldAxis,oldValues,newAxis,options)
%IPM.REMESH.INTERPOLATE Local polynomial interpolation on monotone axes.
%   VNEW = IPM.REMESH.INTERPOLATE(XOLD,VOLD,XNEW) evaluates the degree-five
%   interpolant through a local six-point stencil. OPTIONS.stencilWidth=8
%   selects the degree-seven interpolant used only by the sixth-order path.
%   Interior and boundary targets use the same stencil width, so there is
%   no low-order boundary closure. XOLD and XNEW may be increasing or
%   decreasing and need not be uniform.
%
%   For a matrix VOLD, samples occupy its first dimension by default.  Set
%   OPTIONS.sampleDimension to 2 to interpolate row signals.  If only one
%   matrix dimension matches NUMEL(XOLD), that dimension is inferred.
%
%   OPTIONS.conservation may be 'none' (the default) or 'constant'.  The
%   latter adds one constant to each signal so that its old- and new-grid
%   integrals agree.  Width six uses composite local degree-five quadrature;
%   width eight uses the mapped degree-seven sixth-order norm.  In either
%   case the correction retains the interpolation's formal smooth-grid
%   order. A constant correction is deliberately global and smooth, but it
%   is not range preserving.
%   Conservation requires equal old/new domain endpoints and at least one
%   complete interpolation stencil on the new grid.

if nargin < 4 || isempty(options)
    options = struct();
end
[conservation,sampleDimension,stencilWidth] = parse_options(options);

validateattributes(oldAxis,{'numeric'}, ...
    {'vector','real','finite','nonempty'},mfilename,'oldAxis');
validateattributes(newAxis,{'numeric'}, ...
    {'vector','real','finite','nonempty'},mfilename,'newAxis');
minimumGridNodes = 6;
if stencilWidth == 8
    minimumGridNodes = 9;
end
if numel(oldAxis) < minimumGridNodes
    error('ipm:HighOrderInterpolationGrid', ...
        ['High-order interpolation requires at least as many old-axis ' ...
        'nodes as its selected stencil.']);
end
oldDirection = monotone_direction(oldAxis,'oldAxis');
newDirection = monotone_direction(newAxis,'newAxis');

validateattributes(oldValues,{'numeric'}, ...
    {'2d','nonempty','finite'},mfilename,'oldValues');
numberOfOldNodes = numel(oldAxis);
[signals,isVectorValue,sampleDimension] = orient_signals( ...
    oldValues,numberOfOldNodes,sampleDimension);
signals = double(signals);

axis = double(oldAxis(:));
if oldDirection < 0
    axis = flipud(axis);
    signals = flipud(signals);
end
queries = double(newAxis(:));
if stencilWidth == 6
    rawValues = interpolate_local_degree_five(axis,signals,queries);
else
    rawValues = interpolate_local_degree_seven(axis,signals,queries);
end

numberOfSignals = size(signals,2);
correction = zeros(1,numberOfSignals);
massOld = nan(1,numberOfSignals);
massBefore = nan(1,numberOfSignals);
massAfter = nan(1,numberOfSignals);
sameDomain = domains_match(axis,queries);
canMeasureMass = sameDomain && numel(queries) >= minimumGridNodes;
needsMass = strcmp(conservation,'constant') || nargout > 1;
if needsMass && canMeasureMass
    [queryAscending,rawAscending] = ascending_queries( ...
        queries,rawValues,newDirection);
    if stencilWidth == 8
        oldWeights = ipm.mesh.quadrature6(axis);
        newWeights = ipm.mesh.quadrature6(queryAscending);
    else
        oldWeights = ipm.mesh.quadrature(axis);
        newWeights = ipm.mesh.quadrature(queryAscending);
    end
    massOld = oldWeights.'*signals;
    massBefore = newWeights.'*rawAscending;
end

correctedValues = rawValues;
if strcmp(conservation,'constant')
    if ~canMeasureMass
        error('ipm:HighOrderConservationDomain', ...
            ['Conservative interpolation requires matching domain endpoints ' ...
            'and at least one complete new-axis stencil.']);
    end
    domainLength = axis(end)-axis(1);
    correction = (massOld-massBefore)/domainLength;
    correctedValues = correctedValues+correction;

    % Remove the last matrix-product roundoff without changing the formal
    % correction.  Usually this second update is exactly zero.
    if newDirection > 0
        correctedAscending = correctedValues;
    else
        correctedAscending = flipud(correctedValues);
    end
    residual = massOld-newWeights.'*correctedAscending;
    correction = correction+residual/domainLength;
    correctedValues = correctedValues+residual/domainLength;
end
if nargout > 1 && canMeasureMass
    if newDirection > 0
        correctedAscending = correctedValues;
    else
        correctedAscending = flipud(correctedValues);
    end
    massAfter = newWeights.'*correctedAscending;
end

newValues = restore_shape(correctedValues,isVectorValue, ...
    sampleDimension,size(newAxis));
if nargout > 1
    info = struct( ...
        'stencilWidth',stencilWidth, ...
        'polynomialDegree',stencilWidth-1, ...
        'sampleDimension',sampleDimension, ...
        'conservation',conservation, ...
        'sameDomain',sameDomain, ...
        'extrapolatedNodeCount', ...
            sum(queries < axis(1) | queries > axis(end)), ...
        'massOld',massOld, ...
        'massBeforeCorrection',massBefore, ...
        'massAfterCorrection',massAfter, ...
        'massDefectBefore',massBefore-massOld, ...
        'massDefectAfter',massAfter-massOld, ...
        'constantCorrection',correction, ...
        'rangeViolationBefore',range_violation(signals,rawValues), ...
        'rangeViolationAfter',range_violation(signals,correctedValues));
end
end

function [conservation,sampleDimension,stencilWidth] = parse_options(options)
sampleDimension = [];
stencilWidth = 6;
if ischar(options) || (isstring(options) && isscalar(options))
    conservation = char(options);
elseif isstruct(options) && isscalar(options)
    if isfield(options,'conservation')
        conservation = options.conservation;
    else
        conservation = 'none';
    end
    if isfield(options,'sampleDimension')
        sampleDimension = options.sampleDimension;
        validateattributes(sampleDimension,{'numeric'}, ...
            {'scalar','integer','>=',1,'<=',2},mfilename, ...
            'options.sampleDimension');
    end
    if isfield(options,'stencilWidth')
        stencilWidth = options.stencilWidth;
    end
else
    error('ipm:HighOrderInterpolationOptions', ...
        'Options must be a scalar structure or a conservation-mode string.');
end
if ~(ischar(conservation) || ...
        (isstring(conservation) && isscalar(conservation)))
    error('ipm:HighOrderConservationMode', ...
        'The conservation mode must be ''none'' or ''constant''.');
end
conservation = lower(strrep(char(conservation),'-','_'));
if strcmp(conservation,'global_constant')
    conservation = 'constant';
end
if ~ismember(conservation,{'none','constant'})
    error('ipm:HighOrderConservationMode', ...
        'The conservation mode must be ''none'' or ''constant''.');
end
validateattributes(stencilWidth,{'numeric'}, ...
    {'scalar','integer'},mfilename,'options.stencilWidth');
if ~ismember(stencilWidth,[6,8])
    error('ipm:HighOrderInterpolationStencil', ...
        'OPTIONS.stencilWidth must be 6 or 8.');
end
end

function direction = monotone_direction(axis,name)
if isscalar(axis)
    direction = 1;
    return;
end
increments = diff(axis(:));
if all(increments > 0)
    direction = 1;
elseif all(increments < 0)
    direction = -1;
else
    error('ipm:HighOrderInterpolationAxis', ...
        '%s must be strictly monotone.',name);
end
end

function [signals,isVectorValue,sampleDimension] = orient_signals( ...
        oldValues,numberOfOldNodes,sampleDimension)
isVectorValue = isvector(oldValues) && numel(oldValues) == numberOfOldNodes;
if isVectorValue
    signals = oldValues(:);
    if isempty(sampleDimension)
        if isrow(oldValues)
            sampleDimension = 2;
        else
            sampleDimension = 1;
        end
    end
    return;
end
matchingDimensions = [size(oldValues,1) == numberOfOldNodes, ...
    size(oldValues,2) == numberOfOldNodes];
if isempty(sampleDimension)
    if matchingDimensions(1)
        sampleDimension = 1;
    elseif matchingDimensions(2)
        sampleDimension = 2;
    else
        error('ipm:HighOrderInterpolationValues', ...
            'One oldValues dimension must equal numel(oldAxis).');
    end
elseif ~matchingDimensions(sampleDimension)
    error('ipm:HighOrderInterpolationValues', ...
        'The selected sample dimension must equal numel(oldAxis).');
end
if sampleDimension == 1
    signals = oldValues;
else
    signals = oldValues.';
end
end

function values = interpolate_local_degree_five(axis,signals,queries)
numberOfOldNodes = numel(axis);
numberOfStencils = numberOfOldNodes-5;
centers = zeros(1,numberOfStencils);
scales = zeros(1,numberOfStencils);
normalizedNodes = zeros(6,numberOfStencils);
barycentricWeights = zeros(6,numberOfStencils);
for first = 1:numberOfStencils
    stencil = first:first+5;
    centers(first) = 0.5*(axis(stencil(1))+axis(stencil(end)));
    scales(first) = 0.5*(axis(stencil(end))-axis(stencil(1)));
    nodes = (axis(stencil)-centers(first))/scales(first);
    normalizedNodes(:,first) = nodes;
    for localNode = 1:6
        others = [1:localNode-1,localNode+1:6];
        barycentricWeights(localNode,first) = ...
            1/prod(nodes(localNode)-nodes(others));
    end
end

values = zeros(numel(queries),size(signals,2));
for queryIndex = 1:numel(queries)
    query = queries(queryIndex);
    if query <= axis(1)
        cellIndex = 1;
    elseif query >= axis(end)
        cellIndex = numberOfOldNodes-1;
    else
        cellIndex = find(axis <= query,1,'last');
    end
    first = min(max(cellIndex-2,1),numberOfStencils);
    stencil = first:first+5;
    nodes = normalizedNodes(:,first);
    normalizedQuery = (query-centers(first))/scales(first);
    [nearestDistance,nearestNode] = min(abs(nodes-normalizedQuery));
    if nearestDistance <= 8*eps(max(1,abs(normalizedQuery)))
        values(queryIndex,:) = signals(stencil(nearestNode),:);
    else
        factors = barycentricWeights(:,first)./(normalizedQuery-nodes);
        values(queryIndex,:) = ...
            (factors.'*signals(stencil,:))/sum(factors);
    end
end
end

function values = interpolate_local_degree_seven(axis,signals,queries)
numberOfOldNodes = numel(axis);
numberOfStencils = numberOfOldNodes-7;
centers = zeros(1,numberOfStencils);
scales = zeros(1,numberOfStencils);
normalizedNodes = zeros(8,numberOfStencils);
barycentricWeights = zeros(8,numberOfStencils);
for first = 1:numberOfStencils
    stencil = first:first+7;
    centers(first) = 0.5*(axis(stencil(1))+axis(stencil(end)));
    scales(first) = 0.5*(axis(stencil(end))-axis(stencil(1)));
    nodes = (axis(stencil)-centers(first))/scales(first);
    normalizedNodes(:,first) = nodes;
    for localNode = 1:8
        others = [1:localNode-1,localNode+1:8];
        barycentricWeights(localNode,first) = ...
            1/prod(nodes(localNode)-nodes(others));
    end
end

values = zeros(numel(queries),size(signals,2));
for queryIndex = 1:numel(queries)
    query = queries(queryIndex);
    if query <= axis(1)
        cellIndex = 1;
    elseif query >= axis(end)
        cellIndex = numberOfOldNodes-1;
    else
        cellIndex = find(axis <= query,1,'last');
    end
    first = min(max(cellIndex-3,1),numberOfStencils);
    stencil = first:first+7;
    nodes = normalizedNodes(:,first);
    normalizedQuery = (query-centers(first))/scales(first);
    [nearestDistance,nearestNode] = min(abs(nodes-normalizedQuery));
    if nearestDistance <= 8*eps(max(1,abs(normalizedQuery)))
        values(queryIndex,:) = signals(stencil(nearestNode),:);
    else
        factors = barycentricWeights(:,first)./(normalizedQuery-nodes);
        values(queryIndex,:) = ...
            (factors.'*signals(stencil,:))/sum(factors);
    end
end
end

function answer = domains_match(oldAscending,queries)
queryBounds = [min(queries),max(queries)];
oldBounds = oldAscending([1,end]).';
scale = max([1,abs(oldBounds),abs(queryBounds)]);
answer = all(abs(queryBounds-oldBounds) <= 128*eps(scale));
end

function [queryAscending,valuesAscending] = ascending_queries( ...
        queries,values,newDirection)
if newDirection > 0
    queryAscending = queries;
    valuesAscending = values;
else
    queryAscending = flipud(queries);
    valuesAscending = flipud(values);
end
end

function violation = range_violation(oldValues,newValues)
if ~isreal(oldValues) || ~isreal(newValues)
    violation = nan(1,size(oldValues,2));
    return;
end
oldMinimum = min(oldValues,[],1);
oldMaximum = max(oldValues,[],1);
newMinimum = min(newValues,[],1);
newMaximum = max(newValues,[],1);
violation = max(max(oldMinimum-newMinimum,newMaximum-oldMaximum),0);
end

function values = restore_shape(orientedValues,isVectorValue, ...
        sampleDimension,newAxisSize)
if isVectorValue
    values = reshape(orientedValues,newAxisSize);
elseif sampleDimension == 1
    values = orientedValues;
else
    values = orientedValues.';
end
end
