function fluxDerivative = weno7FluxDerivative( ...
    q,a,J,computationalSpacing,options)
%IPM.FIELD.WENO7FLUXDERIVATIVE Mapped nodal WENO-Z7 flux derivative.
%   DFDX = IPM.FIELD.WENO7FLUXDERIVATIVE(Q,A,J,DXI) approximates
%
%       d(A.*Q)/dx = (1/J) d(A.*Q)/dxi,   J = dx/dxi,
%
%   on a uniform computational grid.  Q stores nodes in columns.  The
%   global Lax--Friedrichs split acts on the mapped state J.*Q, and the
%   interface flux uses four degree-three substencils with the standard
%   WENO-Z7 global indicator tau7=abs(beta0-beta3).  The default nonlinear
%   weight power is q=4, the sufficient choice for seventh-order recovery
%   when both the first and second derivatives vanish.  A caller may set
%   OPTIONS.weightPower to any finite value >=2 for controlled studies.
%
%   A and J may be scalar, the size of Q, one value per line, or one row
%   shared by every line.  OPTIONS accepts lowerBoundary/upperBoundary
%   ('extrapolate' or 'reflect'), extrapolationDegree (6 by default, or 7),
%   epsilon, weightPower, and an optional LF alpha bound.  Four ghost nodes
%   are built on each side.  The scale-aware epsilon is
%   max(epsilon,1/(nNodes-1)^2).

if nargin < 5 || isempty(options)
    options = struct();
end

validateattributes(q,{'single','double'}, ...
    {'2d','real','finite','nonempty'},mfilename,'q',1);
if size(q,2) == 1
    error('ipm:Weno7FdOrientation', ...
        ['Q must store grid nodes in columns. Transpose a column vector ' ...
        'before calling this function.']);
end
validateattributes(computationalSpacing,{'single','double'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'computationalSpacing',4);
if ~isstruct(options) || ~isscalar(options)
    error('ipm:Weno7FdOptions','OPTIONS must be a scalar structure.');
end

[numberOfLines,numberOfNodes] = size(q);
a = expand_to_size(a,numberOfLines,numberOfNodes,'a',q);
J = expand_to_size(J,numberOfLines,numberOfNodes,'J',q);
if any(J(:) <= 0)
    error('ipm:Weno7FdMetric', ...
        'J must be finite and strictly positive at every node.');
end

lowerBoundary = boundary_option(options,'lowerBoundary','extrapolate');
upperBoundary = boundary_option(options,'upperBoundary','extrapolate');
extrapolationDegree = scalar_option(options,'extrapolationDegree',6);
if ~ismember(extrapolationDegree,[6,7])
    error('ipm:Weno7FdExtrapolationDegree', ...
        'OPTIONS.extrapolationDegree must be 6 or 7.');
end
epsilonFloor = scalar_option(options,'epsilon',1e-12);
if ~isfinite(epsilonFloor) || epsilonFloor <= 0
    error('ipm:Weno7FdEpsilon', ...
        'OPTIONS.epsilon must be a finite positive scalar.');
end

needsExtrapolation = strcmp(lowerBoundary,'extrapolate') || ...
    strcmp(upperBoundary,'extrapolate');
minimumNodes = 7;
if needsExtrapolation
    minimumNodes = max(minimumNodes,extrapolationDegree+1);
end
if numberOfNodes < minimumNodes
    error('ipm:Weno7FdGridSize', ...
        ['At least %d nodes are required for the selected WENO7 ' ...
        'boundary extension.'],minimumNodes);
end
normalizedEpsilon = max(epsilonFloor,1/(numberOfNodes-1)^2);
weightPower = scalar_option(options,'weightPower',4);
if ~isfinite(weightPower) || weightPower < 2
    error('ipm:Weno7FdWeightPower', ...
        'OPTIONS.weightPower must be finite and at least 2.');
end
if needsExtrapolation
    [lowerExtrapolationWeights,upperExtrapolationWeights] = ...
        extrapolation_weights(extrapolationDegree);
else
    lowerExtrapolationWeights = [];
    upperExtrapolationWeights = [];
end

minimumAlpha = max(abs(a./J),[],2);
alpha = alpha_option(options,minimumAlpha,numberOfLines,q);

extendedQ = extend_lines(q,lowerBoundary,upperBoundary, ...
    lowerExtrapolationWeights,upperExtrapolationWeights,1);
extendedA = extend_lines(a,lowerBoundary,upperBoundary, ...
    lowerExtrapolationWeights,upperExtrapolationWeights,-1);
extendedJ = extend_lines(J,lowerBoundary,upperBoundary, ...
    lowerExtrapolationWeights,upperExtrapolationWeights,1);

mappedState = extendedJ.*extendedQ;
physicalFlux = extendedA.*extendedQ;
positiveFlux = 0.5*(physicalFlux+alpha.*mappedState);
negativeFlux = 0.5*(physicalFlux-alpha.*mappedState);
physicalNodes = 5:numberOfNodes+4;
positiveScale = max(abs(positiveFlux(:,physicalNodes)),[],2);
negativeScale = max(abs(negativeFlux(:,physicalNodes)),[],2);
commonScale = max(positiveScale,negativeScale);
splitScaleFloor = sqrt(eps(class(q)))*commonScale;
positiveScale = max(positiveScale,splitScaleFloor);
negativeScale = max(negativeScale,splitScaleFloor);
[candidateWeights,smoothnessMatrices] = weno7_data();
linearWeights = cast([1,12,18,4]/35,'like',q);

faceFlux = zeros(numberOfLines,numberOfNodes+1,'like',q);
for face = 0:numberOfNodes
    positiveStencil = positiveFlux(:,face+(1:7));
    negativeStencil = negativeFlux(:,face+(8:-1:2));
    faceFlux(:,face+1) = weno_z7_left( ...
        positiveStencil,positiveScale,normalizedEpsilon,weightPower, ...
        candidateWeights,smoothnessMatrices,linearWeights) + ...
        weno_z7_left(negativeStencil,negativeScale,normalizedEpsilon, ...
        weightPower,candidateWeights,smoothnessMatrices,linearWeights);
end

computationalDerivative = diff(faceFlux,1,2)/computationalSpacing;
fluxDerivative = computationalDerivative./J;
end

function values = expand_to_size(values,numberOfLines,numberOfNodes, ...
        argumentName,prototype)
validateattributes(values,{'single','double'}, ...
    {'2d','real','finite','nonempty'},mfilename,argumentName);
if isscalar(values)
    values = repmat(values,numberOfLines,numberOfNodes);
elseif isequal(size(values),[numberOfLines,numberOfNodes])
    % Already in the required form.
elseif isequal(size(values),[numberOfLines,1])
    values = repmat(values,1,numberOfNodes);
elseif isvector(values) && numel(values) == numberOfNodes
    values = repmat(reshape(values,1,numberOfNodes),numberOfLines,1);
else
    error('ipm:Weno7FdSize', ...
        ['%s must be scalar, the size of Q, a 1-by-nNodes row, or an ' ...
        'nLines-by-1 column.'],upper(argumentName));
end
values = cast(values,'like',prototype);
end

function value = scalar_option(options,name,defaultValue)
if isfield(options,name)
    value = options.(name);
else
    value = defaultValue;
end
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('ipm:Weno7FdOptions', ...
        'OPTIONS.%s must be a real numeric scalar.',name);
end
value = double(value);
end

function value = boundary_option(options,name,defaultValue)
if isfield(options,name)
    value = options.(name);
else
    value = defaultValue;
end
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('ipm:Weno7FdBoundary', ...
        'OPTIONS.%s must be ''extrapolate'' or ''reflect''.',name);
end
value = lower(char(value));
if ~ismember(value,{'extrapolate','reflect'})
    error('ipm:Weno7FdBoundary', ...
        'OPTIONS.%s must be ''extrapolate'' or ''reflect''.',name);
end
end

function alpha = alpha_option(options,minimumAlpha,numberOfLines,prototype)
if ~isfield(options,'alpha') || isempty(options.alpha)
    alpha = minimumAlpha;
    return;
end
alpha = options.alpha;
validateattributes(alpha,{'single','double'}, ...
    {'vector','real','finite','nonnegative'},mfilename,'OPTIONS.alpha');
if isscalar(alpha)
    alpha = repmat(alpha,numberOfLines,1);
elseif numel(alpha) == numberOfLines
    alpha = reshape(alpha,numberOfLines,1);
else
    error('ipm:Weno7FdAlphaSize', ...
        'OPTIONS.alpha must be scalar or have one value per line.');
end
alpha = cast(alpha,'like',prototype);
tolerance = 100*eps(class(prototype)).*max(1,minimumAlpha);
if any(alpha < minimumAlpha-tolerance)
    error('ipm:Weno7FdAlphaBound', ...
        'OPTIONS.alpha must bound max(abs(A./J)) on every line.');
end
end

function extended = extend_lines(values,lowerBoundary,upperBoundary, ...
        lowerExtrapolationWeights,upperExtrapolationWeights, ...
        reflectionParity)
[numberOfLines,numberOfNodes] = size(values);
extended = zeros(numberOfLines,numberOfNodes+8,'like',values);
extended(:,5:numberOfNodes+4) = values;

if strcmp(lowerBoundary,'reflect')
    extended(:,1:4) = reflectionParity*values(:,[5,4,3,2]);
else
    stencilWidth = size(lowerExtrapolationWeights,2);
    extended(:,1:4) = values(:,1:stencilWidth)* ...
        lowerExtrapolationWeights.';
end

if strcmp(upperBoundary,'reflect')
    extended(:,numberOfNodes+(5:8)) = reflectionParity* ...
        values(:,numberOfNodes-[1,2,3,4]);
else
    stencilWidth = size(upperExtrapolationWeights,2);
    extended(:,numberOfNodes+(5:8)) = ...
        values(:,numberOfNodes-stencilWidth+1:numberOfNodes)* ...
        upperExtrapolationWeights.';
end
end

function [lowerWeights,upperWeights] = extrapolation_weights(degree)
persistent lowerCache upperCache
if isempty(lowerCache)
    lowerCache = cell(1,2);
    upperCache = cell(1,2);
end
cacheIndex = degree-5;
if isempty(lowerCache{cacheIndex})
    lowerCache{cacheIndex} = lagrange_weights(0:degree,-4:-1);
    upperCache{cacheIndex} = lagrange_weights(-degree:0,1:4);
end
lowerWeights = lowerCache{cacheIndex};
upperWeights = upperCache{cacheIndex};
end

function weights = lagrange_weights(nodes,targets)
numberOfNodes = numel(nodes);
numberOfTargets = numel(targets);
weights = ones(numberOfTargets,numberOfNodes);
for node = 1:numberOfNodes
    others = [1:node-1,node+1:numberOfNodes];
    weights(:,node) = prod((targets(:)-nodes(others)) ./ ...
        (nodes(node)-nodes(others)),2);
end
end

function reconstructed = weno_z7_left( ...
        values,lineScale,epsilon,weightPower,candidateWeights, ...
        smoothnessMatrices,linearWeights)
safeScale = max(lineScale,sqrt(realmin(class(values))));
values = values./safeScale;

numberOfLines = size(values,1);
candidates = zeros(numberOfLines,4,'like',values);
beta = zeros(numberOfLines,4,'like',values);
for stencilIndex = 1:4
    localValues = values(:,stencilIndex:stencilIndex+3);
    candidates(:,stencilIndex) = localValues * ...
        candidateWeights(stencilIndex,:).';
    matrix = smoothnessMatrices(:,:,stencilIndex);
    beta(:,stencilIndex) = sum((localValues*matrix).*localValues,2);
end
beta = max(beta,0);
tau7 = abs(beta(:,1)-beta(:,4));
alpha = linearWeights.*(1+(tau7./(beta+epsilon)).^weightPower);
reconstructed = safeScale.*sum(alpha.*candidates,2)./sum(alpha,2);
end

function [candidateWeights,smoothnessMatrices] = weno7_data()
persistent storedCandidateWeights storedSmoothnessMatrices
if isempty(storedCandidateWeights)
    % Finite-difference WENO reconstructs a numerical flux whose interface
    % difference is high-order, not the pointwise Lagrange value at the
    % interface.  These are the Balsara--Shu seventh-order candidates.
    referenceCandidateWeights = [ ...
        -1/4,13/12,-23/12,25/12; ...
         1/12,-5/12,13/12,1/4; ...
        -1/12,7/12,7/12,-1/12; ...
         1/4,13/12,-5/12,1/12];
    storedCandidateWeights = zeros(4,4);
    storedSmoothnessMatrices = zeros(4,4,4);
    derivativeGram = zeros(4);
    for derivativeOrder = 1:3
        for leftPower = derivativeOrder:3
            leftCoefficient = factorial(leftPower) / ...
                factorial(leftPower-derivativeOrder);
            for rightPower = derivativeOrder:3
                rightCoefficient = factorial(rightPower) / ...
                    factorial(rightPower-derivativeOrder);
                power = leftPower+rightPower-2*derivativeOrder;
                if mod(power,2) == 0
                    integral = 2*(0.5^(power+1))/(power+1);
                else
                    integral = 0;
                end
                derivativeGram(leftPower+1,rightPower+1) = ...
                    derivativeGram(leftPower+1,rightPower+1) + ...
                    leftCoefficient*rightCoefficient*integral;
            end
        end
    end
    for stencilIndex = 1:4
        cells = (-4+stencilIndex):(-1+stencilIndex);
        vandermonde = zeros(4);
        for cellIndex = 1:4
            for power = 0:3
                vandermonde(cellIndex,power+1) = ...
                    ((cells(cellIndex)+0.5)^(power+1) - ...
                    (cells(cellIndex)-0.5)^(power+1))/(power+1);
            end
        end
        storedCandidateWeights(stencilIndex,:) = ...
            (0.5.^(0:3))/vandermonde;
        matrix = vandermonde'\(derivativeGram/vandermonde);
        storedSmoothnessMatrices(:,:,stencilIndex) = ...
            0.5*(matrix+matrix.');
    end
    if max(abs(storedCandidateWeights-referenceCandidateWeights), ...
            [],'all') > 100*eps
        error('ipm:Weno7InternalCoefficients', ...
            'The generated WENO7 candidates failed their exact checksum.');
    end
end
candidateWeights = storedCandidateWeights;
smoothnessMatrices = storedSmoothnessMatrices;
end
