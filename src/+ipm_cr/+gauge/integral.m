function integral = integral(wallGradient,ops,settings)
%IPM.EVOLVE.WALLGRADIENTINTEGRAL Integrate R_X on the wall over [-L,L].
%   L=wallIntegralRadius is fixed in rescaled X coordinates. The weights
%   integrate the piecewise-linear wall-gradient interpolant exactly, including
%   partial end cells when -L or L is not a mesh node. Weights are rebuilt
%   after remeshing; the coordinate interval and initial reference are fixed.

x = reshape(ops.x,1,[]);
radius = settings.wallIntegralRadius;

left = x(1:end-1);
right = x(2:end);
lower = max(left,-radius);
upper = min(right,radius);
active = upper > lower;
indices = find(active);
cellWidth = right(active)-left(active);
clippedWidth = upper(active)-lower(active);
averageCoordinate = ((lower(active)-left(active)) + ...
    (upper(active)-left(active)))./(2*cellWidth);
weights = zeros(size(x));
weights(indices) = clippedWidth.*(1-averageCoordinate);
weights(indices+1) = weights(indices+1) + ...
    clippedWidth.*averageCoordinate;
support = weights > 0;
supportPoints = nnz(support);

wallGradient = reshape(wallGradient,1,[]);
value = sum(weights.*wallGradient);
absoluteValue = sum(weights.*abs(wallGradient));
condition = value/max(absoluteValue,realmin);
integral = struct('weights',weights,'value',value, ...
    'radius',radius,'center',0,'support',support, ...
    'supportPoints',supportPoints,'condition',condition, ...
    'absoluteValue',absoluteValue,'meanValue',value/(2*radius));
end
