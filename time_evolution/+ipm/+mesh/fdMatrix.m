function D = fdMatrix(axis,derivativeOrder,stencilWidth)
%IPM.MESH.FDMATRIX Local-polynomial finite-difference matrix on any monotone axis.
%   D = IPM.MESH.FDMATRIX(X,M,WIDTH) differentiates nodal data M times.  Every
%   row is exact for polynomials through degree WIDTH-1, including the
%   one-sided boundary closures.  WIDTH=7 therefore gives sixth-order first
%   derivatives and at least fifth-order second derivatives on a smooth grid.

if nargin < 3
    stencilWidth = 7;
end
validateattributes(axis,{'numeric'}, ...
    {'vector','real','finite','increasing'},mfilename,'axis');
validateattributes(derivativeOrder,{'numeric'}, ...
    {'scalar','integer','nonnegative'},mfilename,'derivativeOrder');
validateattributes(stencilWidth,{'numeric'}, ...
    {'scalar','integer','>=',derivativeOrder+1},mfilename,'stencilWidth');

axis = axis(:);
numberOfNodes = numel(axis);
if stencilWidth > numberOfNodes
    error('ipm:FiniteDifferenceStencil', ...
        'The stencil width cannot exceed the number of axis nodes.');
end

rows = zeros(numberOfNodes*stencilWidth,1);
columns = zeros(numberOfNodes*stencilWidth,1);
values = zeros(numberOfNodes*stencilWidth,1);
halfWidth = floor(stencilWidth/2);
cursor = 0;
for node = 1:numberOfNodes
    first = min(max(node-halfWidth,1),numberOfNodes-stencilWidth+1);
    stencil = first:first+stencilWidth-1;
    offsets = axis(stencil)-axis(node);
    scale = max(abs(offsets));
    if scale == 0
        error('ipm:FiniteDifferenceScale','A finite-difference stencil collapsed.');
    end
    normalized = offsets/scale;
    moments = normalized'.^((0:stencilWidth-1)');
    target = zeros(stencilWidth,1);
    target(derivativeOrder+1) = ...
        factorial(derivativeOrder)/scale^derivativeOrder;
    coefficients = moments\target;
    slots = cursor+(1:stencilWidth);
    rows(slots) = node;
    columns(slots) = stencil;
    values(slots) = coefficients;
    cursor = cursor+stencilWidth;
end
D = sparse(rows,columns,values,numberOfNodes,numberOfNodes);
end
