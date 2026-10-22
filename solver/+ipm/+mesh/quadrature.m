function weights = quadrature(axis,stencilWidth)
%IPM.MESH.QUADRATURE Composite local polynomial quadrature weights.
%   The default six-node local interpolant retains the certified fourth-order
%   path. STENCILWIDTH=8 integrates a local degree-seven interpolant for the
%   sixth-order path.

if nargin < 2
    stencilWidth = 6;
end

validateattributes(axis,{'numeric'}, ...
    {'vector','real','finite','increasing'},mfilename,'axis');
validateattributes(stencilWidth,{'numeric'}, ...
    {'scalar','integer'},mfilename,'stencilWidth');
if ~ismember(stencilWidth,[6,8])
    error('ipm:HighOrderQuadratureStencil', ...
        'The quadrature stencil width must be 6 or 8.');
end
if numel(axis) < stencilWidth
    error('ipm:HighOrderQuadratureGrid', ...
        ['High-order quadrature requires at least as many axis nodes as ' ...
        'its stencil width.']);
end
wasRow = isrow(axis);
axis = axis(:);
numberOfNodes = numel(axis);
weights = zeros(numberOfNodes,1);
for cellIndex = 1:numberOfNodes-1
    if stencilWidth == 6
        first = min(max(cellIndex-2,1),numberOfNodes-stencilWidth+1);
    else
        first = min(max(cellIndex-3,1),numberOfNodes-stencilWidth+1);
    end
    stencil = first:first+stencilWidth-1;
    center = (axis(cellIndex)+axis(cellIndex+1))/2;
    scale = max(abs(axis(stencil)-center));
    nodes = (axis(stencil)-center)/scale;
    lower = (axis(cellIndex)-center)/scale;
    upper = (axis(cellIndex+1)-center)/scale;
    moments = nodes'.^((0:stencilWidth-1)');
    integralMoments = scale*((upper.^(1:stencilWidth) - ...
        lower.^(1:stencilWidth))./(1:stencilWidth))';
    localWeights = moments\integralMoments;
    weights(stencil) = weights(stencil)+localWeights;
end
if wasRow
    weights = weights';
end
end
