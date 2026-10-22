function weights = quadrantQuadrature(x,stencilWidth)
%IPM.MESH.QUADRANTQUADRATURE Positive-half norm with even virtual stencils.
% Integrates only cells [0,H], using the same symmetric local polynomial
% closure as the full-axis norm. It stores only positive-axis weights.
if nargin<2,stencilWidth=6;end
validateattributes(x,{'numeric'},{'vector','real','finite','increasing'});
assert(x(1)==0 && numel(x)>=stencilWidth && ...
    any(stencilWidth==[6,8]),'ipm:QuadrantQuadratureAxis', ...
    'A positive axis and a six- or eight-node stencil are required.');
wasRow=isrow(x);x=x(:);n=numel(x);weights=zeros(n,1);
virtualCount=2*n-1;offset=stencilWidth/2-1;
for cellIndex=1:n-1
    virtualCell=n+cellIndex-1;
    first=min(max(virtualCell-offset,1),virtualCount-stencilWidth+1);
    stencil=first:first+stencilWidth-1;
    indices=abs(stencil-n)+1;
    nodes=x(indices);nodes(stencil<n)=-nodes(stencil<n);
    center=(x(cellIndex)+x(cellIndex+1))/2;
    scale=max(abs(nodes-center));
    normalized=(nodes-center)/scale;
    moments=normalized(:)'.^((0:stencilWidth-1)');
    lower=(x(cellIndex)-center)/scale;
    upper=(x(cellIndex+1)-center)/scale;
    integrals=scale*((upper.^(1:stencilWidth)- ...
        lower.^(1:stencilWidth))./(1:stencilWidth))';
    local=moments\integrals;
    weights=weights+accumarray(indices(:),local,[n,1],@sum,0);
end
if wasRow,weights=weights';end
end
