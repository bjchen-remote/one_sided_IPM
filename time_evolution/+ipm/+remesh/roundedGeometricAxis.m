function [axis,info]=roundedGeometricAxis(referenceAxis,coreWidth,targetCells,roundingCells)
%IPM.REMESH.ROUNDEDGEOMETRICAXIS Equalized positive axis with balanced log growth.
% Cell log-width grows quadratically across the first roundingCells cells,
% then linearly with a constant interface ratio across the outer region.
% A scalar bisection chooses the smallest growth rate meeting the requested
% fractional core-cell count at the same node count and box size.
v=referenceAxis(:);
validateattributes(v,{'numeric'},{'vector','real','finite','increasing'});
assert(numel(v)>=17&&v(1)==0&&v(end)>0,'ipm:RoundedGeometricReference');
validateattributes(coreWidth,{'numeric'},{'scalar','real','finite','positive','<',v(end)});
validateattributes(targetCells,{'numeric'},{'scalar','real','finite','positive','<',numel(v)-1});
validateattributes(roundingCells,{'numeric'},{'scalar','integer','>=',2,'<',numel(v)-1});
n=numel(v)-1;H=v(end);j=(0:n-1)';
shape=j-roundingCells/2;
inside=j<=roundingCells;
shape(inside)=j(inside).^2/(2*roundingCells);
make=@(s)axis_from_slope(s,shape,H);
lower=0;upper=log(1.08);
upperAxis=make(upper);
if interval_count(upperAxis,coreWidth)<targetCells
    error('ipm:RoundedGeometricInfeasible', ...
        'The requested positive-axis core count requires adjacent ratio above 1.08.');
end
for k=1:55
    middle=.5*(lower+upper);
    if interval_count(make(middle),coreWidth)<targetCells
        lower=middle;
    else
        upper=middle;
    end
end
% Give the fractional-count gate a small numerical margin without changing
% the design objective at a visible scale.
axis=make(upper+64*eps(upper));axis(1)=0;axis(end)=H;
width=diff(axis);
info=struct('kind','one_sided_rounded_geometric_equalizer_v1', ...
    'nodeCount',numel(v),'coreWidth',coreWidth,'targetCells',targetCells, ...
    'achievedCoreCells',interval_count(axis,coreWidth), ...
    'roundingCells',roundingCells,'logWidthGrowth',upper, ...
    'maximumAdjacentRatio',max(exp(abs(diff(log(width))))), ...
    'outerInterfaceCount',n-roundingCells-1, ...
    'sameBoxSameNodeCount',true);
end

function axis=axis_from_slope(s,shape,H)
u=s*(shape-max(shape));w=exp(u);w=H*w/sum(w);
axis=[0;cumsum(w)];axis(end)=H;
end

function count=interval_count(axis,coreWidth)
w=diff(axis);
count=sum(max(0,min(axis(2:end),coreWidth)-axis(1:end-1))./w);
end
