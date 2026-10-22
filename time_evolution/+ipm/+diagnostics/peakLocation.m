function location = peakLocation(values,x,peakIndex)
%IPM.DIAGNOSTICS.PEAKLOCATION Quadratic subgrid location of a resolved line maximum.

values = values(:)';
x = x(:)';
if nargin < 3
    [~,peakIndex] = max(values);
end
location = x(peakIndex);
if peakIndex <= 1 || peakIndex >= numel(x)
    return;
end
indices = peakIndex-1:peakIndex+1;
localX = x(indices)-x(peakIndex);
coefficient = polyfit(localX,values(indices),2);
if coefficient(1) >= 0 || abs(coefficient(1)) <= eps || ...
        any(~isfinite(coefficient))
    return;
end
candidate = x(peakIndex)-coefficient(2)/(2*coefficient(1));
location = min(x(indices(end)),max(x(indices(1)),candidate));
end
