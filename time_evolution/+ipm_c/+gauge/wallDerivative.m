function q=wallDerivative(values,Dx)
%WALLDERIVATIVE Apply the original wall derivative to local value differences.
% Each stencil annihilates constants before multiplication; large platforms
% cannot produce a derivative from the floating point row sum of Dx.
persistent previous rows columns weights
if isempty(previous)||~isequal(previous,Dx)
    [rows,columns,weights]=find(Dx);previous=Dx;
end
values=values(:);
assert(numel(values)==size(Dx,2)&&all(isfinite(values)), ...
    'polar:WallDerivative','The original wall data must be finite.');
q=accumarray(rows,weights.*(values(columns)-values(rows)), ...
    [size(Dx,1),1],@sum,0).';
end
