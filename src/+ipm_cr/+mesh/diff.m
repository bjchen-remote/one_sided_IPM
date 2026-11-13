function d=diff(rho,g)
%DIFFERENTIATE Local polynomial derivatives, annihilating constants first.
d=struct();
if isfield(g.core,'stair')
    d.core=ipm_cr.mesh.staircase('diff',rho.core,g.core);
else
    d.core=struct('x',apply(g.core.Dx,rho.core')','y',apply(g.core.Dy,rho.core));
end
fr=apply(g.polar.Dr,rho.polar);ft=apply(g.polar.Dt,rho.polar')';
d.polar=struct('x',fr.*g.polar.cosTheta-ft.*g.polar.sinTheta./g.polar.r, ...
    'y',fr.*g.polar.sinTheta+ft.*g.polar.cosTheta./g.polar.r);
end

function result=apply(D,field)
[rows,cols,weights]=find(D);
map=sparse(rows,(1:numel(rows))',weights,size(D,1),numel(rows));
result=map*(field(cols,:)-field(rows,:));
end
