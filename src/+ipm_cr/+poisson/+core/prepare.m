function map=prepare(c,x,y)
%PREPARE Geometry-only evaluation matrices; no solution values are retained.
assert(isequal(size(x),size(y)),'ipm_cr:core:QuerySize', ...
    'x and y must have matching sizes.');
map=struct('x',x(:),'y',y(:),'xNodes',c.x,'yNodes',c.y, ...
    'Bx',ipm_cr.internal.bary(c.x,x),'By',ipm_cr.internal.bary(c.y,y));
end
