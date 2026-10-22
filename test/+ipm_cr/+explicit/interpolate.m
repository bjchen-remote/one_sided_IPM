function [values,map]=interpolate(fields,g,x,y)
%INTERPOLATE Local degree-six values with one active-core ownership rule.
% FIELDS.core is y-by-x; FIELDS.polar is r-by-theta, without a polar origin.
% Inside the active Cartesian rectangle use Cartesian data; elsewhere use
% the annulus. Full arrays, including synchronized halos, supply stencils.
% Optional MAP contains geometry-only sparse core/polar maps for reuse:
% reshape(map.core*fields.core(:)+map.polar*fields.polar(:),map.shape).
% Invalid queries return NaN. Only tiny wall/outer-radius roundoff is repaired.
assert(isstruct(g)&&all(isfield(g,{'core','polar'})), ...
    'ipm_cr:ExplicitInterpolationGrid','Supply the hybrid transport geometry.');
c=g.core;p=g.polar;xAxis=c.x(:)';yAxis=c.y(:);r=p.r(:);theta=p.theta(:)';
assert(isstruct(fields)&&isscalar(fields)&&all(isfield(fields,{'core','polar'}))&& ...
    isnumeric(fields.core)&&isreal(fields.core)&&isequal(size(fields.core),[numel(yAxis),numel(xAxis)])&& ...
    isnumeric(fields.polar)&&isreal(fields.polar)&&isequal(size(fields.polar),[numel(r),numel(theta)])&& ...
    all(isfinite(fields.core),'all')&&all(isfinite(fields.polar),'all'), ...
    'ipm_cr:ExplicitInterpolationField','Use complete finite native core/polar scalar fields.');
assert(isnumeric(x)&&isreal(x)&&isnumeric(y)&&isreal(y)&&isequal(size(x),size(y)), ...
    'ipm_cr:ExplicitInterpolationQueries','Use matching real query arrays.');
assert(all(isfield(c,{'activeRows','activeColumns'})), ...
    'ipm_cr:ExplicitInterpolationGrid','The geometry must identify the active Cartesian rectangle.');
columns=c.activeColumns;rows=c.activeRows;
bounds=[xAxis(columns(1)),xAxis(columns(end)),yAxis(rows(1)),yAxis(rows(end))];
key=struct('x',xAxis,'y',yAxis,'r',r,'theta',theta,'bounds',bounds);
persistent retained
hit=[];
for k=1:numel(retained)
    if isequaln(retained{k}.key,key)&&isequaln(retained{k}.x,x)&&isequaln(retained{k}.y,y)
        hit=k;break;
    end
end
if isempty(hit)
    map=prepare_map(key,x,y);
    entry=struct('key',key,'x',x,'y',y,'map',map);
    if isempty(retained),retained={entry};else,retained{end+1}=entry;end
    % Retain only geometry/query maps, never a sampled field or value.
    if numel(retained)>4,retained=retained(end-3:end);end
else
    map=retained{hit}.map;
end
values=map.core*fields.core(:)+map.polar*fields.polar(:);
values(~map.valid)=NaN;values=reshape(values,map.shape);
end

function map=prepare_map(key,x,y)
shape=size(x);x=x(:);y=y(:);n=numel(x);width=7;
assert(min([numel(key.x),numel(key.y),numel(key.r),numel(key.theta)])>=width&& ...
    all(isfinite([key.x(:);key.y(:);key.r(:);key.theta(:)]))&& ...
    all(diff(key.x)>0)&&all(diff(key.y)>0)&&all(diff(key.r)>0)&& ...
    all(diff(key.theta)>0)&&key.r(1)>0, ...
    'ipm_cr:ExplicitInterpolationAxes','Finite increasing axes with seven nodes and positive annular radii are required.');
wallTolerance=8*eps(max(1,abs(x)));
snap=isfinite(x)&isfinite(y)&y<0&y>=-wallTolerance;y(snap)=0;
radius=hypot(x,y);finite=isfinite(x)&isfinite(y)&y>=0;
withinDisk=radius<=key.r(end)*(1+16*eps);
coreOwner=finite&withinDisk&x>=key.bounds(1)&x<=key.bounds(2)& ...
    y>=key.bounds(3)&y<=key.bounds(4);
polarOwner=finite&withinDisk&~coreOwner&radius>=key.r(1)*(1-16*eps);
core=spalloc(n,numel(key.x)*numel(key.y),nnz(coreOwner)*width^2);
polar=spalloc(n,numel(key.r)*numel(key.theta),nnz(polarOwner)*width^2);
indices=find(coreOwner);
if ~isempty(indices)
    local=ipm_cr.carrier.tensorMap(key.y,key.x,y(indices),x(indices),width);
    [row,column,value]=find(local);core=sparse(indices(row),column,value,n,size(core,2));
end
indices=find(polarOwner);
if ~isempty(indices)
    rq=max(key.r(1),min(key.r(end),radius(indices)));
    tq=max(key.theta(1),min(key.theta(end),atan2(y(indices),x(indices))));
    local=ipm_cr.carrier.tensorMap(key.r,key.theta,rq,tq,width);
    [row,column,value]=find(local);polar=sparse(indices(row),column,value,n,size(polar,2));
end
valid=coreOwner|polarOwner;
rowNorm=full(sum(abs(core),2)+sum(abs(polar),2));
if isempty(rowNorm),operatorNorm=0;else,operatorNorm=max(rowNorm);end
map=struct('core',core,'polar',polar,'shape',shape,'valid',valid, ...
    'coreOwner',coreOwner,'polarOwner',polarOwner,'allValid',all(valid), ...
    'activeCoreBounds',key.bounds,'interpolationWidth',width, ...
    'operatorInfinityNorm',operatorNorm,'fieldCacheRetained',false);
end
