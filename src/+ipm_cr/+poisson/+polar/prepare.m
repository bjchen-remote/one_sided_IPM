function map=prepare(p,x,y,requested)
%PREPARE Cache query ownership and interpolation, never sampled fields.
if nargin<4,requested={'psi','u1','u2','laplacian','dt','dtheta'};end
if ischar(requested),requested={requested};end
assert(isequal(size(x),size(y)),'ipm_cr:polar:QuerySize', ...
    'x and y must have matching sizes.');
shape=size(x);x=x(:);y=y(:);radius=hypot(x,y);
wallTolerance=128*eps(max(1,radius));
valid=isfinite(x)&isfinite(y)&y>=-wallTolerance& ...
    radius>=p.rInner-128*eps(p.rInner)&radius<=p.rOuter+128*eps(p.rOuter);
theta=min(pi,max(0,atan2(max(y,0),x)));
t=log(max(p.rInner,min(p.rOuter,radius)));
first=any(ismember(requested,{'u1','u2','dtheta'}));
second=ismember('laplacian',requested);
panels=cell(p.nPanels,1);
for j=1:p.nPanels
    if j==p.nPanels
        selected=valid&t>=p.logEdges(j)&t<=p.logEdges(j+1);
    else
        selected=valid&t>=p.logEdges(j)&t<p.logEdges(j+1);
    end
    indices=find(selected);
    Ba=ipm_cr.internal.bary(p.theta,theta(selected));
    entry=struct('indices',indices, ...
        'Br',ipm_cr.internal.bary(p.panelNodes{j},t(selected)),'Ba',Ba);
    if first,entry.angularFirst=Ba*p.Dtheta;end
    if second,entry.angularSecond=Ba*p.Dtheta2;end
    panels{j}=entry;
end
map=struct('x',x,'y',y,'shape',shape,'radius',radius,'theta',theta, ...
    'valid',valid,'panels',{panels},'thetaNodes',p.theta,'logEdges',p.logEdges, ...
    'radialDegree',p.degree);
end
