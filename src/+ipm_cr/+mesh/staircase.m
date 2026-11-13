function varargout=staircase(action,varargin)
%STAIRCASE Active Cartesian mask, native FEM support and active-only donors.
switch action
    case 'build',varargout={build(varargin{:})};
    case 'contains',varargout={contains_points(varargin{:})};
    case 'diff',varargout={differentiate(varargin{:})};
    case 'donors',varargout={donors(varargin{:})};
    case 'ray',varargout={ray_limit(varargin{:})};
    otherwise,error('ipm_cr:StaircaseAction','Unknown staircase operation.');
end
end

function c=build(c,cover,h,w)
a=ipm_cr.poisson.fem.axis(c.x,5,8,'native',false);
b=ipm_cr.poisson.fem.axis(c.y,5,8,'native',false);
x=c.x;y=c.y;
touch=hypot(y(1:end-1),min(abs(x(1:end-1)),abs(x(2:end))))<cover;
required=false(size(c.X));
for row=0:1
    for col=0:1
        i=(1:size(touch,1))+row;j=(1:size(touch,2))+col;
        required(i,j)=required(i,j)|touch;
    end
end
% The active mask is set by transport overlap. FEM support encloses it
% by whole native elements, including the Cartesian transport halo.
active=extend_strips(required,w,c.originIndex);
assert(all(c.active(active)),'ipm_cr:StairCoverage', ...
    'The staircase and its directional stencils must fit the Cartesian interior.');
padding=conv2(double(active),ones(1,2*h+1),'same')>0| ...
    conv2(double(active),ones(2*h+1,1),'same')>0;
elements=false(b.elementCount,a.elementCount);support=false(size(c.X));
for j=1:a.elementCount
    cols=a.connectivity{j};
    for i=1:b.elementCount
        rows=b.connectivity{i};
        if any(padding(rows,cols),'all')
            elements(i,j)=true;support(rows,cols)=true;
        end
    end
end
cells=active(1:end-1,1:end-1)&active(2:end,1:end-1)& ...
    active(1:end-1,2:end)&active(2:end,2:end);
assert(all(cells(touch))&&isequal(active,fliplr(active)));
top=zeros(1,c.nx-c.originIndex);
for j=1:numel(top)
    row=find(cells(:,c.originIndex+j-1),1,'last');
    if ~isempty(row),top(j)=y(row+1);end
end
c.active=active;c.ghost=~active;c.activeIndices=find(active);c.ghostIndices=find(~active);
c.activeRows=find(any(active,2))';c.activeColumns=find(any(active,1));
c.activeXBounds=x(c.activeColumns([1,end]));c.activeYBounds=y(c.activeRows([1,end]));
c.stair=struct('edges',x(c.originIndex:end),'top',top,'elements',elements, ...
    'support',support,'cover',cover,'degree',5,'quadOrder',8);
c.stair.dx=directional(c,'x',w);c.stair.dy=directional(c,'y',w);
end

function hi=ray_limit(s,theta)
% In each vertical slab the outer radius is bounded by the right edge
% and by the top step. The union is star-shaped about the wall origin.
theta=theta(:);
co=max(abs(cos(theta)),realmin);si=max(sin(theta),realmin);
lo=s.edges(1:end-1)./co;
upper=min(s.edges(2:end)./co,s.top./si);
upper(upper<lo|s.top==0)=0;
hi=max(upper,[],2);
end

function mask=extend_strips(mask,w,center)
d=(w-1)/2;
while true
    rows=find(sum(mask,2)>0&sum(mask,2)<w);
    cols=find(sum(mask,1)>0&sum(mask,1)<w);
    if isempty(rows)&&isempty(cols),return;end
    mask(rows,center+(-d:d))=true;mask(1:w,cols)=true;
end
end

function yes=contains_points(s,x,y)
ax=abs(x);edge=s.edges;tol=64*eps(max(1,ax));
ax=max(0,min(edge(end),ax-tol));
k=discretize(ax,edge);
k(isnan(k))=1;
height=reshape(s.top(k),size(y));
yes=isfinite(x)&isfinite(y)&abs(x)<=edge(end)+tol& ...
    y>=0&height>0&y<=height+64*eps(max(1,height));
end

function M=donors(c,x,y,w)
xAxis=c.x;yAxis=c.y;
pack=zeros(size(c.X));pack(c.activeIndices)=1:numel(c.activeIndices);
rows=zeros(numel(x),w*w);cols=rows;values=rows;
tol=64*eps(max(1,max(abs([x(:);y(:)]))));
for q=1:numel(x)
    i=find(xAxis(1:end-w+1)<=x(q)+tol&xAxis(w:end)>=x(q)-tol);
    j=find(yAxis(1:end-w+1)<=y(q)+tol&yAxis(w:end)>=y(q)-tol);
    [~,u]=sort(abs(xAxis(i+(w-1)/2)-x(q)));i=i(u);
    [~,v]=sort(abs(yAxis(j+(w-1)/2)-y(q)));j=j(v);
    found=false;
    for a=i(:)'
        ii=a+(0:w-1);
        for b=j(:)'
            jj=b+(0:w-1);ids=pack(jj,ii);
            if all(ids>0,'all'),found=true;break;end
        end
        if found,break;end
    end
    assert(found,'ipm_cr:StairDonor','The polar halo must have a bracketed active 7-by-7 tensor stencil.');
    wx=ipm_cr.mesh.localMap(xAxis(ii),x(q),w);
    wy=ipm_cr.mesh.localMap(yAxis(jj),y(q),w);
    z=full(wy'*wx);
    rows(q,:)=q;cols(q,:)=ids(:);values(q,:)=z(:);
end
M=sparse(rows(:),cols(:),values(:),numel(x),numel(c.activeIndices));
end

function P=directional(c,kind,w)
rows=[];cols=[];values=[];mask=c.active;index=reshape(1:numel(mask),size(mask));
if strcmp(kind,'x'),axis=c.x;mask=mask';index=index';else,axis=c.y;end
for j=1:size(mask,2)
    i=find(mask(:,j));if isempty(i),continue;end
    assert(numel(i)>=w&&all(diff(i)==1));
    [u,v,d]=find(ipm_cr.mesh.fd(axis(i),1,w));ids=index(i,j);
    rows=[rows;ids(u)];cols=[cols;ids(v)];values=[values;d]; %#ok<AGROW>
end
P=struct('i',rows,'j',cols,'D',sparse(rows,(1:numel(rows))',values,numel(mask),numel(rows)));
end

function d=differentiate(u,c)
P=c.stair.dx;v=u(:);x=reshape(P.D*(v(P.j)-v(P.i)),size(u));
P=c.stair.dy;y=reshape(P.D*(v(P.j)-v(P.i)),size(u));
d=struct('x',x,'y',y);
end
