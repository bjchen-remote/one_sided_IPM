function varargout=staircase(action,varargin)
%STAIRCASE Active Cartesian mask, native FEM support and active-only donors.
switch action
    case 'build',varargout={build(varargin{:})};
    case 'contains',varargout={contains_points(varargin{:})};
    case 'partition',varargout={partition(varargin{:})};
    case 'donors',varargout={donors(varargin{:})};
    case 'ray',varargout={ray_limit(varargin{:})};
    otherwise,error('ipm_cr:StaircaseAction','Unknown staircase operation.');
end
end

function c=build(c,cover,h,w)
a=partition(c.x);
b=partition(c.y);
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
left=floor((w-1)/2);right=ceil((w-1)/2);
for q=1:numel(x)
    i=find(xAxis(1:end-w+1)<=x(q)+tol&xAxis(w:end)>=x(q)-tol);
    j=find(yAxis(1:end-w+1)<=y(q)+tol&yAxis(w:end)>=y(q)-tol);
    xCenter=(xAxis(i+left)+xAxis(i+right))/2;
    yCenter=(yAxis(j+left)+yAxis(j+right))/2;
    [~,u]=sort(abs(xCenter-x(q)));i=i(u);
    [~,v]=sort(abs(yCenter-y(q)));j=j(v);
    found=false;
    for a=i(:)'
        ii=a+(0:w-1);
        for b=j(:)'
            jj=b+(0:w-1);ids=pack(jj,ii);
            if all(ids>0,'all'),found=true;break;end
        end
        if found,break;end
    end
    assert(found,'ipm_cr:StairDonor', ...
        'The query must have a bracketed active %d-by-%d tensor stencil.',w,w);
    wx=ipm_cr.mesh.interpolate.local(xAxis(ii),x(q),w);
    wy=ipm_cr.mesh.interpolate.local(yAxis(jj),y(q),w);
    z=full(wy'*wx);
    rows(q,:)=q;cols(q,:)=ids(:);values(q,:)=z(:);
end
M=sparse(rows(:),cols(:),values(:),numel(x),numel(c.activeIndices));
end

function out=partition(x)
% The original collar partition defines coverage on the native axes.
% Its cell endpoints select whole B-spline knot cells for the weak form.
x=double(x(:));p=5;n=numel(x)-1;ne=floor(n/p);extra=n-ne*p;
assert(ne>=1&&extra<=ne,'ipm_cr:FEMNativeGrid');
degree=p*ones(ne,1);
if extra>0
    degree(floor(((1:extra)-.5)*ne/extra)+1)=p+1;
end
k=(1:7)';off=k./sqrt(4*k.^2-1);[~,L]=eig(diag(off,1)+diag(off,-1));
degree=stable_partition(x,p,sort(diag(L)),degree);
edges=[1;1+cumsum(degree)];
conn=arrayfun(@(j) edges(j):edges(j+1),(1:ne)','UniformOutput',false);
out=struct('breaks',x(edges),'connectivity',{conn},'elementCount',ne);
end

function degrees=stable_partition(x,p,gq,initial)
% Moving only the p+1 elements changes neither the nodes nor their number.
% Minimize the worst Gauss-point Lebesgue factor by dynamic programming;
% ties retain the most evenly distributed placement. No source/solution
% value or measured PDE error participates in this geometric choice.
ne=numel(initial);extra=sum(initial-p);oldExtra=find(initial>p);
cost=inf(ne+1,extra+1);distance=cost;choice=zeros(ne+1,extra+1);
cost(1,1)=0;distance(1,1)=0;
for e=0:ne-1
    for used=0:min(e,extra)
        if ~isfinite(cost(e+1,used+1)),continue;end
        first=1+e*p+used;
        for take=0:1
                after=used+take;
            if after>extra||extra-after>ne-e-1,continue;end
            value=local_condition(x,first,p+take,gq);
            worst=max(cost(e+1,used+1),value);
            displacement=distance(e+1,used+1);
            if take,displacement=displacement+abs((e+1)-oldExtra(after));end
            previous=cost(e+2,after+1);tolerance=1e-12*max(1,worst);
            if ~isfinite(previous)||worst<previous-tolerance|| ...
                    (abs(worst-previous)<=tolerance&&displacement<distance(e+2,after+1))
                cost(e+2,after+1)=worst;distance(e+2,after+1)=displacement;
                choice(e+2,after+1)=take;
            end
        end
    end
end
degrees=p*ones(ne,1);used=extra;
for e=ne:-1:1
    take=choice(e+1,used+1);degrees(e)=p+take;used=used-take;
end
assert(used==0&&sum(degrees)==numel(x)-1,'ipm_cr:FEMNativeGrid','Native element partition failed.');
end
function [score,gapRatio]=local_condition(x,first,degree,gq)
values=x(first:first+degree);h=values(end)-values(1);
reference=2*(values-values(1))/h-1;reference([1,end])=[-1,1];
B=geometry_polynomials(reference,gq);
score=norm(B,inf);gap=diff(reference);gapRatio=max(gap)/min(gap);
end


function [B,D]=geometry_polynomials(nodes,x)
validateattributes(nodes,{'numeric'},{'vector','real','finite','nonempty'});
validateattributes(x,{'numeric'},{'real','finite'});
nodes=double(nodes(:));x=double(x(:));n=numel(nodes);
assert(n>=2&&all(diff(nodes)>0),'ipm_cr:FEMLocalNodes', ...
    'Use distinct increasing local interpolation nodes.');
B=ones(numel(x),n);D=zeros(numel(x),n);
% Product-rule recurrence differentiates the Lagrange polynomial itself.
% It remains regular at interpolation nodes and never forms a finite-
% difference derivative or a pointwise second-derivative operator.
for i=1:n
    value=ones(size(x));derivative=zeros(size(x));
    for j=1:n
        if j==i,continue;end
        denominator=nodes(i)-nodes(j);factor=(x-nodes(j))/denominator;
        derivative=derivative.*factor+value/denominator;
        value=value.*factor;
    end
    B(:,i)=value;D(:,i)=derivative;
end
end
