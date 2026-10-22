function cart=cartesian(data,box,resolution,options)
%CARTESIAN Resample one saved polar slice onto a rectangular Cartesian grid.
% c=ipm_c.cartesian(d,[xmin xmax ymin ymax],[nx ny],options)
% d is a selected slice returned by ipm_c.data. resolution may be one integer.
% options.fields defaults to all available core fields; coordinates is
% 'native' (default) or 'physical'; method is 'cubic' (default) or 'linear'.
% Cubic uses local P3 Hermite/Bernstein finite elements on the actual polar
% triangles. Vertex values and shared vertex gradients determine every edge;
% a local cubic least-squares centroid reconstruction determines the bubble.
% It is C0 across triangle edges, not a C1 element and not shape preserving.
% Linear is the P1 finite element alternative without interpolation overshoot.
% The repeated origin row is one vertex. No extrapolation is performed across
% the wall or beyond the straight outer polygon, including its circular sliver.
% Physical coordinates map x=(X-X_shift)/C_l, y=Y/C_l. Physical rho,
% gradients, streamfunction and Darcy velocity receive their proper units.
% F/rho_tau, Q/rho_xt and vx/vy retain CANONICAL units, even when plotted at
% physical positions; metadata.fieldUnits records this explicitly.
% Optional options.gradients.<field>.X/Y provide native nodal derivatives;
% otherwise saved rho_X/Y or seven-point polar derivatives are used.
if nargin<4,options=struct();end
if nargin<3||isempty(resolution),resolution=[401,401];end
validateattributes(data,{'struct'},{'scalar'});
validateattributes(box,{'numeric'},{'real','finite','vector','numel',4});
box=double(box(:)');assert(box(1)<box(2)&&box(3)<box(4), ...
    'ipm_c:CartesianBox','Use [xmin xmax ymin ymax] with increasing bounds.');
if isscalar(resolution),resolution=[resolution,resolution];end
validateattributes(resolution,{'numeric'},{'real','finite','integer','vector','numel',2,'>=',2});
resolution=double(resolution(:)');
coordinates=char(option(options,'coordinates','native'));
method=char(option(options,'method','cubic'));
assert(ismember(coordinates,{'native','physical'}),'ipm_c:CartesianCoordinates', ...
    'coordinates must be native or physical.');
assert(ismember(method,{'cubic','linear'}),'ipm_c:CartesianMethod', ...
    'method must be cubic or linear.');
d=selectedData(data);r=double(d.r(:));theta=double(d.theta(:)');
nr=numel(r);nt=numel(theta);
assert(nr>=4&&nt>=4&&r(1)==0&&all(isfinite(r))&&all(diff(r)>0)&& ...
    all(isfinite(theta))&&all(diff(theta)>0)&&theta(1)==0&& ...
    abs(theta(end)-pi)<=32*eps(pi),'ipm_c:CartesianGeometry', ...
    'The selected frame must contain its actual increasing half-disk axes.');
ct=cos(theta);st=sin(theta);ct([1,end])=[1,-1];st([1,end])=0;
X=r*ct;Y=r*st;
if mod(nt,2)==1&&abs(theta((nt+1)/2)-pi/2)<=32*eps(pi)
    X(:,(nt+1)/2)=0;
end
scale=option(d,'scale',struct('logC_l',0,'logC_omega',0,'X_shift',0));
clength=exp(scale.logC_l);comega=exp(scale.logC_omega);
shift=option(scale,'X_shift',0);
if strcmp(coordinates,'physical')
    assert(isfinite(clength)&&clength>0&&isfinite(comega)&&comega>0&&isfinite(shift), ...
        'ipm_c:CartesianPhysicalScale','The selected physical scales must be representable.');
end
x=linspace(box(1),box(2),resolution(1));y=linspace(box(3),box(4),resolution(2))';
[qx,qy]=meshgrid(x,y);
if strcmp(coordinates,'physical'),nativeX=clength*qx+shift;nativeY=clength*qy;
else,nativeX=qx;nativeY=qy;end
[inside,triangle,vertices,bary,cellR,cellT]=locate(nativeX,nativeY,r,theta,X,Y);
nodeX=[0;reshape(X(2:end,:),[],1)];nodeY=[0;reshape(Y(2:end,:),[],1)];
names=option(options,'fields',{});
if isempty(names)
    candidates={'rho','rho_X','rho_Y','psi','u1','u2','vx','vy','rho_tau','Q'};
    names=candidates(cellfun(@(name) isfield(d,name)&&all(isfinite(d.(name)),'all'),candidates));
elseif ischar(names)||isstring(names),names=cellstr(names);end
assert(iscell(names)&&all(cellfun(@ischar,names)),'ipm_c:CartesianFields', ...
    'fields must contain field names.');
assert(~isempty(names),'ipm_c:CartesianFields','No requested fields are available.');
names=unique(names,'stable');nfields=numel(names);nnode=numel(nodeX);
values=zeros(nnode,nfields);gradX=values;gradY=values;offsets=zeros(1,nfields);
units=struct();gradientSource=struct();
for k=1:nfields
    name=names{k};canonical=canonicalName(name);
    if strcmp(canonical,'rho_t')
        assert(strcmp(coordinates,'physical'),'ipm_c:CartesianRhoT', ...
            'rho_t denotes the physical Eulerian time derivative; select physical coordinates.');
        f=physicalRhoT(d,X,Y,clength,comega);
        unit='physical d(rho)/dt at fixed physical x,y';factor=1;
    else
        assert(isfield(d,canonical),'ipm_c:CartesianMissingField', ...
            'The selected slice does not contain %s.',canonical);
        f=double(d.(canonical));[factor,unit]=fieldUnits(canonical,coordinates,clength,comega);
        if strcmp(canonical,'rho')&&isfield(d,'rho_centered')
            f=double(d.rho_centered);offsets(k)=d.rho_origin*factor;
        elseif strcmp(canonical,'rho_tau')&&isfield(d,'rho_tau_centered')
            f=double(d.rho_tau_centered);offsets(k)=d.rho_origin_tau*factor;
        end
    end
    assert(isequal(size(f),[nr,nt])&&isreal(f)&&all(isfinite(f),'all'), ...
        'ipm_c:CartesianField','%s must be a finite real nr-by-ntheta field.',name);
    assert(all(f(1,:)==f(1,1)),'ipm_c:CartesianOrigin', ...
        '%s must have one identical physical-origin value.',name);
    values(:,k)=uniqueNodes(f)*factor;units.(name)=unit;
    if strcmp(method,'cubic')
        [fx,fy,source]=fieldGradient(f,d,canonical,r,theta,X,Y,options);
        gradX(:,k)=uniqueNodes(fx)*factor;gradY(:,k)=uniqueNodes(fy)*factor;
        gradientSource.(name)=source;
    end
end
cart=struct('kind','ipm_c_cartesian_slice_v1','x',x,'y',y,'X',qx,'Y',qy, ...
    'mask',inside,'box',box,'resolution',resolution,'coordinates',coordinates, ...
    'scale',scale,'selected',option(data,'selected',struct()));
if isfield(data,'meta')&&isfield(data.meta,'datasetFile'),cart.source=data.meta.datasetFile;end
index=find(inside);vertexQuery=vertices;lambda=bary;
if isempty(index)
    reconstruction=struct('cubicTriangles',0,'quadraticTriangles',0,'constantFallbackTriangles',0);
else
    [used,first,which]=unique(triangle,'stable'); %#ok<ASGLU>
    vertexUsed=vertices(first,:);nused=numel(used);
    if strcmp(method,'cubic')
        [weights,reconstruction]=centroidWeights(vertexUsed,cellR(first),cellT(first), ...
            r,theta,nodeX,nodeY);
        [weightRows,weightNodes,weightValues]=find(weights);
        center=zeros(nused,nfields);
        for field=1:nfields
            reference=values(vertexUsed(weightRows,1),field);
            center(:,field)=accumarray(weightRows, ...
                weightValues.*(values(weightNodes,field)-reference),[nused,1],@sum,0);
        end
    else,reconstruction=struct('cubicTriangles',0,'quadraticTriangles',0,'constantFallbackTriangles',0);end
end
for k=1:nfields
    out=nan(size(qx));
    if ~isempty(index)
        f=values(:,k);reference=f(vertexQuery(:,1));
        a=zeros(size(reference));b=f(vertexQuery(:,2))-reference;c=f(vertexQuery(:,3))-reference;
        l1=lambda(:,1);l2=lambda(:,2);l3=lambda(:,3);
        if strcmp(method,'linear'),z=a.*l1+b.*l2+c.*l3;
        else
            v=vertexUsed;vertexReference=f(v(:,1));f1=zeros(size(vertexReference));
            f2=f(v(:,2))-vertexReference;f3=f(v(:,3))-vertexReference;
            gx=gradX(:,k);gy=gradY(:,k);
            b210=f1+direction(gx,gy,v(:,1),v(:,2),nodeX,nodeY)/3;
            b201=f1+direction(gx,gy,v(:,1),v(:,3),nodeX,nodeY)/3;
            b120=f2+direction(gx,gy,v(:,2),v(:,1),nodeX,nodeY)/3;
            b021=f2+direction(gx,gy,v(:,2),v(:,3),nodeX,nodeY)/3;
            b102=f3+direction(gx,gy,v(:,3),v(:,1),nodeX,nodeY)/3;
            b012=f3+direction(gx,gy,v(:,3),v(:,2),nodeX,nodeY)/3;
            bubble=(27*center(:,k)-(f1+f2+f3)-3*(b210+b201+b120+b021+b102+b012))/6;
            z=a.*l1.^3+b.*l2.^3+c.*l3.^3+3*( ...
                b210(which).*l1.^2.*l2+b201(which).*l1.^2.*l3+ ...
                b120(which).*l1.*l2.^2+b021(which).*l2.^2.*l3+ ...
                b102(which).*l1.*l3.^2+b012(which).*l2.*l3.^2)+ ...
                6*bubble(which).*l1.*l2.*l3;
        end
        z=z+reference+offsets(k);
        assert(all(isfinite(z)),'ipm_c:CartesianNonfinite', ...
            'Interpolation of %s is not representable.',names{k});
        out(index)=z;
    end
    cart.(names{k})=out;
end
cart.metadata=struct('method',method,'element',elementName(method), ...
    'continuity','C0; gradients need not agree across triangle edges', ...
    'originIsSingleVertex',true,'extrapolation','none', ...
    'domainMask','intersection of upper half-disk and actual straight-sided polar triangle mesh', ...
    'fieldUnits',units,'gradientSource',gradientSource,'centroidReconstruction',reconstruction, ...
    'nativeR',r,'nativeTheta',theta,'nativeRadius',r(end), ...
    'canonicalDerivativeWarning','rho_tau/F, Q/rho_xt and vx/vy remain canonical quantities at either coordinate choice');
end

function d=selectedData(data)
d=data;
for name={'current','frame'}
    if isfield(data,name{1})&&isstruct(data.(name{1}))
        nested=data.(name{1});fields=fieldnames(nested);
        for k=1:numel(fields),d.(fields{k})=nested.(fields{k});end
    end
end
if ~isfield(d,'scale')&&isfield(d,'selected')&&isfield(d.selected,'scale'),d.scale=d.selected.scale;end
if ~isfield(d,'rho_X')&&isfield(d,'rho_x'),d.rho_X=d.rho_x;end
if ~isfield(d,'rho_Y')&&isfield(d,'rho_y'),d.rho_Y=d.rho_y;end
if ~isfield(d,'rho_tau')&&isfield(d,'F'),d.rho_tau=d.F;end
if ~isfield(d,'Q')&&isfield(d,'rho_xt'),d.Q=d.rho_xt;end
end

function name=canonicalName(name)
switch name
    case 'rho_x',name='rho_X';case 'rho_y',name='rho_Y';
    case 'F',name='rho_tau';case 'rho_xt',name='Q';
end
end

function [factor,unit]=fieldUnits(name,coordinates,cl,cw)
factor=1;unit='native canonical field';
switch name
    case 'rho',unit='native R';if strcmp(coordinates,'physical'),factor=1/cw;unit='physical rho';end
    case {'rho_X','rho_Y'},unit='native derivative dR/dX or dR/dY';
        if strcmp(coordinates,'physical'),factor=cl/cw;unit='physical derivative d(rho)/dx or d(rho)/dy';end
    case {'u1','u2','vmag'},unit='native Darcy velocity U';
        if strcmp(coordinates,'physical'),factor=1/cw;unit='physical Darcy velocity u';end
    case 'psi',unit='native streamfunction Psi';
        if strcmp(coordinates,'physical'),factor=1/(cl*cw);unit='physical streamfunction psi';end
    case {'vx','vy'},unit='canonical composite transport velocity d(X,Y)/d(tau)';
    case 'rho_tau',unit='canonical dR/dtau at fixed native X,Y';
    case 'Q',unit='canonical d/dX(dR/dtau); not physical rho_xt';
end
end

function f=physicalRhoT(d,X,Y,cl,cw)
for name={'rho','rho_X','rho_Y','rho_tau'}
    assert(isfield(d,name{1}),'ipm_c:CartesianPhysicalRhoT','rho_t requires %s.',name{1});
end
rates=option(d,'selected',struct());
if isfield(rates,'rates'),rates=rates.rates;end
for name={'c_l','c_omega','c_r'}
    if isfield(d,name{1})&&isscalar(d.(name{1})),rates.(name{1})=d.(name{1});end
    assert(isfield(rates,name{1})&&isscalar(rates.(name{1}))&&isfinite(rates.(name{1})), ...
        'ipm_c:CartesianPhysicalRhoT','rho_t requires the selected scalar %s.',name{1});
end
if all(isfield(d,{'rho_centered','rho_tau_centered','transportOrigin'}))
    % a_tau-c_omega*a is the small transport origin T(W)(0), retained
    % explicitly rather than recovered from two large amplitude terms.
    fixedCoordinateRate=d.rho_tau_centered+d.transportOrigin-rates.c_omega*d.rho_centered;
else
    fixedCoordinateRate=d.rho_tau-rates.c_omega*d.rho;
end
f=(cl/cw^2)*(fixedCoordinateRate+(rates.c_l*X+rates.c_r).*d.rho_X+ ...
    rates.c_l*Y.*d.rho_Y);
end

function [mask,triangle,v,lambda,radial,angular]=locate(qx,qy,r,theta,X,Y)
% A polar sector's straight radial chords satisfy q.normal=r*cos(dtheta/2).
% This supplies exact structured triangle location without a global Delaunay
% triangulation, whose scales would span roughly twelve orders of magnitude.
shape=size(qx);nr=numel(r);nt=numel(theta);qx=qx(:);qy=qy(:);
a=atan2(qy,qx);mask=qy>=0&isfinite(qx)&isfinite(qy)&hypot(qx,qy)<=r(end);
angular=discretize(a,theta);angular(a==theta(end))=nt-1;
mask=mask&~isnan(angular);ids=find(mask);j=angular(ids);
mid=(theta(j)+theta(j+1))'/2;half=(theta(j+1)-theta(j))'/2;
qr=(qx(ids).*cos(mid)+qy(ids).*sin(mid))./cos(half);
ir=discretize(qr,r);ir(qr==r(end))=nr-1;
valid=~isnan(ir);mask(ids(~valid))=false;ids=ids(valid);j=j(valid);ir=ir(valid);
map=[ones(1,nt);reshape(2:1+(nr-1)*nt,nr-1,nt)];
v=zeros(numel(ids),3);triangle=zeros(numel(ids),1);
fan=ir==1;
v(fan,:)=[ones(sum(fan),1),map(sub2ind([nr,nt],2*ones(sum(fan),1),j(fan))), ...
    map(sub2ind([nr,nt],2*ones(sum(fan),1),j(fan)+1))];
triangle(fan)=j(fan);
other=find(~fan);i=ir(other);jj=j(other);
aidx=map(sub2ind([nr,nt],i,jj));bidx=map(sub2ind([nr,nt],i+1,jj));
cidx=map(sub2ind([nr,nt],i+1,jj+1));didx=map(sub2ind([nr,nt],i,jj+1));
nx=[0;reshape(X(2:end,:),[],1)];ny=[0;reshape(Y(2:end,:),[],1)];
cross=(nx(cidx)-nx(aidx)).*(qy(ids(other))-ny(aidx))- ...
    (ny(cidx)-ny(aidx)).*(qx(ids(other))-nx(aidx));
second=cross>=0;v(other,:)=[aidx,bidx,cidx];v(other(second),:)=[aidx(second),cidx(second),didx(second)];
triangle(other)=(nt-1)+2*((jj-1)*(nr-2)+(i-2))+1+second;
ax=nx(v(:,1));ay=ny(v(:,1));bx=nx(v(:,2));by=ny(v(:,2));cx=nx(v(:,3));cy=ny(v(:,3));
den=(bx-ax).*(cy-ay)-(by-ay).*(cx-ax);
l2=((qx(ids)-ax).*(cy-ay)-(qy(ids)-ay).*(cx-ax))./den;
l3=((bx-ax).*(qy(ids)-ay)-(by-ay).*(qx(ids)-ax))./den;
lambda=[1-l2-l3,l2,l3];
valid=all(lambda>=-256*eps,2)&all(lambda<=1+256*eps,2)&all(isfinite(lambda),2);
mask(ids(~valid))=false;v=v(valid,:);lambda=lambda(valid,:);triangle=triangle(valid);
radial=ir(valid);angular=j(valid);
lambda=max(0,min(1,lambda));lambda=lambda./sum(lambda,2);
mask=reshape(mask,shape);
end

function [fx,fy,source]=fieldGradient(f,d,name,r,theta,X,Y,options)
if isfield(options,'gradients')&&isfield(options.gradients,name)
    provided=options.gradients.(name);fx=provided.X;fy=provided.Y;source='provided native nodal gradients';
elseif strcmp(name,'rho')&&isfield(d,'rho_X')&&isfield(d,'rho_Y')
    fx=d.rho_X;fy=d.rho_Y;source='saved native rho_X and rho_Y';
else
    geometry=gradientGeometry(r,theta,X,Y);
    fr=geometry.Dr*(f-f(1,1));ft=(f-f(:,1))*geometry.Dtheta';
    inverse=1./r(2:end);ct=cos(theta);st=sin(theta);ct([1,end])=[1,-1];st([1,end])=0;
    fx=zeros(size(f));fy=fx;
    fx(2:end,:)=fr(2:end,:).*ct-(ft(2:end,:).*inverse).*st;
    fy(2:end,:)=fr(2:end,:).*st+(ft(2:end,:).*inverse).*ct;
    residual=f(geometry.rows,:)-f(1,1);coeff=geometry.originMap*residual(:);
    fx(1,:)=coeff(1)/geometry.scale;fy(1,:)=coeff(2)/geometry.scale;
    source='seven-point polar finite differences and one Cartesian cubic origin fit';
end

assert(isequal(size(fx),size(f))&&isequal(size(fy),size(f))&& ...
    all(isfinite(fx),'all')&&all(isfinite(fy),'all'),'ipm_c:CartesianGradient', ...
    'Native nodal gradients of %s must be finite and match its field shape.',name);
assert(all(fx(1,:)==fx(1,1))&&all(fy(1,:)==fy(1,1)), ...
    'ipm_c:CartesianOriginGradient','The physical-origin Cartesian gradient must be shared.');
end

function geometry=gradientGeometry(r,theta,X,Y)
% Cache only disposable plotting operators; actual axes key every reuse.
persistent previous
if ~isempty(previous)&&isequal(previous.r,r)&&isequal(previous.theta,theta)
    geometry=previous;return
end
width=min([7,numel(r),numel(theta)]);rows=2:min(5,numel(r));scale=r(rows(end));
xx=X(rows,:)/scale;yy=Y(rows,:)/scale;xx=xx(:);yy=yy(:);
design=[xx,yy,xx.^2,xx.*yy,yy.^2,xx.^3,xx.^2.*yy,xx.*yy.^2,yy.^3];
[Q,R]=qr(design,0);
geometry=struct('r',r,'theta',theta,'Dr',ipm.mesh.fdMatrix(r,1,width), ...
    'Dtheta',ipm.mesh.fdMatrix(theta,1,width),'rows',rows,'scale',scale,'originMap',R\Q');
previous=geometry;
end

function [weights,info]=centroidWeights(v,ir,jt,r,theta,nx,ny)
nr=numel(r);nt=numel(theta);n=size(v,1);maxEntries=64*n;
ii=zeros(maxEntries,1);jj=ii;vv=ii;position=0;counts=zeros(1,3);
map=[ones(1,nt);reshape(2:1+(nr-1)*nt,nr-1,nt)];
for k=1:n
    radial=max(2,ir(k)-2):min(nr,max(5,ir(k)+3));
    angular=max(1,jt(k)-2):min(nt,jt(k)+3);
    ids=map(radial,angular);ids=ids(:);
    if ir(k)<=3,ids=[1;ids];end
    centerX=mean(nx(v(k,:)));centerY=mean(ny(v(k,:)));
    angle=(theta(jt(k))+theta(jt(k)+1))/2;
    dx=nx(ids)-centerX;dy=ny(ids)-centerY;
    xx=dx*cos(angle)+dy*sin(angle);yy=-dx*sin(angle)+dy*cos(angle);
    xx=xx/max(abs(xx));yy=yy/max(abs(yy));
    design=[ones(size(xx)),xx,yy,xx.^2,xx.*yy,yy.^2, ...
        xx.^3,xx.^2.*yy,xx.*yy.^2,yy.^3];
    [Q,R]=qr(design,0);kind=1;
    if size(R,1)<10||rcond(R)<1e-10
        design=design(:,1:6);[Q,R]=qr(design,0);kind=2;
    end
    if size(R,1)<6||rcond(R)<1e-12
        ids=v(k,:)';w=ones(3,1)/3;kind=3;
    else
        e=zeros(size(R,1),1);e(1)=1;w=Q*(R'\e);
        % Re-enforce constants despite the last few least-squares ulps.
        w=w/sum(w);
    end
    counts(kind)=counts(kind)+1;take=position+(1:numel(ids));
    ii(take)=k;jj(take)=ids;vv(take)=w;position=position+numel(ids);
end
weights=sparse(ii(1:position),jj(1:position),vv(1:position),n,numel(nx));
info=struct('cubicTriangles',counts(1),'quadraticTriangles',counts(2), ...
    'constantFallbackTriangles',counts(3));
end

function result=direction(gx,gy,a,b,x,y)
result=gx(a).*(x(b)-x(a))+gy(a).*(y(b)-y(a));
end

function out=uniqueNodes(f)
out=[f(1,1);reshape(f(2:end,:),[],1)];
end

function name=elementName(method)
if strcmp(method,'cubic'),name='P3 Hermite/Bernstein triangle with local cubic centroid reconstruction';
else,name='P1 barycentric triangle';end
end

function value=option(options,name,fallback)
if isfield(options,name),value=options.(name);else,value=fallback;end
end
