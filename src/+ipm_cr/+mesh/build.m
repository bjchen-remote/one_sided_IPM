function g=build(options)
%GEOMETRY Tensor Cartesian and origin-free mapped polar transport grids.
% Both overlapping transport interiors evolve. Three Cartesian side/top
% collars and three inward annular rows are interpolation ghosts only.
% preset='medium'/'coarse' loads portable physical axes, with coarse exactly
% every other medium node. axes supplies x/y/r/theta and positive metric
% steps dxMetric/dyMetric/drMetric for remeshing; no state is transferred.
if nargin<1,options=struct();end
assert(isstruct(options)&&isscalar(options),'ipm_cr:ExplicitGeometry', ...
    'Geometry options must be a scalar structure.');
allowed={'grid','coreHalfWidth','coreHeight','coreSpacing','polarStart', ...
    'ghostWidth','interpolationWidth','radialRedistribution','coreMap','preset','axes'};
assert(all(ismember(fieldnames(options),allowed)),'ipm_cr:ExplicitGeometry', ...
    'An unknown explicit transport geometry option was supplied.');
timer=tic;
explicit=isfield(options,'axes')||isfield(options,'preset');
halo=option(options,'ghostWidth',3);width=option(options,'interpolationWidth',7);
validateattributes(halo,{'numeric'},{'scalar','integer','>=',3});
validateattributes(width,{'numeric'},{'scalar','integer','>=',3});
if explicit
    [axes,axisSource]=explicit_axes(options);
    x=axes.x;y=axes.y;intervalX=numel(x)-1;intervalY=numel(y)-1;
    a=x(end);height=y(end);mapped=true;
    grid=struct('nr',numel(axes.r),'ntheta',numel(axes.theta), ...
        'rmax',axes.r(end),'firstSpacing',axes.drMetric(1), ...
        'radialMap','explicit_axes','stencilWidth',7);
    start=axes.r(1);
else
grid=struct('nr',385,'ntheta',128,'rmax',1e10,'firstSpacing',.03, ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01, ...
    'stencilWidth',7);
if isfield(options,'grid')
    assert(isstruct(options.grid)&&isscalar(options.grid),'ipm_cr:ExplicitGeometry', ...
        'grid must contain native ipm_cr.mesh.polar overrides.');
    names=fieldnames(options.grid);
    for k=1:numel(names),grid.(names{k})=options.grid.(names{k});end
end
a=option(options,'coreHalfWidth',1.2);height=option(options,'coreHeight',1.2);
h=option(options,'coreSpacing',.03);start=option(options,'polarStart',.72);
validateattributes(a,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(height,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(h,{'numeric'},{'vector','real','finite','positive'});
assert(ismember(numel(h),[1,2]),'ipm_cr:ExplicitGeometry','Use one or two Cartesian spacings.');
if isscalar(h),h=[h,h];end
validateattributes(start,{'numeric'},{'scalar','real','finite','positive'});
map=option(options,'coreMap',struct('enabled',false));
assert(isstruct(map)&&isscalar(map)&&all(ismember(fieldnames(map), ...
    {'enabled','kind','intervals','minSpacing'})), ...
    'ipm_cr:ExplicitGeometry','Unknown coreMap input.');
mapped=option(map,'enabled',false);validateattributes(mapped,{'logical'},{'scalar'});
if mapped
    axes=ipm_cr.mesh.tensorAxis(a,height,map);
    intervalX=axes.coreMap.intervals(1);intervalY=axes.coreMap.intervals(2);
    x=axes.x;y=axes.y;
else
    % Keep the original uniform constructor and its arithmetic unchanged.
    intervalX=round(2*a/h(1));intervalY=round(height/h(2));
    assert(abs(intervalX*h(1)-2*a)<=128*eps(max(1,2*a))&& ...
        abs(intervalY*h(2)-height)<=128*eps(max(1,height)), ...
        'ipm_cr:ExplicitGeometry','The Cartesian dimensions must be integer multiples of their spacing.');
    x=linspace(-a,a,intervalX+1);y=linspace(0,height,intervalY+1)';
    if mod(intervalX,2)==0,x(intervalX/2+1)=0;end
end
end
[X,Y]=meshgrid(x,y);nx=numel(x);ny=numel(y);
ix=(halo+1):(nx-halo);iy=1:(ny-halo);
assert(numel(ix)>=max(width,7)&&numel(iy)>=max(width,7),'ipm_cr:ExplicitGeometry', ...
    'The active Cartesian axes must fit the interpolation and D7 stencils.');
active=false(ny,nx);active(iy,ix)=true;
core=struct('x',x,'y',y,'X',X,'Y',Y,'nx',nx,'ny',ny, ...
    'dx',2*a/intervalX,'dy',height/intervalY,'hx',2*a/intervalX,'hy',height/intervalY, ...
    'Dx',active_derivative(x,ix),'Dy',active_derivative(y,iy), ...
    'active',active,'ghost',~active,'activeRows',iy,'activeColumns',ix, ...
    'activeIndices',find(active),'ghostIndices',find(~active), ...
    'activeXBounds',x(ix([1,end])),'activeYBounds',y(iy([1,end])), ...
    'ghostWidth',halo,'physicalWallRow',1);
[~,core.originIndex]=min(abs(x));
assert(x(core.originIndex)==0,'ipm_cr:ExplicitGeometry','An exact Cartesian origin is required.');
if mapped
    core.Jx=axes.Jx;core.Jy=axes.Jy;
    core.dxMetric=axes.dxMetric;core.dyMetric=axes.dyMetric;
    core.dxi=axes.dxi;core.deta=axes.deta;core.xi=axes.xi;core.eta=axes.eta;
    core.coreMap=axes.coreMap;core.mapParameters=axes.mapParameters;
    % Scalar compatibility fields denote minimum mapped metric spacing.
    % Every mapped transport/CFL operation below uses the complete vectors.
    core.dx=min(core.dxMetric);core.dy=min(core.dyMetric);
    core.hx=core.dx;core.hy=core.dy;
end

if explicit
    r=axes.r;theta=axes.theta;drMetric=axes.drMetric;
    co=cos(theta);si=sin(theta);co([1,end])=[1,-1];si([1,end])=0;
    if mod(numel(theta),2)==1
        middle=(numel(theta)+1)/2;co(middle)=0;si(middle)=1;
    end
    original=struct('cosTheta',co,'sinTheta',si);
    fullX=r*co;fullY=r*si;fullNr=numel(r);selected=(1:fullNr)';
    redistribute=false;
else
original=ipm_cr.mesh.polar(grid);
redistribution=option(options,'radialRedistribution',struct());
assert(isstruct(redistribution)&&isscalar(redistribution)&&all(ismember(fieldnames(redistribution), ...
    {'enabled','factor','preserveRadius','fullRefinementRadius'})), ...
    'ipm_cr:ExplicitGeometry','Unknown radial redistribution option.');
redistribute=option(redistribution,'enabled',false);
validateattributes(redistribute,{'numeric','logical'},{'scalar','real','finite'});
assert(any(redistribute==[0,1]),'ipm_cr:ExplicitGeometry','Redistribution enabled must be a logical flag.');
fullR=original.r;fullNr=original.nr;
fullDrMetric=original.radialMetric/(original.nr-1);
fullX=original.X;fullY=original.Y;
if redistribute
    [fullR,fullDrMetric,redistributionInfo]=ipm_cr.mesh.radial(original,redistribution);
    fullNr=numel(fullR);fullX=fullR*original.cosTheta;fullY=fullR*original.sinTheta;
    copied=redistributionInfo.preservedPrefixCount;
    fullX(1:copied,:)=original.X(1:copied,:);fullY(1:copied,:)=original.Y(1:copied,:);
end
first=find(fullR>=start-64*eps(max(1,start)),1,'first');
assert(~isempty(first)&&first>1&&fullNr-first+1>halo+width, ...
    'ipm_cr:ExplicitGeometry','The selected annulus must omit the origin and retain its halo/interior.');
selected=(first:fullNr)';r=fullR(selected);theta=original.theta;
drMetric=fullDrMetric(selected);
end
nr=numel(r);nt=numel(theta);
assert(nr>halo+width&&nr-halo>=7&&nt>=7,'ipm_cr:ExplicitGeometry', ...
    'The annulus must retain its halo and complete derivative/donor stencils.');
Jr=drMetric*(nr-1);
polarActive=false(nr,nt);polarActive((halo+1):end,:)=true;
Dt=ipm_cr.mesh.fd(theta,1,7);
polar=struct('r',r,'theta',theta,'X',fullX(selected,:), ...
    'Y',fullY(selected,:),'nr',nr,'ntheta',nt,'nt',nt, ...
    'Jr',Jr,'radialMetric',Jr,'drMetric',drMetric, ...
    'radialLogicalSpacing',1/(nr-1),'dtheta',pi/(nt-1),'angularSpacing',pi/(nt-1), ...
    'Dr',active_derivative(r,(halo+1):nr),'Dt',Dt,'Dtheta',Dt, ...
    'cosTheta',original.cosTheta,'sinTheta',original.sinTheta, ...
    'active',polarActive,'ghost',~polarActive,'activeRows',(halo+1):nr, ...
    'activeColumns',1:nt,'activeIndices',find(polarActive), ...
    'ghostIndices',find(~polarActive),'innerHaloRows',1:halo, ...
    'ghostWidth',halo,'originalIndices',selected,'originalNr',fullNr, ...
    'requestedInnerHaloRadius',start,'innerHaloRadius',r(1), ...
    'innerActiveRadius',r(halo+1),'domain','half_annulus');
polar.Y(:,[1,end])=0;

% Build the donor axes from ACTIVE values only. Neither interpolation map
% can read another ghost, so both collars can be exchanged simultaneously.
cx=X(core.ghostIndices);cy=Y(core.ghostIndices);
coreQueryR=hypot(cx,cy);coreQueryTheta=atan2(cy,cx);
assert(min(coreQueryR)>=r(halo+1)&&max(coreQueryR)<=r(end), ...
    'ipm_cr:ExplicitOverlap','Every Cartesian ghost must lie in the active annulus.');
pr=polar.X(polar.ghostIndices);py=polar.Y(polar.ghostIndices);
assert(all(pr>=x(ix(1))&pr<=x(ix(end)))&& ...
    all(py>=y(iy(1))&py<=y(iy(end))), ...
    'ipm_cr:ExplicitOverlap','Every annular inward ghost must lie in the active Cartesian rectangle.');
coreFromPolar=ipm_cr.mesh.tensorMap(r(polar.activeRows),theta, ...
    coreQueryR,coreQueryTheta,width);
polarFromCore=ipm_cr.mesh.tensorMap(y(iy),x(ix),py,pr,width);
if redistribute
    % Preserve the physical donors and weights, not sparse column strides:
    % the latter change when the complete annular radial dimension changes.
    [~,donorColumns]=find(coreFromPolar);
    donorRows=polar.activeRows(mod(donorColumns-1,numel(polar.activeRows))+1);
    assert(all(selected(donorRows)<=redistributionInfo.preservedPrefixCount)&& ...
        selected(halo)<=redistributionInfo.preservedPrefixCount&& ...
        selected(halo+1)<=redistributionInfo.copiedMetricPrefixCount, ...
        'ipm_cr:RadialRedistributionProtection', ...
        'The protected radius must retain complete overlap donor stencils and the inner active metric.');
end
assert(size(coreFromPolar,2)==nnz(polarActive)&&size(polarFromCore,2)==nnz(active), ...
    'ipm_cr:ExplicitOverlap','The exchange map must read the declared active donor packing.');
exchange=struct('coreFromPolar',coreFromPolar,'polarFromCore',polarFromCore, ...
    'coreDifferencePlan',difference_plan(coreFromPolar), ...
    'polarDifferencePlan',difference_plan(polarFromCore), ...
    'coreGhostIndices',core.ghostIndices,'polarGhostIndices',polar.ghostIndices, ...
    'coreDonorIndices',core.activeIndices,'polarDonorIndices',polar.activeIndices, ...
    'coreMapInfinityNorm',norm(coreFromPolar,inf), ...
    'polarMapInfinityNorm',norm(polarFromCore,inf), ...
    'interpolationWidth',width,'maximumSupportNodes',width^2, ...
    'donorsAreActiveOnly',true,'simultaneousExchange',true);
g=struct('kind','overlapping_uniform_core_mapped_annulus_explicit_transport_v1', ...
    'core',core,'polar',polar,'exchange',exchange,'gridOptions',grid, ...
    'ghostWidth',halo,'interpolationWidth',width, ...
    'activeDofCount',nnz(active)+nnz(polarActive), ...
    'activeOverlapEvolvedTwice',true,'polarOriginPresent',false, ...
    'buildSeconds',toc(timer));
g.geometryKey=struct('kind',g.kind,'x',x,'y',y,'r',r,'theta',theta, ...
    'drMetric',drMetric,'ghostWidth',halo,'interpolationWidth',width, ...
    'physicalDerivativeClosure','same_patch_active_D7_v1');
if mapped
    g.kind='overlapping_mapped_tensor_annulus_explicit_transport_v1';
    g.geometryKey.kind=g.kind;g.geometryKey.coreMap=core.coreMap;
    g.geometryKey.Jx=core.Jx;g.geometryKey.Jy=core.Jy;
    g.geometryKey.dxMetric=core.dxMetric;g.geometryKey.dyMetric=core.dyMetric;
    g.geometryKey.logicalSpacing=[core.dxi,core.deta];
    g.geometryKey.mapParameters=core.mapParameters;
end
if redistribute
    redistributionInfo.annularFirstFullIndex=first;
    redistributionInfo.retainedNr=nr;
    redistributionInfo.overlapDonorsProtected=true;
    g.radialRedistribution=redistributionInfo;
end
if explicit,g.axisSource=axisSource;end
end

function D=active_derivative(axis,active)
% Physical q and q_tau use the same independent patch degrees of freedom.
% Advective stencils still use the exchanged collars. Differentiating a
% synchronized RHS across those collars would mix the two patches' distinct
% truncation errors; use degree-six one-sided closures at artificial edges.
D=ipm_cr.mesh.fd(axis,1,7);
D(active,:)=0;
D(active,active)=ipm_cr.mesh.fd(axis(active),1,7);
end

function [a,source]=explicit_axes(o)
% One physical-axis boundary for portable presets and remesh candidates.
assert(~(isfield(o,'axes')&&isfield(o,'preset')),'ipm_cr:ExplicitAxes', ...
    'Use either preset or axes, not both.');
assert(~any(isfield(o,{'coreHalfWidth','coreHeight','coreSpacing','polarStart', ...
    'coreMap','radialRedistribution'})),'ipm_cr:ExplicitAxes', ...
    'Physical axes cannot be combined with coordinate-generator options.');
names={'x','y','r','theta','dxMetric','dyMetric','drMetric'};
if isfield(o,'preset')
    name=char(validatestring(o.preset,{'medium','coarse'}));
    data=jsondecode(fileread(fullfile(ipm_cr.root(),'data','mesh.json')));
    assert(strcmp(data.format,'ipm_cr_physical_axes_v1'),'ipm_cr:MeshPreset', ...
        'Unsupported portable mesh data format.');
    a=data.axes;stride=1+strcmp(name,'coarse');
    for k=1:4
        key=names{k};assert(mod(numel(a.(key))-1,stride)==0,'ipm_cr:MeshPreset', ...
            'Exact halving requires even interval counts in every direction.');
        a.(key)=a.(key)(1:stride:end);
    end
    for k=5:7,key=names{k};a.(key)=stride*a.(key)(1:stride:end);end
    source=struct('kind','portable_preset','preset',name, ...
        'stride',stride,'sourceGeometrySHA256',data.sourceGeometrySHA256, ...
        'coordinatesAreExactSubsequence',true);
else
    a=o.axes;
    source=struct('kind','explicit_axes','coordinatesAreExactSubsequence',false);
end
assert(isstruct(a)&&isscalar(a)&&all(isfield(a,names))&& ...
    all(ismember(fieldnames(a),names)),'ipm_cr:ExplicitAxes', ...
    'Supply exactly x/y/r/theta and dxMetric/dyMetric/drMetric.');
for k=1:4
    key=names{k};validateattributes(a.(key),{'double'},{'vector','real','finite','increasing'});
    assert(numel(a.(key))>=7,'ipm_cr:ExplicitAxes','Each physical axis must fit a seven-point stencil.');
end
a.x=reshape(a.x,1,[]);a.y=a.y(:);a.r=a.r(:);a.theta=reshape(a.theta,1,[]);
assert(a.x(1)==-a.x(end)&&any(a.x==0)&&a.y(1)==0&&a.y(end)>0&&a.r(1)>0, ...
    'ipm_cr:ExplicitAxes','Use a symmetric Cartesian box with exact origin and positive annular radii.');
assert(a.theta(1)==0&&a.theta(end)==pi&& ...
    max(abs(diff(a.theta)-pi/(numel(a.theta)-1)))<=64*eps(pi), ...
    'ipm_cr:ExplicitAxes','The angular axis must be uniform from zero to pi.');
counts=[numel(a.x),numel(a.y),numel(a.r)];
for k=1:3
    key=names{k+4};validateattributes(a.(key),{'double'}, ...
        {'vector','real','finite','positive','numel',counts(k)});
end
a.dxMetric=reshape(a.dxMetric,1,[]);a.dyMetric=a.dyMetric(:);a.drMetric=a.drMetric(:);
if isfield(o,'grid')
    assert(isstruct(o.grid)&&isscalar(o.grid)&&all(ismember(fieldnames(o.grid),{'rmax'})), ...
        'ipm_cr:ExplicitAxes','Only the coupled rmax check is allowed with physical axes.');
    if isfield(o.grid,'rmax')
        assert(o.grid.rmax==a.r(end),'ipm_cr:DomainRadius', ...
            'The explicit outer axis must match the requested physical radius.');
    end
end
a.dxi=2/(counts(1)-1);a.deta=1/(counts(2)-1);
a.xi=linspace(-1,1,counts(1));a.eta=linspace(0,1,counts(2))';
a.Jx=a.dxMetric/a.dxi;a.Jy=a.dyMetric/a.deta;
% Keep the same multiplication identity used by transport validation.
a.dxMetric=a.Jx*a.dxi;a.dyMetric=a.Jy*a.deta;
a.coreMap=struct('enabled',true,'kind','explicit_axes','intervals',counts(1:2)-1, ...
    'minSpacing',[min(a.dxMetric),min(a.dyMetric)]);
a.mapParameters=struct('originMetricSpacing', ...
    [a.dxMetric(find(a.x==0,1)),a.dyMetric(1)], ...
    'averagePhysicalSpacing',[2*a.x(end)/(counts(1)-1),a.y(end)/(counts(2)-1)], ...
    'logicalSpacing',[a.dxi,a.deta],'logicalBounds',[-1,1,0,1], ...
    'definition','supplied_positive_metric_step_equals_J_times_logical_spacing');
end

function value=option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end

function p=difference_plan(M)
% Evaluate each transfer row around its strongest local donor. Constants
% then produce exactly zero differences, independently of row-sum rounding.
[row,column,weight]=find(M);[~,anchor]=max(abs(M),[],2);
p=struct('row',row,'column',column,'anchor',full(anchor), ...
    'sum',sparse(row,(1:numel(row))',weight,size(M,1),numel(row)));
end
