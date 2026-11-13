function [state,context,info]=remesh(state,context,settings,force)
%REMESH Shape-triggered mapped core and annulus reconstruction.
% Rebuild the mapped patches, transfer active fields and evaluate one complete RHS.
% Calibrate the wall-connected main 25% and 80% q regions on the first check.
% A 25% change from either last-remesh anchor triggers mesh adaptation.
% P=remesh('config') returns the serializable controller defaults.
% R=remesh('region',S,G,P) measures the levels in P, or 25% when omitted.
% G=remesh('geometry',G0,EXTENTS0,EXTENTS) builds geometry without state transfer.
% [T,A]=remesh('transfer',S,GOLD,GNEW) uses active-only degree-six transfer.
% Fixed-frame redistribution continues each patch's own endpoint polynomial,
% preserving its independent overlap representation as active boundaries move;
% T has no valid context or RHS and must never be advanced before rebinding.
if ischar(state)||(isstring(state)&&isscalar(state))
    switch char(state)
        case 'config'
            if nargin<2,context=struct();end
            state=policy(context);context=[];info=[];
        case 'region'
            if nargin<4,force=struct('level',.25);end
            state=region(context,settings,policy(force));context=[];info=[];
        case 'geometry'
            memory=struct('initialMesh',mesh_calibration(context),'initialExtents',settings);
            [axes,info]=propose(force,memory);
            state=rebuild_geometry(context,axes);context=[];
        case 'transfer',[state,context]=transfer(context,settings,force);info=[];
        otherwise,error('ipm_cr:RemeshAction','Use config, region, geometry, or transfer.');
    end
    return
end
if nargin<3
    settings=struct();
    if isa(context,'ipm_cr.context')&&isfield(context.opts,'remesh'),settings=context.opts.remesh;end
end
if nargin<4,force=false;end
p=policy(settings);
tau=state.scale.canonicalTime;g=context.mesh;
memory=controller(state);
nLevels=numel(p.level);
info=struct('status','not_triggered','reason','disabled','remeshed',false, ...
    'timeAdvanced',0,'forced',force,'tau',tau,'step',state.step, ...
    'checked',false,'triggered',false,'directionalChange',zeros(nLevels,2), ...
    'triggerMask',false(nLevels,2),'directionalRateEstimate',nan(nLevels,2), ...
    'checkCadence',p.checkEveryTau, ...
    'independentGreenAuditPerformed',false,'continuumAccuracyCertified',false);
if ~p.enabled&&~force,return;end
if ~isempty(memory.referenceShape)&&~force&&tau<memory.nextCheckTau
    info.reason='check_cadence';return
end
shape=region(state,g,p);
require_regions(shape);
extents=vertcat(shape.triggerExtents);
info.checked=true;info.region=shape;memory.checks=memory.checks+1;
if isempty(memory.referenceShape)
    memory.referenceShape=small_shape(shape);
    memory.initialExtents=shape(1).extents;
    memory.initialMesh=mesh_calibration(g);
    memory.lastShape=memory.referenceShape;memory.lastCheckTau=tau;
    memory.nextCheckTau=tau+p.checkEveryTau;
    if ~force,info.reason='initial_calibration';finish_controller();return;end
end
info.referenceShape=memory.referenceShape;
reference=vertcat(memory.referenceShape.triggerExtents);
info.directionalChange=abs(extents./reference-1);
info.triggerMask=extents>=reference*(1+p.threshold)|extents<=reference*(1-p.threshold);
interval=tau-memory.lastCheckTau;cadence=p.checkEveryTau;
if interval>0
    info.directionalRateEstimate=abs(extents-vertcat(memory.lastShape.triggerExtents))./ ...
        reference/interval;
    speed=max(info.directionalRateEstimate,[],'all');
    if speed>0
        margin=max(0,p.threshold-max(info.directionalChange,[],'all'));
        cadence=max(p.minCheckEveryTau,min(p.checkEveryTau,.8*margin/speed));
    end
end
memory.lastCheckTau=tau;memory.lastShape=small_shape(shape);
memory.nextCheckTau=tau+cadence;info.checkCadence=cadence;
info.triggered=force||any(info.triggerMask,'all');
if ~info.triggered
    info.reason='below_directional_threshold';finish_controller();return
end
timer=tic;
if isempty(p.candidateAxes)
    [axes,info.geometryProposal]=propose(shape(1).extents,memory);
else
    axes=p.candidateAxes;
end
staircase=isfield(g.core,'stair');
if ~staircase&&isequal(axes.x(:),g.core.x(:))&&isequal(axes.y(:),g.core.y(:))&& ...
        isequal(axes.r(:),g.polar.r(:))&&isequal(axes.theta(:),g.polar.theta(:))
    info.reason='geometry_unchanged';finish_controller();return
end
newG=rebuild_geometry(g,axes);
if staircase
    if isequaln(newG.geometryKey,g.geometryKey)
        info.reason='geometry_unchanged';finish_controller();return
    end
    info.geometryProposal.innerRadius=newG.polar.r(1);
    info.geometryProposal.activeRadialLayers=newG.core.stair.activeRadialLayers;
end
[moved,info.transfer]=transfer(state,g,newG,false);
% Actual geometry changes; the original method prescription does not.
newC=ipm_cr.context(newG,context.wall,context.opts);
moved.ctxId=newC.id;moved.ctxKey=newC.key;
[moved.rhs,moved.W]=ipm_cr.rhs(moved,newC);
memory.remeshCount=memory.remeshCount+1;
memory.lastRemeshTau=tau;
afterShape=region(moved,newG,p);
require_regions(afterShape);
memory.referenceShape=small_shape(afterShape);
memory.lastShape=memory.referenceShape;memory.nextCheckTau=tau+p.checkEveryTau;
info.checkCadence=p.checkEveryTau;
info.status='remeshed';info.reason='shape_change';info.remeshed=true;
state=moved;context=newC;
info.totalSeconds=toc(timer);finish_controller();

    function finish_controller()
        state.remesh=memory;info.controller=memory;
    end
end

function new=rebuild_geometry(old,axes)
geometry=struct('axes',axes,'ghostWidth',old.ghostWidth,'interpolationWidth',old.interpolationWidth);
if isfield(old.core,'stair')
    geometry.staircase=true;
    geometry.staircaseFine=old.axisSource.staircaseFine;
    if isfield(old.axisSource,'interfaceAllocation')
        geometry.interfaceAllocation=old.axisSource.interfaceAllocation.measure.plan;
    end
end
new=ipm_cr.mesh.build(geometry);
end

function p=policy(overrides)
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:RemeshOptions','Use scalar remesh settings.');
p=struct('enabled',true,'level',[.25,.80],'threshold',.25,'checkEveryTau',.5, ...
    'minCheckEveryTau',.005,'candidateAxes',[]);
names=fieldnames(overrides);
assert(all(ismember(names,fieldnames(p))),'ipm_cr:RemeshOptions','Unknown remesh option.');
for k=1:numel(names),p.(names{k})=overrides.(names{k});end
validateattributes(p.enabled,{'logical'},{'scalar'});
skip={'enabled','candidateAxes','level'};
allNames=fieldnames(p);
for k=1:numel(allNames)
    key=allNames{k};if ismember(key,skip),continue;end
    validateattributes(p.(key),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(p.level,{'numeric'},{'vector','real','finite','>',0,'<',1,'increasing'});
assert(numel(p.level)<=2,'ipm_cr:RemeshOptions','Use one or two increasing superlevel fractions.');
p.level=reshape(p.level,1,[]);
assert(p.minCheckEveryTau<=p.checkEveryTau,'ipm_cr:RemeshOptions', ...
    'minCheckEveryTau must not exceed checkEveryTau.');
assert(isempty(p.candidateAxes)||(isstruct(p.candidateAxes)&&isscalar(p.candidateAxes)), ...
    'ipm_cr:RemeshOptions','candidateAxes must be empty or an explicit-axis structure.');
end

function c=controller(s)
if isfield(s,'remesh')
    c=s.remesh;
    return
end
c=struct('format','ipm_cr_mapped_remesh_controller_v2','lastRemeshTau',-realmax, ...
    'lastCheckTau',-realmax,'nextCheckTau',-realmax,'lastShape',[], ...
    'referenceShape',[],'initialExtents',[],'checks',0,'remeshCount',0, ...
    'initialMesh',[]);
end

function initial=mesh_calibration(g)
c=g.core;grid=g.gridOptions;
stride=1;factor=1;
if isfield(g,'axisSource')&&strcmp(g.axisSource.kind,'portable_preset')
    % Full radial map used to construct the reference mesh in the paper.
    stride=g.axisSource.stride;factor=2;
    grid=struct('nr',385,'ntheta',numel(g.polar.theta),'rmax',g.polar.r(end), ...
        'firstSpacing',.03,'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01);
elseif strcmp(grid.radialMap,'explicit_axes')
    grid.radialMap='uniform_core';
    grid.firstSpacing=min(min(g.polar.drMetric),grid.rmax/(4*(grid.nr-1)));
    grid.coreRadius=min(max(2,2*grid.firstSpacing),grid.rmax/4);
elseif isfield(g,'radialRedistribution')
    factor=g.radialRedistribution.factor;
end
if isfield(c,'dxMetric')
    spacing=[c.dxMetric(c.originIndex),c.dyMetric(1)];
else
    spacing=[c.dx,c.dy];
end
initial=struct('bounds',[c.x(end),c.y(end)], ...
    'coreMap',struct('intervals',[c.nx-1,c.ny-1],'minSpacing',spacing), ...
    'theta',g.polar.theta,'radialGrid',grid,'radialStride',stride,'tailFactor',factor);
if isfield(c,'stair')
    initial.stairAxes=struct('x',c.x,'y',c.y,'dxMetric',c.dxMetric,'dyMetric',c.dyMetric);
    if isfield(g.axisSource,'interfaceAllocation')
        allocation=g.axisSource.interfaceAllocation;
        initial.stairAxes.x=allocation.referenceX;
        initial.stairAxes.dxMetric=allocation.referenceDx;
        if isfield(allocation,'referenceY')
            initial.stairAxes.y=allocation.referenceY;
            initial.stairAxes.dyMetric=allocation.referenceDy;
        end
    end
end
end

function s=small_shape(s)
s=rmfield(s,'crossings');
end

function require_regions(shapes)
for shape=shapes
    assert(shape.valid,'ipm_cr:RemeshRegion', ...
        'The %.0f%% wall-connected level set is undefined: %s.',100*shape.level,shape.reason);
end
end

function shapes=region(s,g,p)
% Piecewise-linear superlevel geometry of the original active-owner field.
% Both levels share one triangulation covering the lowest superlevel set.
assert(isstruct(s)&&isfield(s,'rhs')&&isfield(s.rhs,'q')&& ...
    isequal(size(s.rhs.q.core),size(g.core.X))&& ...
    isequal(size(s.rhs.q.polar),size(g.polar.X)), ...
    'ipm_cr:RemeshRegion','Use native q arrays matching both patches.');
q=s.rhs.q;ca=g.core.active;pa=g.polar.active;
cx=g.core.X(ca);cy=g.core.Y(ca);px=g.polar.X(pa);py=g.polar.Y(pa);
qc=q.core(ca);qp=q.polar(pa);[pc,ic]=max(qc);[pp,ip]=max(qp);
shape=struct('valid',false,'reason','nonpositive_peak','level',NaN, ...
    'definition','positive_peak_connected_component_projection_v1', ...
    'method','active_owner_local_linear_triangles_subcell_crossings', ...
    'peak',max(pc,pp),'peakXY',[NaN,NaN],'ownedPeak',NaN,'ownedPeakXY',[NaN,NaN], ...
    'peakDifference',NaN,'peakRelativeDifference',NaN,'threshold',NaN, ...
    'xBounds',[NaN,NaN],'yBounds',[NaN,NaN],'widths',[NaN,NaN],'extents',[NaN,NaN], ...
    'triggerExtents',[NaN,NaN], ...
    'monitorRadius',NaN,'monitorNodes',0,'triangleCount',0,'componentNodes',0, ...
    'otherPositiveNodes',0,'crossings',zeros(0,2),'wallConnected',false, ...
    'globalActiveCount',numel(qc)+numel(qp),'ghostDonorsUsed',false, ...
    'overlapOwner','core_inside_closed_active_rectangle');
if pc>=pp,shape.peakXY=[cx(ic),cy(ic)];else,shape.peakXY=[px(ip),py(ip)];end
shapes=fill_levels(shape,p.level);
if ~(shape.peak>0),return;end
threshold=p.level(1)*shape.peak;
radius=max(hypot(cx,cy));above=qp>=threshold;
if any(above),radius=max(radius,max(hypot(px(above),py(above))));end
rr=g.polar.r(g.polar.activeRows);
outer=find(rr>radius+64*eps(max(1,radius)),1);
if isempty(outer)
    shape.reason='level_set_not_enclosed_by_active_domain';
    shapes=fill_levels(shape,p.level);return
end
shape.monitorRadius=rr(outer);
inside=px>=min(cx)&px<=max(cx)&py>=min(cy)&py<=max(cy);
if isfield(g.core,'stair')
    inside=ipm_cr.mesh.staircase('contains',g.core.stair,px,py);
    shape.overlapOwner='core_inside_closed_active_staircase';
end
take=~inside&hypot(px,py)<=shape.monitorRadius+64*eps(max(1,shape.monitorRadius));
points=[cx,cy;px(take),py(take)];values=[qc;qp(take)];
[points,keep]=unique(points,'rows','stable');values=values(keep);
shape.monitorNodes=size(points,1);[shape.ownedPeak,seed]=max(values);
shape.ownedPeakXY=points(seed,:);shape.peakDifference=shape.peak-shape.ownedPeak;
shape.peakRelativeDifference=shape.peakDifference/shape.peak;
triangles=delaunay(points(:,1),points(:,2));shape.triangleCount=size(triangles,1);
[edges,~,edgeId]=unique(sort([triangles(:,[1,2]);triangles(:,[2,3]);triangles(:,[3,1])],2),'rows');
triangleEdges=reshape(edgeId,[],3);
shapes=fill_levels(shape,p.level);
for k=1:numel(shapes)
    shapes(k)=measure_level(shapes(k),points,values,seed,edges,triangleEdges);
end
end

function shapes=fill_levels(shape,levels)
shapes=repmat(shape,size(levels));
for k=1:numel(shapes)
    shapes(k).level=levels(k);
    shapes(k).threshold=levels(k)*shape.peak;
end
end

function shape=measure_level(shape,points,values,seed,edges,triangleEdges)
threshold=shape.threshold;
if shape.ownedPeak<threshold,shape.reason='peak_absent_from_owner_field';return;end
positive=values>=threshold;joined=edges(positive(edges(:,1))&positive(edges(:,2)),:);
n=size(points,1);adjacency=sparse([joined(:,1);joined(:,2)],[joined(:,2);joined(:,1)],1,n,n);
component=false(n,1);component(seed)=true;queue=zeros(n,1);queue(1)=seed;head=1;tail=1;
while head<=tail
    neighbors=find(adjacency(:,queue(head)));head=head+1;
    fresh=neighbors(~component(neighbors));component(fresh)=true;
    queue(tail+(1:numel(fresh)))=fresh;tail=tail+numel(fresh);
end
shape.componentNodes=nnz(component);shape.otherPositiveNodes=nnz(positive&~component);
shape.wallConnected=any(component&points(:,2)==0);
if ~shape.wallConnected,shape.reason='main_component_detached_from_wall';return;end
isCrossed=xor(component(edges(:,1)),component(edges(:,2)));
crossed=edges(isCrossed,:);
a=crossed(:,1);b=crossed(:,2);
fraction=(threshold-values(a))./(values(b)-values(a));
crossings=nan(size(edges,1),2);
crossings(isCrossed,:)=points(a,:)+fraction.*(points(b,:)-points(a,:));
shape.crossings=unique(crossings(isCrossed,:),'rows','stable');
support=[shape.crossings;points(component,:)];
shape.xBounds=[min(support(:,1)),max(support(:,1))];
shape.yBounds=[min(support(:,2)),max(support(:,2))];
shape.widths=[diff(shape.xBounds),diff(shape.yBounds)];
shape.extents=[max(abs(shape.xBounds)),shape.yBounds(2)];
height=peak_section(crossings,triangleEdges,isCrossed,shape.ownedPeakXY);
if isempty(height),shape.reason='peak_normal_section_not_enclosed';return;end
shape.triggerExtents=[shape.extents(1),height];
shape.valid=true;shape.reason='ok';
end

function height=peak_section(crossings,triangleEdges,isCrossed,peak)
% Every mixed triangle cuts the contour in a segment. Intersect that
% segment with the wall-normal ray and take the first exit above the peak.
cut=triangleEdges(any(isCrossed(triangleEdges),2),:);
cut(~isCrossed(cut))=0;
cut=sort(cut,2,'descend');
a=crossings(cut(:,1),:);b=crossings(cut(:,2),:);
x=peak(1);tol=16*eps(max(1,abs(x)));
between=x>=min(a(:,1),b(:,1))-tol&x<=max(a(:,1),b(:,1))+tol;
dx=b(:,1)-a(:,1);slanted=between&dx~=0;vertical=between&dx==0;
t=max(0,min(1,(x-a(slanted,1))./dx(slanted)));
y=[a(slanted,2)+t.*(b(slanted,2)-a(slanted,2));a(vertical,2);b(vertical,2)];
height=min(y(y>peak(2)));
end

function [axes,info]=propose(extents,memory)
initial=memory.initialMesh;
ratio=extents./memory.initialExtents;
bounds=initial.bounds.*ratio;
if isfield(initial,'stairAxes')
    core=initial.stairAxes;
    core.x=ratio(1)*core.x;core.dxMetric=ratio(1)*core.dxMetric;
    core.y=ratio(2)*core.y;core.dyMetric=ratio(2)*core.dyMetric;
else
    core=ipm_cr.mesh.tensorAxis(bounds(1),bounds(2),initial.coreMap);
end

grid=initial.radialGrid;
innerScale=min([ratio,1]);
grid.firstSpacing=grid.firstSpacing*innerScale;
if isfield(grid,'coreRadius'),grid.coreRadius=grid.coreRadius*innerScale;end
full=ipm_cr.mesh.polar(grid);
r=full.r;dr=full.radialMetric/(full.nr-1);
if initial.tailFactor>1
    preserve=max(4*innerScale,2.5*hypot(bounds(1),bounds(2)));
    transition=min(grid.rmax,max(20*innerScale,4*preserve));
    [r,dr]=ipm_cr.mesh.radial(full,struct('enabled',true,'factor',initial.tailFactor, ...
        'preserveRadius',preserve,'fullRefinementRadius',transition));
end
% Keep the common outer endpoint when halving the regenerated full axis.
stride=initial.radialStride;keep=fliplr(numel(r):-stride:1);
r=r(keep);dr=stride*dr(keep);
requested=.6*min(bounds);
first=find(r>=requested,1);
if isfield(initial,'stairAxes')
    % Let the geometry constructor align the active staircase and annulus.
    first=find(r>0,1);requested=NaN;
end
axes=struct('x',core.x,'y',core.y,'r',r(first:end),'theta',initial.theta, ...
    'dxMetric',core.dxMetric,'dyMetric',core.dyMetric,'drMetric',dr(first:end));
info=struct('initialBounds',initial.bounds,'initialExtents',memory.initialExtents, ...
    'measuredExtents',extents,'bounds',bounds,'requestedInnerRadius',requested, ...
    'innerRadius',axes.r(1),'precedingRadius',r(max(1,first-1)), ...
    'outerRadius',r(end),'radialStride',stride,'fullRadialNodes',numel(r));
end

function [out,audit]=transfer(s,old,new,allowContinuation)
if nargin<4,allowContinuation=true;end
out=s;key='W';if isfield(s,'eta'),key='eta';end
[out.(key),audit]=transfer_field(s.(key),old,new,allowContinuation);
if strcmp(key,'eta'),out.W=struct('core',zeros(size(new.core.X)),'polar',zeros(size(new.polar.X)));end
out.ctxId='';out.ctxKey=[];
out.rhs=struct('c_l',s.rhs.c_l,'c_r',s.rhs.c_r);
audit.integratedField=key;audit.timeAdvanced=0;audit.originalReferencesPreserved=true;
end

function [out,audit]=transfer_field(input,old,new,allowContinuation)
out=struct();counts=struct();maxNorm=0;
staircase=isfield(old.core,'stair');
for patch={'core','polar'}
    key=patch{1};target=new.(key);indices=target.activeIndices;
    x=target.X(indices);y=target.Y(indices);
    if strcmp(key,'polar')
        % Native axes avoid hypot roundoff outside the outer circle.
        [ir,it]=ind2sub(size(target.X),indices);
        r=target.r(ir);theta=target.theta(it);theta=theta(:);
    else
        r=hypot(x,y);theta=atan2(y,x);
    end
    coreCovered=covered(x,old.core.activeXBounds)&covered(y,old.core.activeYBounds);
    if staircase,coreCovered=ipm_cr.mesh.staircase('contains',old.core.stair,x,y);end
    polarCovered=covered(r,[old.polar.innerActiveRadius,old.polar.r(end)])& ...
        covered(theta,[old.polar.theta(1),old.polar.theta(end)]);
    fixedFrame=same_frame(old.(key),target,key);
    continueLocal=allowContinuation&&fixedFrame&&(~staircase||strcmp(key,'polar'));
    if strcmp(key,'core'),coveredBySame=coreCovered;else,coveredBySame=polarCovered;end
    same=coveredBySame;
    if continueLocal,same(:)=true;end
    if strcmp(key,'core'),coreOwner=same;else,coreOwner=~same;end
    % Standalone transfers retain fixed-frame continuation when requested.
    % Mapped remeshing uses the same patch only within its active coverage.
    values=zeros(size(x));
    for owner={'core','polar'}
        name=owner{1};take=coreOwner;if strcmp(name,'polar'),take=~take;end
        if ~any(take),continue;end
        p=old.(name);rows=p.activeRows;columns=p.activeColumns;
        if strcmp(name,'core')&&staircase
            M=ipm_cr.mesh.staircase('donors',p,x(take),y(take),7);
            donor=input.core(p.activeIndices);
        elseif strcmp(name,'core')
            M=ipm_cr.mesh.tensorMap(p.y(rows),p.x(columns),y(take),x(take),7,continueLocal);
            donor=input.core(rows,columns);
        else
            M=ipm_cr.mesh.tensorMap(p.r(rows),p.theta(columns),r(take),theta(take),7,continueLocal);
            donor=input.polar(rows,columns);
        end
        values(take)=difference_apply(M,donor(:));maxNorm=max(maxNorm,norm(M,inf));
    end
    out.(key)=zeros(size(target.X));out.(key)(indices)=values;
    counts.(key)=struct('samePatch',nnz(same),'crossPatch',nnz(~same), ...
        'samePatchBoundaryExtension',nnz(same&~coveredBySame),'fixedPatchFrame',fixedFrame);
end
out=ipm_cr.mesh.fill(out,new);
audit=struct('kind','independent_active_patch_local7_same_time_transfer_v1', ...
    'routes',counts,'mapInfinityNorm',maxNorm,'ghostValuesRead',false, ...
    'activeOverlapProjected',false,'polynomialDegree',6, ...
    'fixedFrameContinuationAllowed',allowContinuation, ...
    'continuumAccuracyCertified',false);
end

function yes=same_frame(a,b,patch)
% Complete axes keep the patch frame; active endpoints may move with their
% third interior nodes. Counts and full endpoints identify redistribution,
% as distinct from cross-method import or patch-domain replacement.
if strcmp(patch,'core'),first='y';second='x';else,first='r';second='theta';end
yes=isequal(size(a.X),size(b.X))&& ...
    isequal(a.(first)([1,end]),b.(first)([1,end]))&& ...
    isequal(a.(second)([1,end]),b.(second)([1,end]));
end

function yes=covered(value,bounds)
% Use each endpoint's local ulp, as in localMap's roundoff repair.
yes=value>=bounds(1)-64*eps(max(1,abs(bounds(1))))& ...
    value<=bounds(2)+64*eps(max(1,abs(bounds(2))));
end

function v=difference_apply(M,donor)
[row,column,weight]=find(M);[~,anchor]=max(abs(M),[],2);anchor=full(anchor);
base=donor(anchor);v=base+accumarray(row,weight.*(donor(column)-base(row)),[size(M,1),1]);
end
