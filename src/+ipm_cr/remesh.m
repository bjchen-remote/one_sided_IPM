function [state,context,info]=remesh(state,context,settings,force)
%REMESH Shape-triggered source-driven same-time core redistribution.
% One fixed-reference redistribution, active transfer and complete RHS.
% First calibrate the wall-connected main 25% q region without interpolation.
% Later checks adapt only after a 20% directional change from the last remesh.
% P=remesh('config') returns the serializable controller defaults.
% M=remesh('monitor',S,G) inspects active core resolution, without a solve.
% R=remesh('region',S,G,P) measures the main wall-connected q superlevel set.
% [T,A]=remesh('transfer',S,GOLD,GNEW) uses active-only degree-six transfer.
% Fixed-frame redistribution continues each patch's own endpoint polynomial,
% preserving its independent overlap representation as active boundaries move;
% T has no valid context or RHS and must never be advanced before rebinding.
if ischar(state)||(isstring(state)&&isscalar(state))
    switch char(state)
        case 'config'
            if nargin<2,context=struct();end
            state=policy(context);context=[];info=[];
        case 'monitor',state=monitor(context,settings);context=[];info=[];
        case 'region'
            if nargin<4,force=struct();end
            state=region(context,settings,policy(force));context=[];info=[];
        case 'transfer',[state,context]=transfer(context,settings,force);info=[];
        otherwise,error('ipm_cr:RemeshAction','Use config, monitor, region, or transfer.');
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
memory=controller(state,g,p);
info=struct('status','not_triggered','reason','disabled','remeshed',false, ...
    'timeAdvanced',0,'forced',force,'tau',tau,'step',state.step, ...
    'checked',false,'triggered',false,'directionalChange',[0,0], ...
    'directionalRateEstimate',[NaN,NaN],'checkCadence',p.checkEveryTau, ...
    'independentGreenAuditPerformed',false,'continuumAccuracyCertified',false);
if ~p.enabled&&~force,return;end
if ~isempty(memory.referenceShape)&&~force&&tau<memory.nextCheckTau
    info.reason='check_cadence';return
end
shape=region(state,g,p);
assert(shape.valid,'ipm_cr:RemeshRegion','The main wall-connected region is undefined: %s.',shape.reason);
info.checked=true;info.region=shape;memory.checks=memory.checks+1;
if isempty(memory.referenceShape)
    memory.referenceShape=small_shape(shape);
    memory.lastShape=memory.referenceShape;memory.lastCheckTau=tau;
    memory.nextCheckTau=tau+p.checkEveryTau;
    if ~force,info.reason='initial_calibration';finish_controller();return;end
end
info.referenceShape=memory.referenceShape;
info.directionalChange=abs(shape.extents./memory.referenceShape.extents-1);
interval=tau-memory.lastCheckTau;cadence=p.checkEveryTau;
if interval>0
    info.directionalRateEstimate=abs(shape.extents-memory.lastShape.extents)./ ...
        memory.referenceShape.extents/interval;
    speed=max(info.directionalRateEstimate);
    if speed>0
        margin=max(0,p.threshold-max(info.directionalChange));
        cadence=max(p.minCheckEveryTau,min(p.checkEveryTau,.8*margin/speed));
    end
end
memory.lastCheckTau=tau;memory.lastShape=small_shape(shape);
memory.nextCheckTau=tau+cadence;info.checkCadence=cadence;
info.triggered=force||max(info.directionalChange)>=p.threshold;
if ~info.triggered
    info.reason='below_directional_threshold';finish_controller();return
end
info.monitor=monitor(state,g);
timer=tic;
if isempty(p.candidateAxes),axes=propose(g,info.monitor,p,memory.referenceCore);
else,axes=p.candidateAxes;end
if isequal(axes.x(:),g.core.x(:))&&isequal(axes.y(:),g.core.y(:))&& ...
        isequal(axes.r(:),g.polar.r(:))&&isequal(axes.theta(:),g.polar.theta(:))
    info.reason='geometry_unchanged';finish_controller();return
end
newG=ipm_cr.mesh.build(struct('axes',axes, ...
    'ghostWidth',g.ghostWidth,'interpolationWidth',g.interpolationWidth));
[moved,info.transfer]=transfer(state,g,newG);
% Actual geometry changes; the original method prescription does not.
newC=ipm_cr.context(newG,context.wall,context.opts);
moved.ctxId=newC.id;moved.ctxKey=newC.key;
[moved.rhs,moved.W]=ipm_cr.rhs(moved,newC);
info.monitorAfter=monitor(moved,newG);
memory.remeshCount=memory.remeshCount+1;
memory.lastRemeshTau=tau;
afterShape=region(moved,newG,p);
assert(afterShape.valid,'ipm_cr:RemeshRegion','The new main wall-connected region is undefined: %s.',afterShape.reason);
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

function p=policy(overrides)
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:RemeshOptions','Use scalar remesh settings.');
p=struct('enabled',true,'level',.25,'threshold',.2,'checkEveryTau',.5, ...
    'minCheckEveryTau',.005,'monitorStrength',2,'blend',.25,'maxNodeMotion',.3,'candidateAxes',[]);
names=fieldnames(overrides);
assert(all(ismember(names,fieldnames(p))),'ipm_cr:RemeshOptions','Unknown remesh option.');
for k=1:numel(names),p.(names{k})=overrides.(names{k});end
validateattributes(p.enabled,{'logical'},{'scalar'});
skip={'enabled','candidateAxes'};
allNames=fieldnames(p);
for k=1:numel(allNames)
    key=allNames{k};if ismember(key,skip),continue;end
    validateattributes(p.(key),{'numeric'},{'scalar','real','finite','positive'});
end
assert(p.blend<=1,'ipm_cr:RemeshOptions','A convex redistribution requires blend<=1.');
assert(p.level<1,'ipm_cr:RemeshOptions','The superlevel fraction must be less than one.');
assert(p.minCheckEveryTau<=p.checkEveryTau,'ipm_cr:RemeshOptions', ...
    'minCheckEveryTau must not exceed checkEveryTau.');
assert(isempty(p.candidateAxes)||(isstruct(p.candidateAxes)&&isscalar(p.candidateAxes)), ...
    'ipm_cr:RemeshOptions','candidateAxes must be empty or an explicit-axis structure.');
end

function c=controller(s,g,p)
if isfield(s,'remesh')
    c=s.remesh;
    if ~isfield(c,'referenceShape'),c.referenceShape=[];c.lastCheckTau=-realmax;end
    if ~isfield(c,'checks'),c.checks=0;end
    if ~isfield(c,'lastShape'),c.lastShape=c.referenceShape;end
    if ~isfield(c,'nextCheckTau'),c.nextCheckTau=c.lastCheckTau+p.checkEveryTau;end
else
    c=struct('format','ipm_cr_remesh_controller_v5','lastRemeshTau',-realmax, ...
        'lastCheckTau',-realmax,'nextCheckTau',-realmax,'lastShape',[], ...
        'referenceShape',[],'checks',0,'remeshCount',0, ...
        'referenceCore',struct('x',g.core.x,'y',g.core.y));
end
end

function s=small_shape(s)
s=rmfield(s,'crossings');
end

function shape=region(s,g,p)
% Piecewise-linear superlevel geometry of the original active-owner field.
% The threshold sees BOTH active copies; the core owns its closed rectangle.
% Enclose every above-level island plus one exterior polar ring, so distant
% subthreshold tails do not spoil the near-field Delaunay coordinate scale.
assert(isstruct(s)&&isfield(s,'rhs')&&isfield(s.rhs,'q')&& ...
    isequal(size(s.rhs.q.core),size(g.core.X))&& ...
    isequal(size(s.rhs.q.polar),size(g.polar.X)), ...
    'ipm_cr:RemeshRegion','Use native q arrays matching both patches.');
q=s.rhs.q;ca=g.core.active;pa=g.polar.active;
cx=g.core.X(ca);cy=g.core.Y(ca);px=g.polar.X(pa);py=g.polar.Y(pa);
qc=q.core(ca);qp=q.polar(pa);[pc,ic]=max(qc);[pp,ip]=max(qp);
shape=struct('valid',false,'reason','nonpositive_peak','level',p.level, ...
    'definition','positive_peak_connected_component_projection_v1', ...
    'method','active_owner_local_linear_triangles_subcell_crossings', ...
    'peak',max(pc,pp),'peakXY',[NaN,NaN],'ownedPeak',NaN,'ownedPeakXY',[NaN,NaN], ...
    'peakDifference',NaN,'peakRelativeDifference',NaN,'threshold',NaN, ...
    'xBounds',[NaN,NaN],'yBounds',[NaN,NaN],'widths',[NaN,NaN],'extents',[NaN,NaN], ...
    'monitorRadius',NaN,'monitorNodes',0,'triangleCount',0,'componentNodes',0, ...
    'otherPositiveNodes',0,'crossings',zeros(0,2),'wallConnected',false, ...
    'globalActiveCount',numel(qc)+numel(qp),'ghostDonorsUsed',false, ...
    'overlapOwner','core_inside_closed_active_rectangle');
if pc>=pp,shape.peakXY=[cx(ic),cy(ic)];else,shape.peakXY=[px(ip),py(ip)];end
if ~(shape.peak>0),return;end
threshold=p.level*shape.peak;shape.threshold=threshold;
radius=max(hypot(cx,cy));above=qp>=threshold;
if any(above),radius=max(radius,max(hypot(px(above),py(above))));end
rr=g.polar.r(g.polar.activeRows);
outer=find(rr>radius+64*eps(max(1,radius)),1);
if isempty(outer),shape.reason='level_set_not_enclosed_by_active_domain';return;end
shape.monitorRadius=rr(outer);
inside=px>=min(cx)&px<=max(cx)&py>=min(cy)&py<=max(cy);
take=~inside&hypot(px,py)<=shape.monitorRadius+64*eps(max(1,shape.monitorRadius));
points=[cx,cy;px(take),py(take)];values=[qc;qp(take)];
[points,keep]=unique(points,'rows','stable');values=values(keep);
shape.monitorNodes=size(points,1);[shape.ownedPeak,seed]=max(values);
shape.ownedPeakXY=points(seed,:);shape.peakDifference=shape.peak-shape.ownedPeak;
shape.peakRelativeDifference=shape.peakDifference/shape.peak;
if shape.ownedPeak<threshold,shape.reason='peak_absent_from_owner_field';return;end
triangles=delaunay(points(:,1),points(:,2));shape.triangleCount=size(triangles,1);
edges=unique(sort([triangles(:,[1,2]);triangles(:,[2,3]);triangles(:,[3,1])],2),'rows');
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
crossed=edges(xor(component(edges(:,1)),component(edges(:,2))),:);
a=crossed(:,1);b=crossed(:,2);
fraction=(threshold-values(a))./(values(b)-values(a));
shape.crossings=unique(points(a,:)+fraction.*(points(b,:)-points(a,:)),'rows','stable');
support=[shape.crossings;points(component,:)];
shape.xBounds=[min(support(:,1)),max(support(:,1))];
shape.yBounds=[min(support(:,2)),max(support(:,2))];
shape.widths=[diff(shape.xBounds),diff(shape.yBounds)];
shape.extents=[max(abs(shape.xBounds)),shape.yBounds(2)];
shape.valid=true;shape.reason='ok';
end

function m=monitor(s,g)
% Both source projections use the original active core, never ghost values.
c=g.core;rows=c.activeRows;cols=c.activeColumns;
x=c.x(cols);y=c.y(rows);q=s.rhs.q.core(rows,cols);
ix=unique([1:2:numel(x),numel(x)]);iy=unique([1:2:numel(y),numel(y)]);
width=min([7,numel(ix),numel(iy)]);
mx=ipm_cr.mesh.localMap(x(ix),x,width);my=ipm_cr.mesh.localMap(y(iy),y,width);
ex=abs(q-q(:,ix)*mx');ey=abs(q-my*q(iy,:));
scale=max(max(abs(q),[],'all'),realmin);
direction=[max(ex,[],'all'),max(ey,[],'all')]/scale;
% Smooth magnitude plus unresolved content defines the density monitor.
px=max(abs(q),[],1)+max(ex,[],1);py=max(abs(q),[],2)+max(ey,[],2);
px=smooth_profile(px(:));py=smooth_profile(py(:));
m=struct('indicator',max(direction),'directional',direction, ...
    'sourceInfinity',scale,'activeX',x,'activeY',y,'profileX',px,'profileY',py, ...
    'definition','active-core every-other-node local polynomial source discrepancy', ...
    'independentContinuumError',false,'ghostValuesRead',false);
end

function v=smooth_profile(v)
v=[v(1);v(1);v(:);v(end);v(end)];v=conv(v,[1;4;6;4;1]/16,'valid');
v=v/max(max(v),realmin);
end

function a=propose(g,m,p,reference)
c=g.core;a=axes_of(g);
px=interp1(m.activeX,m.profileX,reference.x,'linear','extrap');
py=interp1(m.activeY,m.profileY,reference.y,'linear','extrap');
px=max(0,min(1,px));py=max(0,min(1,py));
origin=c.originIndex;
left=redistribute(reference.x(1:origin),c.x(1:origin),px(1:origin),p);
right=redistribute(reference.x(origin:end),c.x(origin:end),px(origin:end),p);
x=[left(1:end-1),right];y=redistribute(reference.y',c.y',py',p)';
if isequal(x,c.x)&&isequal(y,c.y),return;end
a.x=x;a.y=y;
% D7 supplies numerical coordinate Jacobians, not analytic sinh metrics.
dxi=2/(numel(x)-1);deta=1/(numel(y)-1);
a.dxMetric=(ipm_cr.mesh.fd(linspace(-1,1,numel(x)),1,7)*x')'*dxi;
a.dyMetric=ipm_cr.mesh.fd(linspace(0,1,numel(y)),1,7)*y*deta;
end

function out=redistribute(reference,current,profile,p)
x=reference(:)';current=current(:)';profile=reshape(profile,1,[]);n=numel(x);
% The fixed reference prevents repeated multiplication of mesh density.
t=linspace(0,1,n);weight=1+p.monitorStrength*profile;
mass=cumtrapz(t,weight);mass=mass/mass(end);
target=interp1(mass,x,t,'pchip');motion=target-current;
spacing=[diff(current(1:2)),(current(3:end)-current(1:end-2))/2,diff(current(end-1:end))];
ratio=max(abs(motion)./spacing);fraction=p.blend;
if ratio>0,fraction=min(fraction,p.maxNodeMotion/ratio);end
% One convex step preserves ordering of both current and target axes.
out=current+fraction*motion;out([1,end])=x([1,end]);
end

function a=axes_of(g)
c=g.core;p=g.polar;
if isfield(c,'dxMetric'),dx=c.dxMetric;dy=c.dyMetric;
else,dx=c.dx*ones(1,c.nx);dy=c.dy*ones(c.ny,1);end
a=struct('x',c.x,'y',c.y,'r',p.r,'theta',p.theta, ...
    'dxMetric',dx,'dyMetric',dy,'drMetric',p.drMetric);
end

function [out,audit]=transfer(s,old,new)
out=s;key='W';if isfield(s,'eta'),key='eta';end
[out.(key),audit]=transfer_field(s.(key),old,new);
if strcmp(key,'eta'),out.W=struct('core',zeros(size(new.core.X)),'polar',zeros(size(new.polar.X)));end
out.ctxId='';out.ctxKey=[];
out.rhs=struct('c_l',s.rhs.c_l,'c_r',s.rhs.c_r);
audit.integratedField=key;audit.timeAdvanced=0;audit.originalReferencesPreserved=true;
end

function [out,audit]=transfer_field(input,old,new)
out=struct();counts=struct();maxNorm=0;
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
    polarCovered=covered(r,[old.polar.innerActiveRadius,old.polar.r(end)])& ...
        covered(theta,[old.polar.theta(1),old.polar.theta(end)]);
    fixedFrame=same_frame(old.(key),target,key);
    if strcmp(key,'core'),coveredBySame=coreCovered;else,coveredBySame=polarCovered;end
    same=coveredBySame;
    if fixedFrame,same(:)=true;end
    if strcmp(key,'core'),coreOwner=same;else,coreOwner=~same;end
    % A fixed-frame remesh changes coordinates of independent patch DOFs,
    % not their chart. Continue the nearest local polynomial into the newly
    % exposed strip; cross charts only when the patch frame itself changes.
    values=zeros(size(x));
    for owner={'core','polar'}
        name=owner{1};take=coreOwner;if strcmp(name,'polar'),take=~take;end
        if ~any(take),continue;end
        p=old.(name);rows=p.activeRows;columns=p.activeColumns;
        if strcmp(name,'core')
            M=ipm_cr.mesh.tensorMap(p.y(rows),p.x(columns),y(take),x(take),7,fixedFrame);
        else
            M=ipm_cr.mesh.tensorMap(p.r(rows),p.theta(columns),r(take),theta(take),7,fixedFrame);
        end
        donor=input.(name)(rows,columns);
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
    'boundaryContinuation','same-patch nearest local degree-six polynomial on fixed-frame redistribution', ...
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
