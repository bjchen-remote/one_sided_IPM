function [state,context,info]=remesh(state,context,settings,force)
%REMESH Shape-triggered mapped core and annulus reconstruction.
% Rebuild the mapped patches, transfer active fields and evaluate one complete RHS.
% Calibrate the wall-connected main 25% and 80% q regions on the first check.
% A 25% change from either last-remesh anchor triggers mesh adaptation.
if nargin<3
    settings=struct();
    if isa(context,'ipm_cr.context')&&isfield(context.opts,'remesh'),settings=context.opts.remesh;end
end
if nargin<4,force=false;end
p=ipm_cr.mesh.region.policy(settings);
tau=state.scale.canonicalTime;g=context.mesh;
memory=controller(state);
nLevels=numel(p.level);
info=struct('status','not_triggered','reason','disabled','remeshed',false, ...
    'timeAdvanced',0,'forced',force,'tau',tau,'step',state.step, ...
    'checked',false,'triggered',false,'directionalChange',zeros(nLevels,2), ...
    'triggerMask',false(nLevels,2),'directionalRateEstimate',nan(nLevels,2), ...
    'checkCadence',p.checkEveryTau, ...
    'independentGreenAuditPerformed',false,'continuumAccuracyCertified',false);
if ~isempty(p.stopTau)&&tau>=p.stopTau
    info.reason='past_stop_tau';return
end
if ~p.enabled&&~force,return;end
if ~isempty(memory.referenceShape)&&~force&&tau<memory.nextCheckTau
    info.reason='check_cadence';return
end
shape=ipm_cr.mesh.region.measure(state,g,p);
require_regions(shape);
extents=vertcat(shape.triggerExtents);
info.checked=true;info.region=shape;memory.checks=memory.checks+1;
if isempty(memory.referenceShape)
    memory.referenceShape=small_shape(shape);
    memory.initialExtents=shape(1).extents;
    memory.initialMesh=ipm_cr.mesh.axes.reference(g);
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
    [axes,info.geometryProposal]=ipm_cr.mesh.axes.propose(shape(1).extents,memory);
else
    axes=p.candidateAxes;
end
staircase=isfield(g.core,'stair');
if ~staircase&&isequal(axes.x(:),g.core.x(:))&&isequal(axes.y(:),g.core.y(:))&& ...
        isequal(axes.r(:),g.polar.r(:))&&isequal(axes.theta(:),g.polar.theta(:))
    info.reason='geometry_unchanged';finish_controller();return
end
newG=ipm_cr.mesh.axes.rebuild(g,axes);
if staircase
    if isequaln(newG.geometryKey,g.geometryKey)
        info.reason='geometry_unchanged';finish_controller();return
    end
    info.geometryProposal.innerRadius=newG.polar.r(1);
    info.geometryProposal.activeRadialLayers=newG.core.stair.activeRadialLayers;
end
[moved,info.transfer]=ipm_cr.mesh.transfer(state,g,newG,false);
% Actual geometry changes; the original method prescription does not.
newC=ipm_cr.context(newG,context.wall,context.opts);
moved.ctxId=newC.id;moved.ctxKey=newC.key;
[moved.rhs,moved.W]=ipm_cr.rhs(moved,newC);
memory.remeshCount=memory.remeshCount+1;
memory.lastRemeshTau=tau;
afterShape=ipm_cr.mesh.region.measure(moved,newG,p);
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

function s=small_shape(s)
s=rmfield(s,'crossings');
end

function require_regions(shapes)
for shape=shapes
    assert(shape.valid,'ipm_cr:RemeshRegion', ...
        'The %.0f%% wall-connected level set is undefined: %s.',100*shape.level,shape.reason);
end
end
