function [rhoNew,opsNew,info] = directLevelSetAdapt(rho,ops,flow,config,amplitudeScale)
%IPM.REMESH.DIRECTLEVELSETADAPT One inverse-geometry quadrant remesh transaction.
% The measured 90% wall-core/front intervals and vertical core width are
% the only geometric inputs. There is one axis pair and one native transfer.
assert(ipm.mesh.isQuadrant(ops) && strcmp(ops.spatialDiscretization,'high_order'), ...
    'ipm:DirectLevelSetMode', ...
    'Direct level-set remeshing currently requires a high-order quadrant.');
if nargin<5,amplitudeScale=1;end
rhoNew=rho;opsNew=ops;
info=struct('applied',false,'count',ops.remeshCount, ...
    'kind','quadrant_direct_level_set_v1','status','not_triggered', ...
    'oldPeakSpacing',NaN,'newPeakSpacing',NaN, ...
    'oldWallSpacing',ops.y(2)-ops.y(1), ...
    'newWallSpacing',ops.y(2)-ops.y(1), ...
    'targetCells',[NaN,NaN],'actualCells',[NaN,NaN], ...
    'predictedCells',[NaN,NaN,NaN],'feature',struct(), ...
    'xQuality',struct(),'yQuality',struct(), ...
    'massRelativeDefect',NaN,'relativePeakJump',NaN, ...
    'relativeRangeViolation',NaN,'geometryProposals',0,'nativeTransfers',0, ...
    'geometryFailureIdentifier','','geometryFailureMessage','');
view=struct('rho',rho,'x',ops.x,'y',ops.y,'Dx',ops.Dx, ...
    'source',flow.source,'trusted',true);
feature=ipm.diagnostics.meshFeatureIntervals(view);
info.feature=feature;
[targets,triggerCells]=ipm.remesh.directLevelSetTargets(ops);
targetX=targets(1);targetY=targets(2);
targetFront=min(20,max(8,round(.6*targetX)));
info.targetCells=[targetX,targetY];
info.actualCells=feature.actualCoreCells;
info.oldPeakSpacing=local_spacing(ops.x,feature.coreCenter);
info.newPeakSpacing=info.oldPeakSpacing;
if all(feature.actualCoreCells>=triggerCells)
    return
end
if ~isfinite(feature.yCoreWidth) || feature.yCoreWidth<=0 || ...
        feature.coreCenter<=0 || any(diff(feature.coreInterval)<=0) || ...
        any(diff(feature.frontInterval)<=0)
    info.status='invalid_level_set_geometry';return
end

core=feature.coreInterval;front=feature.frontInterval;
spacing=min(diff(core)/targetX,diff(front)/targetFront)/1.15;
left=min(core(1),front(1));right=max(core(2),front(2));
fineCells=max(8,ceil(diff(core)/spacing));
fraction=(feature.coreCenter-left)/(right-left);
leftCells=max(1,min(fineCells-1,round(fineCells*fraction)));
fraction=leftCells/fineCells;
ratioCap=min(config.remesh.remeshMaximumCellRatio,1.20);
rounding=max(8,ceil(log(ratioCap)/.012));
rounding=min(rounding,max(8,floor((ops.nx-1-fineCells)/2)-2));
if ~(isfinite(spacing) && spacing>0 && left>0 && ...
        right<ops.x(end) && fineCells+2*(rounding+2)<=ops.nx-1)
    info.status='axis_capacity_or_boundary';return
end
info.geometryProposals=1;
try
    [x,xConstruction]=ipm.remesh.corePatchAxis( ...
        ops.rescaling.transportAnchorX,feature.coreCenter,ops.x(end), ...
        ops.nx,spacing,fineCells,fraction,rounding,.25,true);
    [y,yConstruction]=ipm.remesh.roundedGeometricAxis( ...
        ops.y,feature.yCoreWidth,1.15*targetY, ...
        min(8,max(2,floor((ops.ny-1)/4))));
catch exception
    if direct_geometry_rejection(exception)
        identifier=exception.identifier;
        if isempty(identifier)
            identifier='unidentified_geometry_failure';
        end
        info.status=['geometry_rejected:' identifier];
        info.geometryFailureIdentifier=exception.identifier;
        info.geometryFailureMessage=exception.message;
        return
    end
    rethrow(exception)
end
info.xConstruction=xConstruction;
info.yConstruction=yConstruction;
qx=ipm.mesh.quality(x,7,ipm.mesh.quadrantQuadrature(x));
qy=ipm.mesh.quality(y,7,ipm.mesh.quadrature(y));
info.xQuality=qx;info.yQuality=qy;
info.predictedCells=[interval_count(x,core), ...
    interval_count(y,[0,feature.yCoreWidth]),interval_count(x,front)];
if qx.maximumAdjacentCellRatio>ratioCap || ...
        qy.maximumAdjacentCellRatio>ratioCap || ...
        ~qx.quadratureWeightsStrictlyPositive || ...
        ~qy.quadratureWeightsStrictlyPositive || ...
        min(qx.minimumStencilRcond,qy.minimumStencilRcond)<1e-9 || ...
        any(info.predictedCells(1:2)<.75*[targetX,targetY]) || ...
        info.predictedCells(3)<.75*targetFront || ...
        ~any(x==ops.rescaling.transportAnchorX)
    info.status='geometric_quality_rejected';return
end
if isequal(x,ops.x) && isequal(y,ops.y)
    info.status='same_axes';return
end
proposal=struct('x',x,'y',y);
info.nativeTransfers=1;
[candidateRho,candidateOps,metrics]= ...
    ipm.remesh.transfer(rho,ops,config,proposal);
oldMass=sum(rho.*ops.integrationWeights,'all');
newMass=sum(candidateRho.*candidateOps.integrationWeights,'all');
info.massRelativeDefect=abs(newMass-oldMass)/max(abs(oldMass),eps);
info.relativePeakJump=abs(metrics.peak-flow.trackedWallPeak)/ ...
    max(abs(flow.trackedWallPeak),eps);
info.relativeRangeViolation=metrics.relativeRangeViolation;
if info.massRelativeDefect>5e-12 || ...
        info.relativePeakJump>config.remesh.remeshPeakChangeTolerance || ...
        metrics.rangeViolation/amplitudeScale> ...
            config.diagnostics.rangeStopTolerance || ...
        metrics.xCore<.75*targetX || metrics.yCore<.75*targetY
    info.status='native_transfer_rejected';return
end
candidateOps.remeshCount=ops.remeshCount+1;
rhoNew=candidateRho;opsNew=candidateOps;
info.applied=true;info.count=candidateOps.remeshCount;
info.status='accepted';
info.newPeakSpacing=local_spacing(x,feature.coreCenter);
info.newWallSpacing=y(2)-y(1);
end

function count=interval_count(axis,bounds)
v=axis(:)';
count=sum(max(0,min(v(2:end),bounds(2))- ...
    max(v(1:end-1),bounds(1)))./diff(v));
end

function spacing=local_spacing(axis,point)
[~,index]=min(abs(axis-point));
left=max(index-1,1);right=min(index+1,numel(axis));
spacing=(axis(right)-axis(left))/(right-left);
end

function answer=direct_geometry_rejection(exception)
identifiers={'ipm:RoundedAxisInfeasible','ipm:RoundedGeometricInfeasible', ...
    'ipm:CorePatchAxisInfeasible'};
legacyMessages={'No feasible direct rounded-log cell split.', ...
    'Insufficient rounded-branch cell budget.', ...
    'No feasible rounded-log cell split.'};
answer=any(strcmp(exception.identifier,identifiers)) || ...
    any(strcmp(exception.message,legacyMessages));
end
