function [candidates,report,cache]=plannedAxisPairs(view,referenceAxes,anchor,policy,cache)
%IPM.REMESH.PLANNEDAXISPAIRS Bounded no-LU, no-transfer axis proposals.
% policy is already resolved by config; this function never resolves config,
% modifies a source/reference, builds a flow, or certifies a native transfer.
feature=ipm.diagnostics.meshFeatureIntervals(view);
if nargin<5,cache=struct();end
assert(isstruct(cache)&&isscalar(cache),'ipm:AutonomousMeshAxisCache');
validate_policy(policy);x=view.x(:)';y=view.y(:);validate_reference(referenceAxes,x,y,anchor,policy);
context=struct('sourceX',x,'sourceY',y,'feature',feature, ...
    'anchor',anchor,'policy',policy);
if isfield(cache,'context')
    assert(isequaln(cache.context,context),'ipm:AutonomousMeshAxisCache', ...
        'Axis trials may only be reused within the same source and policy.');
else
    assert(isempty(fieldnames(cache)),'ipm:AutonomousMeshAxisCache');
    cache.context=context;
end
newNx=numel(referenceAxes.x);newNy=numel(referenceAxes.y);
sameN=newNx==numel(x)&&newNy==numel(y);
search=policy.search;limits=policy.qualityLimits;
meshLimits=rmfield(limits,{'minWeightToControlWidth','maxWeightToControlWidth'});
core=feature.coreInterval;front=feature.frontInterval;
spacing=min(diff(core)/policy.targetCoreCells(1),diff(front)/policy.minimumFrontCells)/search.xPadding;
report=struct('kind','ipm_bounded_planned_axis_pairs_v1','policy',policy,'feature',feature, ...
    'nodeCount',[newNx,newNy],'xLimits',x([1,end]),'yLimits',y([1,end])', ...
    'anchor',anchor,'inputAnchorExact',any(x==anchor)&&any(x==-anchor), ...
    'xTrials',struct([]),'yTrials',struct([]),'selectedPairs',struct([]),'status','planning', ...
    'noLU',true,'noTransfer',true,'noFlow',true,'sameBoxSameNodeCount',sameN,'globalFeasibilityClaim',false);
if isfield(policy,'version')&&any(policy.version == [2,3,4,5])
    report.kind='ipm_bounded_planned_axis_pairs_v2';
    report.actualSourceNodeCount=[numel(x),numel(y)];if policy.version == 3,report.kind='ipm_bounded_planned_axis_pairs_v3';end;if policy.version == 4,report.kind='ipm_bounded_planned_axis_pairs_v4';end;if policy.version == 5,report.kind='ipm_bounded_planned_axis_pairs_v5';end
end
candidates=struct('x',{},'y',{},'unchanged',{},'xIndex',{},'yIndex',{},'predictedCells',{},'quality',{});
xRows=struct([]);xAxes={};yRows=struct([]);yAxes={};
reuseY=false;
if isfield(cache,'yEntries')
    for entry=cache.yEntries
        if isequal(entry.reference,referenceAxes.y)
            yRows=entry.rows;yAxes=entry.axes;reuseY=true;break
        end
    end
end
if ~reuseY
for yTrialIndex=1:numel(search.ySigma)
    sigma=search.ySigma(yTrialIndex);
    if isfield(policy,'densityVersion') && any(policy.densityVersion==[1,2]) && ...
            yTrialIndex==numel(search.ySigma)
        try
            [v,info]=ipm.remesh.roundedGeometricAxis(referenceAxes.y(:), ...
                feature.yCoreWidth,policy.targetCoreCells(2)*search.yPadding,8);
        catch e
            if ~strcmp(e.identifier,'ipm:RoundedGeometricInfeasible'),rethrow(e);end
            v=[];info=struct('kind','one_sided_rounded_geometric_equalizer_v1', ...
                'status','ratio_capacity_exhausted');
        end
    else
        [v,info]=ipm.remesh.equalizedAxis(referenceAxes.y(:)/anchor, ...
            [0,feature.yCoreWidth/anchor],policy.targetCoreCells(2)*search.yPadding, ...
            struct('geometry','positive','focusCenters',0,'sigma',sigma,'power',2,'meshLimits',meshLimits));
        v=v(:)*anchor;v(1)=y(1);v(end)=y(end);
    end
    if isempty(v)
        q=struct();reasons={'geometric_ratio_capacity'};margin=-Inf;count=NaN;
    else
        [q,reasons,margin]=axis_quality(v,limits);
        count=interval_count(v,[0,min(v(end),feature.yCoreWidth)]);
        reasons=with_reasons(reasons,count<policy.targetCoreCells(2)-1e-6,{'y_core_target'});
    end
    row=struct('sigma',sigma,'quality',q,'coreCells',count,'admissible',isempty(reasons), ...
        'reasons',{reasons},'qualityMargin',margin,'equalizerInfo',info);
    yRows=append_row(yRows,row);yAxes{end+1}=v; %#ok<AGROW>
end
    entry=struct('reference',referenceAxes.y,'rows',yRows,'axes',{yAxes});
    if ~isfield(cache,'yEntries')||isempty(cache.yEntries)
        cache.yEntries=entry;
    else
        cache.yEntries(end+1)=entry;
    end
end
reuseX=false;
if isfield(cache,'xEntries')
    for entry=cache.xEntries
        if isequal(entry.reference,referenceAxes.x)
            xRows=entry.rows;xAxes=entry.axes;reuseX=true;break
        end
    end
end
if ~reuseX
for fine=search.fineCells
 for rounding=search.roundingCells
  for fraction=search.coreFineCellFractions
    row=struct('positiveFineCells',fine,'roundingCells',rounding,'coreFineCellFraction',fraction, ...
        'admissible',false,'reasons',{{}},'quality',struct(),'coreCells',NaN,'frontCells',NaN, ...
        'qualityMargin',-Inf,'construction',struct(),'failure',struct());v=[];
    if fine+2*(rounding+2)>(newNx-1)/2
        row.reasons={'rounded_branch_node_budget'};
    elseif feature.coreCenter<=0||core(2)<=core(1)||front(2)<=front(1)
        row.reasons={'unsupported_or_degenerate_core_geometry'};
    else
        try
            [v,construction]=ipm.remesh.corePatchAxis(anchor,feature.coreCenter,x(end),newNx, ...
                spacing,fine,fraction,rounding,search.maximumWarpFraction);
            [q,reasons,margin]=axis_quality(v,limits);
            nx=interval_count(v(v>=0),core);nf=interval_count(v(v>=0),front);
            reasons=with_reasons(reasons,[nx<policy.targetCoreCells(1)-1e-6,nf<policy.minimumFrontCells-1e-6, ...
                ~any(v==anchor)||~any(v==-anchor)],{'x_core_target','x_front_target','exact_anchor'});
            row.quality=q;row.coreCells=nx;row.frontCells=nf;row.qualityMargin=margin;
            row.construction=construction;row.reasons=reasons;row.admissible=isempty(reasons);
        catch e
            assert_geometry_rejection(e);
            row.reasons={'geometry_construction_rejected'};
            row.failure=struct('identifier',e.identifier,'message',e.message,'stack',e.stack);
        end
    end
    xRows=append_row(xRows,row);xAxes{end+1}=v; %#ok<AGROW>
  end
 end
end
    entry=struct('reference',referenceAxes.x,'rows',xRows,'axes',{xAxes});
    if ~isfield(cache,'xEntries')||isempty(cache.xEntries)
        cache.xEntries=entry;
    else
        cache.xEntries(end+1)=entry;
    end
end
report.xTrials=xRows;report.yTrials=yRows;
pairIndex=[];ranks=[];ratios=[];
[xQuality,xReasons]=axis_quality(x,limits);[yQuality,yReasons]=axis_quality(y,limits);
oldCounts=[interval_count(x(x>=0),core),interval_count(y,[0,min(y(end),feature.yCoreWidth)]),interval_count(x(x>=0),front)];
if sameN&&any(x==anchor)&&any(x==-anchor)&&isempty(xReasons)&&isempty(yReasons)&& ...
        all(oldCounts>=[policy.targetCoreCells,policy.minimumFrontCells]-1e-6)
    pairIndex=[0,0];ranks=Inf;ratios=-Inf;
end
for i=find([xRows.admissible])
 for j=find([yRows.admissible])
    pairIndex(end+1,:)=[i,j];ranks(end+1,1)=min(xRows(i).qualityMargin,yRows(j).qualityMargin); %#ok<AGROW>
    ratios(end+1,1)=max(xRows(i).quality.maximumAdjacentCellRatio, ...
        yRows(j).quality.maximumAdjacentCellRatio); %#ok<AGROW>
 end
end
report.admittedAxisPairCount=size(pairIndex,1);
if isempty(pairIndex)
    if ~any([yRows.admissible]),report.status='registered_y_family_infeasible';
    else,report.status='registered_x_family_infeasible';end
    return
end
if isfield(policy,'version')&&policy.version==5
    [~,order]=sortrows([ratios(:),-ranks(:),(1:numel(ranks))'],[1,2,3]);
    report.ranking='minimax adjacent ratio, then quality margin; qualified keep first';
else
    [~,order]=sort(ranks,'descend');
    report.ranking='weakest dimensionless quality margin, stable enumeration ties; qualified keep first';
end
for k=1:min(search.maximumPairCandidates,numel(order))
    ij=pairIndex(order(k),:);unchanged=all(ij==0);
    if unchanged
        vx=x;vy=y;qx=xQuality;qy=yQuality;counts=oldCounts;
    else
        vx=xAxes{ij(1)};vy=yAxes{ij(2)};qx=xRows(ij(1)).quality;qy=yRows(ij(2)).quality;
        counts=[xRows(ij(1)).coreCells,yRows(ij(2)).coreCells,xRows(ij(1)).frontCells];
    end
    row=struct('x',vx,'y',vy,'unchanged',unchanged,'xIndex',ij(1),'yIndex',ij(2), ...
        'predictedCells',counts,'quality',struct('x',qx,'y',qy));
    candidates=append_row(candidates,row);
    report.selectedPairs=append_row(report.selectedPairs,struct('xIndex',ij(1),'yIndex',ij(2), ...
        'unchanged',unchanged,'predictedCells',counts,'rankingMargin',ranks(order(k))));
end
report.status='axes_proposed_require_actual_transfer';
end

function validate_policy(p)
assert(isstruct(p)&&isscalar(p)&&all(isfield(p,{'targetCoreCells','minimumFrontCells','qualityLimits','search'})), ...
    'ipm:AutonomousMeshPolicy','Pass the resolved autonomous mesh policy.');
validateattributes(p.targetCoreCells,{'numeric'},{'row','numel',2,'integer','>=',21});
validateattributes(p.minimumFrontCells,{'numeric'},{'scalar','integer','>=',20});
q=p.qualityLimits;names={'maxAdjacentCellRatio','maxLogSpacingCurvature','minStencilRcond', ...
    'minQuadratureWeightRatio','minWeightToControlWidth','maxWeightToControlWidth'};
assert(isstruct(q)&&isscalar(q)&&isempty(setxor(fieldnames(q),names)));
assert(q.maxAdjacentCellRatio<=1.08&&q.maxLogSpacingCurvature<=.01&&q.minStencilRcond>=1e-9&& ...
    q.minQuadratureWeightRatio>=1e-8&&q.minWeightToControlWidth>=.35&&q.maxWeightToControlWidth<=1.65);
s=p.search;required={'fineCells','roundingCells','coreFineCellFractions','ySigma','maximumAxisCandidates', ...
    'maximumPairCandidates','xPadding','yPadding','maximumWarpFraction'};
assert(isstruct(s)&&isscalar(s)&&all(isfield(s,required)));
validateattributes(s.fineCells,{'numeric'},{'row','integer','>=',8});
validateattributes(s.roundingCells,{'numeric'},{'row','integer','>=',2});
validateattributes(s.coreFineCellFractions,{'numeric'},{'row','finite','>',0,'<',1});
validateattributes(s.ySigma,{'numeric'},{'row','finite','positive'});
validateattributes(s.xPadding,{'numeric'},{'scalar','finite','>=',1});
validateattributes(s.yPadding,{'numeric'},{'scalar','finite','>=',1});
validateattributes(s.maximumWarpFraction,{'numeric'},{'scalar','finite','>',0,'<',.9});
validateattributes(s.maximumAxisCandidates,{'numeric'},{'scalar','integer','positive','<=',500});
validateattributes(s.maximumPairCandidates,{'numeric'},{'scalar','integer','positive','<=',3});
assert(numel(s.fineCells)*numel(s.roundingCells)*numel(s.coreFineCellFractions)<=s.maximumAxisCandidates&&numel(s.ySigma)<=16, ...
    'ipm:AutonomousMeshBudget','The complete search must fit its predeclared finite budget.');
end

function validate_reference(r,x,y,a,policy)
assert(isstruct(r)&&isscalar(r)&&all(isfield(r,{'x','y'})));
validateattributes(a,{'numeric'},{'scalar','finite','positive','<',x(end)});
rx=r.x(:)';ry=r.y(:);
if isfield(policy,'version')&&any(policy.version == [2,3,4,5])
    validateattributes(rx,{'numeric'},{'vector','real','finite','increasing'});
    validateattributes(ry,{'numeric'},{'vector','real','finite','increasing'});
    assert(numel(rx)>=7&&mod(numel(rx),2)==1&&numel(ry)>=7, ...
        'ipm:AutonomousMeshReference','Target references require an odd x axis and at least seven nodes per axis.');
    assert(isfield(policy,'nodeFamily')&&isfield(policy.nodeFamily,'maximumTotalNodes')&& ...
        numel(rx)*numel(ry)<=policy.nodeFamily.maximumTotalNodes, ...
        'ipm:AutonomousMeshResourceCap','Target node count exceeds its explicit resource cap.');
else
    validateattributes(rx,{'numeric'},{'real','finite','increasing','numel',numel(x)});
    validateattributes(ry,{'numeric'},{'real','finite','increasing','numel',numel(y)});
end
assert(isequal(rx,-fliplr(rx))&&ry(1)==0&&isequal(x([1,end]),rx([1,end]))&&isequal(y([1,end]),ry([1,end])), ...
    'ipm:AutonomousMeshReference','The immutable reference must preserve actual units and domain endpoints.');
end

function [q,reasons,margin]=axis_quality(axis,limits)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
v=[q.maximumAdjacentCellRatio,q.maximumLogSpacingCurvature,q.minimumStencilRcond, ...
    q.minimumQuadratureWeightRatio,q.minimumQuadratureWeightToControlWidthRatio,q.maximumQuadratureWeightToControlWidthRatio];
bad=[v(1)>limits.maxAdjacentCellRatio,v(2)>limits.maxLogSpacingCurvature,v(3)<limits.minStencilRcond, ...
    v(4)<limits.minQuadratureWeightRatio,v(5)<limits.minWeightToControlWidth,v(6)>limits.maxWeightToControlWidth];
names={'adjacent_ratio','log_spacing_curvature','stencil_rcond','global_quadrature','local_quadrature_min','local_quadrature_max'};
reasons=names(bad);if ~q.quadratureWeightsStrictlyPositive,reasons{end+1}='nonpositive_quadrature';end
margin=min([limits.maxAdjacentCellRatio/v(1),limits.maxLogSpacingCurvature/max(v(2),realmin), ...
    v(3)/limits.minStencilRcond,v(4)/limits.minQuadratureWeightRatio, ...
    v(5)/limits.minWeightToControlWidth,limits.maxWeightToControlWidth/v(6)]);
end

function count=interval_count(axis,interval)
v=axis(:)';count=sum(max(0,min(v(2:end),interval(2))-max(v(1:end-1),interval(1)))./diff(v));
end
function reasons=with_reasons(reasons,failed,names)
reasons=[reasons,names(failed)];
end
function rows=append_row(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function assert_geometry_rejection(e)
known={'No feasible rounded-log cell split.','Insufficient rounded-branch cell budget.'};
fromFactory=~isempty(e.stack)&&any(strcmp(e.stack(1).name, ...
    {'ipm.remesh.corePatchAxis','ipm.remesh.roundedAxis','ipm.remesh.roundedAxis>growth', ...
    'corePatchAxis','roundedAxis','roundedAxis>growth'}));
if ~(fromFactory&&strcmp(e.identifier,'MATLAB:assertion:failed'))&&~any(strcmp(e.message,known)),rethrow(e);end
end
