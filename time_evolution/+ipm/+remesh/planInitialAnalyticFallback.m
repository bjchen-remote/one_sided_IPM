function [candidates,evidence] = planInitialAnalyticFallback(config,originalAxes,originalAttempt)
%IPM.REMESH.PLANINITIALANALYTICFALLBACK One fixed observation of original t0.
% Omit ORIGINALATTEMPT for a pure preflight before any solver operator is
% built; supply it for the original bounded failure and controller ledger.
% No incoming rho or checkpoint is accepted. Actual source state is untouched.
id='ipm:InitialObservationContext';
assert(isfield(config.remesh,'initialMeshObservationFallback'),id,'Explicit initial observation policy is required.');
choices=flat(config);rule=ipm.config.initialMeshObservationPolicy(config.remesh.initialMeshObservationFallback,choices);
assert(rule.enabled&&isequaln(rule,config.remesh.initialMeshObservationFallback),id,'The policy must be enabled and canonical.');
p=config.remesh.autonomousMesh;x=originalAxes.x;y=originalAxes.y;
assert(isrow(x)&&iscolumn(y)&&numel(x)==config.grid.nx&&numel(y)==config.grid.ny&& ...
    all(isfinite(x))&&all(isfinite(y))&&all(diff(x)>0)&&all(diff(y)>0)&& ...
    isequal(x,-fliplr(x))&&isequal(x([1,end]),config.grid.xlim)&&y(1)==0&&y(end)==config.grid.ymax, ...
    id,'Immutable initial solver axes must match the configured initial budget and box.');
assert(numel(x)*numel(y)<=p.nodeFamily.maximumTotalNodes,id,'Initial solver node budget is not resource-admitted.');
if isfield(config.grid,'customX')&&~isempty(config.grid.customX)
    assert(isequal(x,config.grid.customX)&&isequal(y,config.grid.customY),id,'Explicit original axes must be unchanged.');
end
preflight=nargin<3||isempty(originalAttempt);
if preflight
    originalAttempt=struct('plan',struct('coreCells',[NaN,NaN]), ...
        'attempts',struct([]));
else
    assert(isstruct(originalAttempt)&&isscalar(originalAttempt)&& ...
        all(isfield(originalAttempt,{'plan','attempts'})),id,'Original bounded failure evidence is required.');
end
plan=originalAttempt.plan;attempts=originalAttempt.attempts;
coreFailure=false;
if ~preflight
assert(plan.initial&&plan.requested&&plan.sourceStep==0&&plan.sourceCanonicalTime==0&& ...
    plan.sourcePhysicalTime==0&&plan.sourceNormalizedTime==0&&plan.sourceRemeshCount==0&& ...
    plan.sourceLevelId==1&&isequal(plan.sourceNodeCount,[numel(x),numel(y)])&& ...
    ~isfield(plan,'observationFallback'),id,'Only an uncommitted original zero-time source is eligible once.');
assert(numel(plan.candidates)<=p.search.maximumPairCandidates&&numel(attempts)<=p.search.maximumPairCandidates,id,'Original stage exceeded its bounded candidate budget.');
capacity=strcmp(plan.stopReason,'autonomous_mesh_axis_capacity')&&isempty(plan.candidates)&&isempty(attempts);
coreFailure=isempty(plan.stopReason)&&~isempty(plan.candidates)&&numel(attempts)==numel(plan.candidates);
if coreFailure
    for k=1:numel(attempts)
        a=attempts(k);c=plan.candidates(k);
        valid=isfield(a,'audit')&&isstruct(a.audit)&&isscalar(a.audit)&& ...
            all(isfield(a.audit,{'kind','passed','reasons','qualityPassed','exactZeroTime','coreCells','leftFrontCells'}));
        if ~valid,coreFailure=false;break;end
        z=a.audit;
        onlyResolution=isequal(z.reasons,{'initial_axis_or_analytic_resolution_gate'})&& ...
            all(isfinite([z.coreCells,z.leftFrontCells]))&& ...
            (any(z.coreCells<p.transactionMinimumCoreCells)||z.leftFrontCells<p.minimumFrontCells);
        coreFailure=coreFailure&&a.candidateIndex==k&&~a.passed&&~z.passed&& ...
            isempty(a.errorIdentifier)&&isempty(a.errorMessage)&& ...
            strcmp(z.kind,'analytic_zero_time_axis_selection_audit')&& ...
            z.qualityPassed&&z.exactZeroTime&&onlyResolution&&any(c.x==1)&&any(c.x==-1);
    end
end
assert(capacity||coreFailure,'ipm:InitialObservationIneligibleFailure', ...
    'Only original finite-family exhaustion or solely actual core/front rejection permits observation fallback.');
end
evidence=struct('version',1,'policy',rule,'used',true,'attemptCount',1, ...
    'selectedPhase','fixed_observation','status','planning', ...
    'originalAxes',originalAxes,'originalAttempt',originalAttempt, ...
    'sourceNodeCount',[numel(x),numel(y)],'sourceCoreCells',plan.coreCells, ...
    'originalFailureClass','finite_family_exhaustion','observation',struct(), ...
    'candidateAudits',struct([]),'admittedOriginalPairIndices',[], ...
    'selectedLocalCandidateIndex',0,'selectedOriginalPairIndex',0, ...
    'noSourceFieldInterpolation',true,'noNativeQualificationClaim',true);
if coreFailure,evidence.originalFailureClass='actual_analytic_core_or_front';end
if preflight
    evidence.used=false;evidence.attemptCount=0;
    evidence.selectedPhase='analytic_preflight';
    evidence.originalFailureClass='none_preflight';
end
positive=observation_axis(config.grid.xlim(2),rule);ox=[-fliplr(positive(2:end)),positive];
oy=observation_axis(config.grid.ymax,rule)';
count=[numel(ox),numel(oy)];assert(prod(count)<=rule.maximumObservationNodes, ...
    'ipm:InitialObservationBudget','Fixed observation grid exceeds the observation budget.');
qx=axis_quality(ox,p);qy=axis_quality(oy,p);
assert(qx.passed&&qy.passed,'ipm:InitialObservationQuality','Observation axes fail original quality gates.');
view=analytic_view(ox,oy,config);features=ipm.diagnostics.meshFeatureIntervals(view);
evidence.observation=struct('x',ox,'y',oy,'nodeCount',count,'xQuality',qx.quality,'yQuality',qy.quality, ...
    'features',features,'datum','original_k8_initialDensity','derivative','paired_seven_point_fdMatrix', ...
    'canonicalTime',0,'physicalTime',0,'normalizedTime',0,'solverBudgetChanged',false);
[pairs,axisReport]=ipm.remesh.plannedAxisPairs(view,originalAxes,1,p);
evidence.axisReport=axisReport;candidates=struct([]);
for k=1:numel(pairs)
    c=pairs(k);trial=analytic_view(c.x,c.y,config);f=ipm.diagnostics.meshFeatureIntervals(trial);
    a=axis_quality(c.x,p);b=axis_quality(c.y,p);family=struct();familyFailure=struct();
    geometry=a.passed&&b.passed&&all(f.actualCoreCells>=p.transactionMinimumCoreCells)&& ...
        f.leftFrontCells>=p.minimumFrontCells&&any(c.x==1)&&any(c.x==-1)&& ...
        isequal([numel(c.x),numel(c.y)],evidence.sourceNodeCount)&& ...
        isequal(c.x([1,end]),x([1,end]))&&isequal(c.y([1,end]),y([1,end]));
    try
        family=ipm.remesh.referenceAxisFamily(c.x,c.y,1,p);
        familyPassed=all([family.members.qualityPassed])&&family.members(1).resourceAdmitted;
    catch e
        familyPassed=false;familyFailure=struct('identifier',e.identifier,'message',e.message);
    end
    row=struct('originalPairIndex',k,'x',c.x,'y',c.y,'features',f, ...
        'xQuality',a.quality,'yQuality',b.quality,'geometryPassed',geometry, ...
        'familyPassed',familyPassed,'familyFailure',familyFailure, ...
        'resourceMask',[],'passed',geometry&&familyPassed,'nativeFlowEvaluated',false);
    if familyPassed,row.resourceMask=[family.members.resourceAdmitted];end
    if isempty(evidence.candidateAudits),evidence.candidateAudits=row;else,evidence.candidateAudits(end+1)=row;end
    if row.passed
        c.nodeFamilyIndex=1;
        if isempty(candidates),candidates=c;else,candidates(end+1)=c;end %#ok<AGROW>
        evidence.admittedOriginalPairIndices(end+1)=k;
    end
end
if isempty(candidates),evidence.status='fixed_observation_capacity_rejected';else,evidence.status='pure_candidates_ready';end
end

function axis=observation_axis(endpoint,p)
inner=p.innerHalfWidth;spacing=p.innerSpacing;ramp=p.tailLogSlopeRampCells;
span=endpoint-inner;maximumLogRatio=log(p.maximumTailRatio);
n=ceil(log1p((p.maximumTailRatio-1)*span/spacing)/maximumLogRatio)+ramp;
assert(isfinite(n)&&n>0&&n<=p.maximumObservationNodes,'ipm:InitialObservationBudget', ...
    'The fixed observation tail exceeds its one-dimensional budget.');
indices=(0:n-1)';exponent=indices-ramp/2;mask=indices<ramp;
exponent(mask)=indices(mask).^2/(2*ramp);
assert(spacing*n<span&&spacing*sum(exp(maximumLogRatio*exponent))>=span, ...
    'ipm:InitialObservationTail','The fixed tail equation has no registered bracket.');
rate=fzero(@(v)spacing*sum(exp(v*exponent))-span,[0,maximumLogRatio]);
tail=inner+cumsum(spacing*exp(rate*exponent));tail(end)=endpoint;
axis=[linspace(0,inner,161),tail'];
assert(all(diff(axis)>0)&&any(axis==1),'ipm:InitialObservationTail','The fixed tail lost monotonicity or anchor.');
end
function view=analytic_view(x,y,c)
[X,Y]=meshgrid(x,y);ops=struct('X',X,'Y',Y,'nx',numel(x),'ny',numel(y),'symmetryMode',c.physics.symmetryMode);
rho=ipm.field.initialDensity(ops,c.physics);Dx=ipm.mesh.fdMatrix(x,1,7);source=rho*Dx';
assert(isreal(rho)&&all(isfinite(rho),'all')&&all(isfinite(source),'all'), ...
    'ipm:InitialObservationDatum','Original analytic samples must be finite.');
view=struct('rho',rho,'x',x,'y',y,'Dx',Dx,'source',source,'trusted',true);
end
function out=axis_quality(x,p)
q=ipm.mesh.quality(x,7,ipm.mesh.quadrature(x));l=p.qualityLimits;
passed=q.quadratureWeightsStrictlyPositive&&q.maximumAdjacentCellRatio<=l.maxAdjacentCellRatio&& ...
    q.maximumLogSpacingCurvature<=l.maxLogSpacingCurvature&&q.minimumStencilRcond>=l.minStencilRcond&& ...
    q.minimumQuadratureWeightRatio>=l.minQuadratureWeightRatio&& ...
    q.minimumQuadratureWeightToControlWidthRatio>=l.minWeightToControlWidth&& ...
    q.maximumQuadratureWeightToControlWidthRatio<=l.maxWeightToControlWidth;
out=struct('quality',q,'passed',passed);
end
function values=flat(config)
values=struct();schema=ipm.config.schema();
for d=1:numel(schema.domainNames)
    group=config.(schema.domainNames{d});names=fieldnames(group);
    for k=1:numel(names),values.(names{k})=group.(names{k});end
end
end
