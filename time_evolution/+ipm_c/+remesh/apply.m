function [state,info]=apply(state,fixture)
%APPLY Complete-area radial remesh: propose, transfer, fresh audit, commit.
% The acceptance/retry transaction follows the practical +ipm architecture.
% The geometry adapter retains uniform theta and uses feature-driven radial
% proposals. Transfer uses the shared high-order interpolation kernel, with
% explicit polar-mass and fixed-wall-integral jump gates rather than the
% Cartesian row/column conservation correction. Rejected proposals retain
% the accepted field/grid/clock/reference values and leave the request pending.
% A second argument is exclusively for explicitly labelled manufactured tests.
if nargin<2,fixture=struct();end
assert(isstruct(fixture)&&isscalar(fixture),'ipm_c:RemeshFixture');
manufactured=~isempty(fieldnames(fixture));
if manufactured
    assert(isfield(fixture,'manufactured')&&isequal(fixture.manufactured,true), ...
        'ipm_c:RemeshFixture','Proposal overrides require an explicitly manufactured test.');
end
assert(all(isfield(state,{'rho','ops','flow','scale','config','step'}))&& ...
    all(isfield(state.ops,{'polarGrid','cart','polarCfg','rescaling'})), ...
    'ipm_c:RemeshState','Use an accepted total-density polar state.');
policy=state.config.remesh;ops=state.ops;g=ops.polarGrid;
assert(all(state.rho(1,:)==0),'ipm_c:CenteredDensityState');
info=struct('kind','polar_complete_area_remesh_transaction_v1','applied',false, ...
    'requested',false,'areaRequested',false,'pending',false,'deferred',false, ...
    'count',option(ops,'remeshCount',0),'status','disabled', ...
    'manufacturedFixture',manufactured,'levels',[.1,.2,.5], ...
    'areaChangeThreshold',option(policy,'areaChangeThreshold',.2), ...
    'maxAreaChange',NaN,'referenceAreas',[],'areas',[],'attempts',struct([]), ...
    'oldPolarMass',NaN,'newPolarMass',NaN, ...
    'transferPolicy','shared_degree_five_interpolation_with_polar_mass_and_wall_integral_jump_gates', ...
    'cartesianConservationCorrectionClaimed',false,'uniformThetaRetained',true);
if ~option(policy,'adaptiveRemesh',true),return;end
if info.count>=option(policy,'maxRemeshes',Inf),info.status='maximum_remeshes';return;end
levels=option(policy,'levels',[.1,.2,.5]);
assert(isequal(levels(:)',[.1,.2,.5]),'ipm_c:RemeshLevels','The current practical case uses [.1 .2 .5].');
[features,geometry]=ipm_c.diagnostics.sourceFeatures(state.flow.source,g,levels);
info.areas=features.areas;info.geometry=geometry;
if ~isfield(ops.rescaling,'referenceLevelSetAreas')
    if features.peak<=0,info.status='no_positive_source';return;end
    if info.count==0&&state.step>0
        initial=ipm.field.initialDensity(ops,ops.polarCfg.physics);
        derivative=ipm_c.field.derivatives(initial,g,struct('poleRegularization',false));
        initialFeatures=ipm_c.diagnostics.sourceFeatures(derivative.rx,g,levels);
    else
        initialFeatures=features;
    end
    assert(initialFeatures.usableForCdf,'ipm_c:AreaRemeshReference','The initial complete source must define positive areas.');
    ops.rescaling.referenceLevelSetAreas=initialFeatures.areas;
    ops.rescaling.referenceLevelSetRadialBounds=initialFeatures.radialProjectionBounds;
end
reference=ops.rescaling.referenceLevelSetAreas;
[areaRequested,change]=ipm_c.diagnostics.areaChange(features.areas,reference,info.areaChangeThreshold);
resolutionRequested=false;
if strcmp(option(ops.polarCfg.grid,'radialMap','sinh'),'uniform_core')&&features.usableForCdf
    requestedSpacing=max(option(policy,'minimumFirstSpacing',0),ops.polarCfg.grid.firstSpacing);
    targetSpacing=max(requestedSpacing,selected_core_radius(policy,ops.polarCfg.grid,features)/(.6*(g.nr-1)));
    resolutionRequested=g.r(2)>targetSpacing*(1+1e-8);
end
requested=areaRequested||resolutionRequested;
observer=option(ops.rescaling,'areaRemeshObserver',struct('checks',0,'failures',0,'retryAt',0,'features',[]));
assert(isstruct(observer)&&isscalar(observer)&&all(isfield(observer,{'checks','failures','retryAt','features'})), ...
    'ipm_c:AreaRemeshObserver','A saved retry observer must preserve its complete counters.');
observer.checks=observer.checks+1;ops.rescaling.areaRemeshObserver=observer;
state.ops=ops;info.referenceAreas=reference;info.maxAreaChange=change;
info.requested=requested;info.areaRequested=areaRequested;info.resolutionRequested=resolutionRequested;info.pending=requested;
info.status='not_triggered';
if ~requested
    observer.failures=0;observer.retryAt=0;observer.features=[];
    state.ops.rescaling.areaRemeshObserver=observer;return;
end
info.status='requested';
if observer.checks<observer.retryAt&&feature_change(features,observer.features)<.02-64*eps(.02)
    info.status='pending_backoff';info.deferred=true;info.retryAfterChecks=observer.retryAt-observer.checks;return;
end
accepted=state;
if ~features.usableForCdf
    [state,info]=reject(accepted,info,features,'ipm_c:AreaRemeshGeometry','Requested complete level sets have degenerate radial projections.');return;
end
try
    if isfield(fixture,'proposalAxis')
        target=fixture.proposalAxis(:);proposal=struct('kind','explicit_manufactured_axis');
    elseif strcmp(option(ops.polarCfg.grid,'radialMap','sinh'),'uniform_core')
        spacing=max(option(policy,'minimumFirstSpacing',0),ops.polarCfg.grid.firstSpacing);
        % Low-amplitude tails may expand while the peak stays narrow. They
        % belong on the stretched part, not in an ever coarser uniform core.
        % All levels still contribute to the complete-area trigger above.
        coreRadius=selected_core_radius(policy,ops.polarCfg.grid,features);
        settings=struct('nr',g.nr,'rmax',g.r(end),'firstSpacing',spacing, ...
            'coreRadius',coreRadius, ...
            'coreTailFraction',option(ops.polarCfg.grid,'coreTailFraction',.01));
        [target,uniformCoreInfo]=ipm_c.mesh.uniformCoreAxis(settings);
        proposal=struct('kind','complete_source_feature_uniform_core_axis_v1', ...
            'sourceFeatureBounds',features.featureBounds,'uniformCoreInfo',uniformCoreInfo, ...
            'coreRule','padded_highest_level_envelope_with_stretched_complete_tail', ...
            'coreLevel',levels(end),'corePadding',option(policy,'uniformCorePadding',2), ...
            'firstSpacing',target(2),'exactCdfQuantilesClaimed',false);
    else
        settings=struct('nr',g.nr,'rmax',g.r(end),'featureBounds',features.featureBounds, ...
            'coreScale',option(policy,'coreScale',ops.polarCfg.physics.gaussianWidthX/2), ...
            'featureMass',option(policy,'featureMass',6),'maximumAdjacentRatio',option(policy,'maximumCellRatio',6));
        [target,proposal]=ipm_c.remesh.radialCdf(settings);
    end
    assert(numel(target)==g.nr&&all(isfinite(target))&&target(1)==0&&target(end)==g.r(end)&&all(diff(target)>0), ...
        'ipm_c:AreaRemeshAxis','A proposal must retain the domain, node counts and origin.');
    [target,spacingFloor]=ipm_c.mesh.firstSpacingFloor(target,option(policy,'minimumFirstSpacing',0));
    if spacingFloor.applied
        originalProposal=proposal;
        if strcmp(originalProposal.kind,'positive_mirrored_gaussian_radial_cdf_v1')
            proposal=struct('kind','positive_radial_cdf_with_smooth_odd_first_spacing_floor_v1', ...
                'cdfBeforeFloor',originalProposal,'firstSpacingFloor',spacingFloor, ...
                'firstSpacing',target(2),'exactCdfQuantilesClaimed',false);
        elseif strcmp(originalProposal.kind,'complete_source_feature_uniform_core_axis_v1')
            proposal=originalProposal;
            proposal.kind='complete_source_feature_uniform_core_with_smooth_odd_first_spacing_floor_v1';
            proposal.uniformCoreBeforeFloor=originalProposal;
            proposal.uniformCoreInfoAppliesTo='axis_before_first_spacing_projection';
            proposal.firstSpacingFloor=spacingFloor;proposal.firstSpacing=target(2);
        else
            proposal=struct('kind','manufactured_axis_with_smooth_odd_first_spacing_floor_v1', ...
                'axisBeforeFloor',originalProposal,'firstSpacingFloor',spacingFloor, ...
                'firstSpacing',target(2),'exactCdfQuantilesClaimed',false);
        end
    else
        proposal.firstSpacingFloor=spacingFloor;
    end
catch problem
    if ~recoverable(problem.identifier),rethrow(problem);end
    [state,info]=reject(accepted,info,features,problem.identifier,problem.message);return;
end
info.proposal=proposal;oldMass=sum((state.rho+state.scale.densityOrigin).*ops.integrationWeights,'all');info.oldPolarMass=oldMass;
oldMeasure=ipm_c.gauge.measure(state.rho,g.r,ops.cart.originalWallOps,struct('interpolation','spline'));
minimumImprovement=option(policy,'minimumImprovement',.05);
reason='ipm_c:AreaRemeshNoChange';detail='The proposal does not move the accepted radial mesh.';
maximumRetries=12;
if manufactured&&isfield(fixture,'maximumRetries')
    maximumRetries=fixture.maximumRetries;
    assert(isscalar(maximumRetries)&&isfinite(maximumRetries)&&maximumRetries>=0&&maximumRetries<=12&& ...
        maximumRetries==fix(maximumRetries),'ipm_c:RemeshFixture');
end
for retry=0:maximumRetries
    damping=2^-retry;
    if damping<minimumImprovement,break;end
    r=g.r+damping*(target-g.r);r([1,end])=g.r([1,end]);
    identity=isequaln(r,g.r);
    if identity&&~(manufactured&&option(fixture,'allowIdentity',false)),break;end
    row=struct('retry',retry,'damping',damping,'passed',false,'identifier','','detail','', ...
        'oldPolarMass',oldMass,'newPolarMass',NaN,'massRelativeJump',NaN, ...
        'wallIntegralRelativeJump',NaN,'wallPeakRelativeJump',NaN, ...
        'poissonLinearBackwardError',NaN,'gaugeResidual',NaN, ...
        'poissonStrongResidual',NaN,'velocityTransferDefect',NaN,'qRateTransferDefect',NaN);
    try
        gaps=diff(r);ratio=max([gaps(2:end)./gaps(1:end-1);gaps(1:end-1)./gaps(2:end)]);
        assert(all(gaps>0)&&isfinite(ratio)&&ratio<=option(policy,'maximumCellRatio',6)*(1+128*eps), ...
            'ipm_c:AreaRemeshAxis','The damped axis violates the adjacent-cell ratio contract.');
        candidateGrid=ipm_c.remesh.radialGeometry(g,r);
        candidateOps=ipm_c.spatial.build(ops.polarCfg,ops.cart,candidateGrid);
        hooks=fieldnames(ops.spatialOperators);
        for hook=1:numel(hooks)
            key=hooks{hook};
            if ~isfield(candidateOps.spatialOperators,key),candidateOps.spatialOperators.(key)=ops.spatialOperators.(key);end
        end
        candidateOps.rescaling=ops.rescaling;candidateOps.remeshCount=info.count;
        candidateRho=transfer_centered(g.r,state.rho,r);
        candidateRho(1,:)=state.rho(1,1);candidateRho(end,:)=state.rho(end,:);
        assert(all(isfinite(candidateRho),'all')&&all(candidateRho(1,:)==candidateRho(1,1)), ...
            'ipm_c:AreaRemeshTransfer','The full transferred field must be finite and retain one physical origin.');
        newMass=sum((candidateRho+state.scale.densityOrigin).*candidateOps.integrationWeights,'all');
        row.newPolarMass=newMass;info.newPolarMass=newMass;
        row.massRelativeJump=abs(newMass-oldMass)/max(abs(oldMass),realmin);
        assert(isfinite(row.massRelativeJump)&&row.massRelativeJump<=option(policy,'massTransferTolerance',1e-3), ...
            'ipm_c:AreaRemeshMass','The polar-weighted total mass jump exceeds the transfer contract.');
        candidateMeasure=ipm_c.gauge.measure(candidateRho,r,ops.cart.originalWallOps,struct('interpolation','spline'));
        row.wallIntegralRelativeJump=abs(candidateMeasure.I-oldMeasure.I)/max(abs(oldMeasure.I),realmin);
        assert(row.wallIntegralRelativeJump<=option(policy,'wallIntegralTransferTolerance',1e-3), ...
            'ipm_c:AreaRemeshWallIntegral','The fixed original-wall integral jump exceeds the transfer contract.');
        row.wallPeakRelativeJump=abs(candidateMeasure.P-oldMeasure.P)/max(abs(oldMeasure.P),realmin);
        assert(isfinite(row.wallPeakRelativeJump)&& ...
            row.wallPeakRelativeJump<=option(policy,'peakTransferTolerance',1e-4), ...
            'ipm_c:AreaRemeshWallPeak','The fixed original-wall gradient peak jump exceeds the transfer contract.');
        [candidateRate,candidateFlow]=ipm.evolve.flow(candidateRho,candidateOps,state.scale);
        audit=candidateFlow.ipmCInfo;
        row.poissonLinearBackwardError=audit.poisson.linearSolveBackwardError;
        row.gaugeResidual=max(abs(audit.gauge.auditResidual));
        row.poissonStrongResidual=audit.poisson.maximumOffGridStrongResidualScaledBySource;
        assert(all(isfinite(candidateRate),'all')&&all(isfinite(candidateFlow.source),'all')&& ...
            all(isfinite(candidateFlow.fullQ),'all')&&isfinite(row.poissonLinearBackwardError)&& ...
            isfinite(row.gaugeResidual)&&row.poissonLinearBackwardError<=1e-10&& ...
            audit.gauge.auditTolerance<=1e-5&&row.gaugeResidual<=audit.gauge.auditTolerance, ...
            'ipm_c:AreaRemeshFullAudit','The transferred field must pass the complete fresh physical RHS and original I/P/C gates.');
        % Compare complete fields at the SAME new coordinates. Equal mass and
        % an accurate linear solve cannot detect a corrupted Poisson load.
        predicted=transfer_centered(g.r,cat(2,state.flow.u1,state.flow.u2,state.flow.fullQ),r);
        nt=g.ntheta;
        row.velocityTransferDefect=max(hypot(candidateFlow.u1-predicted(:,1:nt), ...
            candidateFlow.u2-predicted(:,nt+(1:nt))),[],'all')/ ...
            max(1,max(hypot(state.flow.u1,state.flow.u2),[],'all'));
        row.qRateTransferDefect=max(abs(candidateFlow.fullQ-predicted(:,2*nt+(1:nt))),[],'all')/ ...
            max(1,max(abs(state.flow.fullQ),[],'all'));
        assert(row.velocityTransferDefect<=option(policy,'velocityTransferTolerance',.05)&& ...
            row.qRateTransferDefect<=option(policy,'qRateTransferTolerance',.1), ...
            'ipm_c:AreaRemeshFieldJump','The transferred velocity or complete Q changes beyond the spatial transfer limit.');
        if manufactured&&option(fixture,'rejectAfterFullAudit',false)
            error('ipm_c:AreaRemeshManufacturedRejection','Explicit rollback fixture after a successful complete candidate audit.');
        end
        candidateFeatures=ipm_c.diagnostics.sourceFeatures(candidateFlow.source,candidateGrid,levels);
        assert(candidateFeatures.usableForCdf,'ipm_c:AreaRemeshGeometry','The complete transferred source geometry must remain positive.');
        assert(isequaln(candidateOps.cart.originalWallOps,ops.cart.originalWallOps), ...
            'ipm_c:AreaRemeshReferences','Original wall operators and I/P/C references must remain fixed.');
        row.passed=true;info.attempts=append(info.attempts,row);
        candidateOps.remeshCount=info.count+1;
        candidateOps.rescaling.referenceLevelSetAreas=candidateFeatures.areas;
        candidateOps.rescaling.referenceLevelSetRadialBounds=candidateFeatures.radialProjectionBounds;
        observer.failures=0;observer.retryAt=0;observer.features=[];
        candidateOps.rescaling.areaRemeshObserver=observer;
        candidateOps.rescaling.gaugeRateGuess=[candidateFlow.c_l,candidateFlow.c_r];
        state.rho=candidateRho;state.ops=candidateOps;state.flow=candidateFlow;
        state.rhsCache=ipm.evolve.makeRhsCache(candidateRho,candidateRate,candidateFlow,candidateOps,state.scale);
        info.applied=true;info.pending=false;info.status='applied';info.count=candidateOps.remeshCount;
        info.damping=damping;info.areas=candidateFeatures.areas;info.referenceAreas=candidateFeatures.areas;
        info.oldPolarMass=oldMass;info.newPolarMass=newMass;info.geometryChanged=~identity;
        if isfield(state.config,'output')&&option(state.config.output,'verbose',false)
            fprintf(['ipm_c REMESH tau=%.12g step=%d remeshCount=%d ' ...
                'I=%.17g->%.17g P=%.17g->%.17g peakJump=%.6g strong=%.6g\n'], ...
                state.scale.canonicalTime,state.step,info.count, ...
                oldMeasure.I,candidateMeasure.I,oldMeasure.P,candidateMeasure.P, ...
                row.wallPeakRelativeJump,row.poissonStrongResidual);
        end
        return;
    catch problem
        if ~recoverable(problem.identifier),rethrow(problem);end
        reason=problem.identifier;detail=problem.message;row.identifier=reason;row.detail=detail;
        info.attempts=append(info.attempts,row);
    end
end
[state,info]=reject(accepted,info,features,reason,detail);
end

function [state,info]=reject(state,info,features,identifier,message)
observer=state.ops.rescaling.areaRemeshObserver;observer.failures=min(observer.failures+1,5);
delay=8*2^(observer.failures-1);observer.retryAt=observer.checks+delay;
observer.features=struct('areas',features.areas,'radialProjectionBounds',features.radialProjectionBounds);
state.ops.rescaling.areaRemeshObserver=observer;
info.status='candidate_rejected';info.pending=true;info.rejectionReason=identifier;
info.rejectionDetail=strtrim(regexprep(message,'[\r\n]+',' '));info.retryAfterChecks=delay;
end
function change=feature_change(current,previous)
change=Inf;if isempty(previous),return;end
if all(previous.areas>0)&&all(isfinite(current.radialProjectionBounds),'all')&& ...
        all(isfinite(previous.radialProjectionBounds),'all')
    width=max(diff(previous.radialProjectionBounds,1,2),realmin);
    change=max([max(abs(current.areas./previous.areas-1)); ...
        max(abs(current.radialProjectionBounds-previous.radialProjectionBounds)./width,[],'all')]);
end
end
function yes=recoverable(identifier)
yes=startsWith(identifier,'ipm_c:AreaRemesh')||startsWith(identifier,'ipm_c:Cdf')|| ...
    startsWith(identifier,'polar:UniformCore')|| ...
    any(strcmp(identifier,{'polar:GaugeRootUnresolved','polar:GaugeFullRhsAudit','polar:WallDegenerate', ...
    'polar:BSplineStrongResidual', ...
    'polar:NonfiniteRhs','polar:GridMetric','polar:OriginFitCondition','ipm:HighOrderGridMetric', ...
    'ipm:HighOrderInterpolationGrid','ipm:FiniteDifferenceScale','ipm:UnresolvedWallIntegralGauge', ...
    'ipm:DegenerateWallIntegralGauge'}));
end
function result=append(result,row)
if isempty(result),result=row;else,result(end+1)=row;end
end
function result=option(options,key,fallback)
if isfield(options,key),result=options.(key);else,result=fallback;end
end

function radius=selected_core_radius(policy,settings,features)
outer=max(features.featureBounds(:,2));
peak=option(policy,'uniformCorePadding',2)*features.featureBounds(end,2);
radius=max(settings.coreRadius,min(outer,peak));
end

function result=transfer_centered(oldAxis,values,newAxis)
% Keep the shared degree-five polynomial, but remove each stencil's local
% reference before evaluation so a large constant platform is exact.
oldAxis=oldAxis(:);newAxis=newAxis(:);count=numel(oldAxis);
assert(count>=6,'ipm:HighOrderInterpolationGrid');
cellIndex=discretize(newAxis,oldAxis);cellIndex(newAxis==oldAxis(end))=count-1;
assert(all(isfinite(cellIndex)),'ipm_c:AreaRemeshTransfer', ...
    'The remesh transfer may not extrapolate beyond the accepted domain.');
first=min(max(cellIndex-2,1),count-5);result=zeros(numel(newAxis),size(values,2));
for start=unique(first).'
    stencil=start:start+5;queries=find(first==start);reference=values(start+2,:);
    result(queries,:)=ipm.remesh.interpolate(oldAxis(stencil), ...
        values(stencil,:)-reference,newAxis(queries), ...
        struct('sampleDimension',1,'conservation','none'))+reference;
end
end
