function report=test_ipm_cr_fourth_order(outputDir,options)
%TEST_IPM_CR_FOURTH_ORDER Spatial manufacture and complete-IPM RK4 evidence.
% Pure transport and production advance share time.rk4. The
% temporal experiment calls the actual production advance, including its
% Poisson solve, original gauge, scalar scales and compensated clock.
% Every run remains within the complete physical CFL. Classic RK4 has no
% SSP guarantee; neither the tableau nor a CFL check alone proves stability.
if nargin<2,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir)
    outputDir=fullfile(root,'results','fourth_order');
end
o=resolve(options);outputDir=char(outputDir);
assert(startsWith(outputDir,filesep),'ipm_cr:FourthOutput','Use an absolute output directory.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'fourth_order.json')), ...
    'ipm_cr:FourthPreserve','Choose a fresh output directory.');
timer=tic;
tableau=ipm_cr.time.rk4();
assert(tableau.order==4&&tableau.stages==4,'ipm_cr:FourthTableau', ...
    'This qualification requires the shared classic four-stage RK4 tableau.');
report=struct('kind','hybrid_fourth_order_spatial_and_full_IPM_time_v1', ...
    'status','running','options',o,'tableau',tableau, ...
    'timeIntegrator','classic_rk4','sspGuaranteeClaimed',false, ...
    'configuredPhysicalCfl',.3,'allActiveOverlapInteriorsIncluded',true, ...
    'sourceInventory',source_inventory(),'sourcesUnchanged',false, ...
    'spatial',struct(),'poissonChain',struct(),'temporal',struct(), ...
    'globalFourthOrderQualified',false,'tau200AccuracyQualified',false, ...
    'elapsedSeconds',0,'matlabVersion',version,'computer',computer, ...
    'limitations',{{'Manufactured smooth spatial convergence is not an error certificate for every IPM state or exterior source tail.', ...
    'Temporal self-convergence uses nested nominal maximum steps and resolves original spline-peak knot events; the same complete spatial operator and original IPC definition are held fixed.', ...
    'Classic RK4 is not SSP. Complete stage CFL checks are necessary operational constraints, not a theorem of nonlinear overset stability.', ...
    'A short qualification does not certify the entire 200-tau trajectory or its 50-minute time budget.'}});
persist(outputDir,report);
try
    if o.runSpatial
        report.spatial=spatial_test(o,tableau,outputDir,timer,report);
        report.elapsedSeconds=toc(timer);persist(outputDir,report);
    end
    if o.runPoisson
        report.poissonChain=poisson_test(o);
        report.elapsedSeconds=toc(timer);persist(outputDir,report);
    end
    if o.runTemporal
        report.temporal=temporal_test(o,outputDir,timer,report);
        report.elapsedSeconds=toc(timer);persist(outputDir,report);
    end
    report.sourcesUnchanged=isequaln(report.sourceInventory,source_inventory());
    assert(report.sourcesUnchanged,'ipm_cr:FourthSources','Recorded numerical sources changed during qualification.');
    report.globalFourthOrderQualified=o.runSpatial&&o.runPoisson&&o.runTemporal&& ...
        report.spatial.observedAtLeastFourth&&report.spatial.timeContributionScreenPassed&& ...
        report.poissonChain.observedAtLeastFourth&&report.temporal.fieldFourthOrderObserved&& ...
        report.temporal.coupledScalarFourthOrderObserved&&report.temporal.allScalarFourthOrderObserved;
    if report.globalFourthOrderQualified,report.status='completed_requested_fourth_order_screens';
    else,report.status='completed_measurements_qualification_incomplete';end
catch problem
    % Retain the latest completed local levels if a later run fails.
    savedPath=fullfile(outputDir,'fourth_order.mat');
    if isfile(savedPath)
        saved=load(savedPath,'report');report=saved.report;
    end
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function o=resolve(options)
o=struct('runSpatial',true,'runPoisson',true,'runTemporal',true, ...
    'spatialFinalTau',.6,'spatialStepRateProduct',.03,'spatialHalfStepCheck',true, ...
    'timeFinalTau',.4,'timeBaseSteps',20,'timeLevels',4, ...
    'timeExplicitOptions',struct('geometry',struct('coreSpacing',.06, ...
        'grid',struct('nr',193,'ntheta',64,'firstSpacing',.06))), ...
    'poissonSolverOptions',struct( ...
        'coreDegree',[64,64],'angularDegree',80,'polarDegree',20, ...
        'rOuter',1e10,'requireAccuracy',true));
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(o))), ...
    'ipm_cr:FourthOptions','Unknown fourth-order qualification option.');
names=fieldnames(options);for j=1:numel(names),o.(names{j})=options.(names{j});end
for name={'runSpatial','runPoisson','runTemporal','spatialHalfStepCheck'}
    validateattributes(o.(name{1}),{'logical'},{'scalar'});
end
validateattributes(o.spatialFinalTau,{'numeric'},{'scalar','finite','>=',.6});
validateattributes(o.spatialStepRateProduct,{'numeric'},{'scalar','finite','positive','<=',.05});
validateattributes(o.timeFinalTau,{'numeric'},{'scalar','finite','positive'});
validateattributes(o.timeBaseSteps,{'numeric'},{'scalar','integer','>=',1});
validateattributes(o.timeLevels,{'numeric'},{'scalar','integer','>=',3,'<=',5});
assert(isstruct(o.timeExplicitOptions)&&isscalar(o.timeExplicitOptions)&& ...
    isstruct(o.poissonSolverOptions)&&isscalar(o.poissonSolverOptions),'ipm_cr:FourthOptions');
end

function result=spatial_test(o,tableau,folder,timer,report)
plan=grid_plan();rows=struct([]);queries=common_queries();queries.endpoint=o.spatialFinalTau;
result=struct('kind','smooth_nonradial_translated_Gaussians_RK4_spatial', ...
    'manufacture',manufacture(),'targetTau',o.spatialFinalTau, ...
    'configuredCfl',.3,'additionalMaximumDtRateProduct',o.spatialStepRateProduct, ...
    'tableauSharedWithProduction',true,'activeValuesProjected',false, ...
    'allFiveStageAndEndpointCflChecks',true,'records',rows, ...
    'convergence',struct([]),'observedAtLeastFourth',false, ...
    'timeContributionScreenPassed',false, ...
    'regions',struct('allActive','Both complete native active domains including both copies of their overlap; no far-field mask.', ...
        'fixedCommonNear','Identical fractionally offset physical queries in r<=3 plus physical wall points.', ...
        'interfaceBand','Three active stencil layers beside the artificial ghost collars on each native patch.', ...
        'physicalWall','All active native nodes with Y=0, including the full exterior wall.'));
result.radialIntervals=[plan.nr]-1;
result.radialIntervalRefinementRatios=result.radialIntervals(2:end)./result.radialIntervals(1:end-1);
result.angularIntervals=[plan.ntheta]-1;
result.angularIntervalRefinementRatios=result.angularIntervals(2:end)./result.angularIntervals(1:end-1);
for j=1:numel(plan)
    g=make_geometry(plan(j));f=sample_transport(g,0);v=velocity(g);
    [~,rate]=transport_rhs(f,g,v);
    n=ceil(o.spatialFinalTau*rate/o.spatialStepRateProduct);
    [final,row]=transport_run(f,g,v,o.spatialFinalTau,n,tableau,queries);
    row.level=plan(j);row.ghostCrossings=crossings(g,o.spatialFinalTau);
    row.actualSpacing=spacing_summary(g);
    assert(row.ghostCrossings.annularPeakCrosses&&row.ghostCrossings.corePeakCrosses&& ...
        row.ghostCrossings.wallPeakCrosses,'ipm_cr:FourthCrossing', ...
        'The broad manufactured peaks must cross both artificial interfaces and the wall interface.');
    if o.spatialHalfStepCheck&&(j==1||j==numel(plan))
        [half,halfRow]=transport_run(f,g,v,o.spatialFinalTau,2*n,tableau,queries);
        difference=transport_errors(field_subtract(final,half),g,queries,[],true);
        row.halfStep=struct('dt',halfRow.dt,'steps',halfRow.steps, ...
            'maximumStageCfl',halfRow.maximumStageCfl,'fieldDifference',difference, ...
            'differenceOverExactSpatialError',difference.allActiveInfinity/max(row.errors.allActiveInfinity,realmin), ...
            'thresholdFraction',.1,'screenPassed',difference.allActiveInfinity<=.1*row.errors.allActiveInfinity);
    else,row.halfStep=struct('screenPerformed',false);end
    rows=append(rows,row);result.records=rows;
    report.spatial=result;report.elapsedSeconds=toc(timer);persist(folder,report);
    fprintf('FOURTH SPACE h=%.4g steps=%d CFL=%.6g full/near/band/wall %.3e %.3e %.3e %.3e\n', ...
        plan(j).h,n,row.maximumStageCfl,row.errors.allActiveInfinity, ...
        row.errors.fixedCommonNearInfinity,row.errors.interfaceBandInfinity,row.errors.physicalWallInfinity);
end
keys={'allActiveInfinity','fixedCommonNearInfinity','interfaceBandInfinity','physicalWallInfinity'};
allPass=true;
for name=keys
    errors=arrayfun(@(r)r.errors.(name{1}),rows);
    orders=log2(errors(1:end-1)./errors(2:end));
    item=struct('quantity',name{1},'h',[plan.h],'errors',errors,'orders',orders, ...
        'finestObservedOrder',orders(end),'atLeastFourthObserved',orders(end)>=4);
    allPass=allPass&&item.atLeastFourthObserved;result.convergence=append(result.convergence,item);
end
result.observedAtLeastFourth=allPass;
result.timeContributionScreenPassed=o.spatialHalfStepCheck&&rows(1).halfStep.screenPassed&& ...
    rows(end).halfStep.screenPassed;
end

function [rho,row]=transport_run(initial,g,v,T,n,tableau,queries)
timer=tic;rho=initial;dt=T/n;maximumCfl=0;maximumActiveSyncChange=0;
for step=1:n
    stages=cell(1,tableau.stages);
    for k=1:tableau.stages
        stage=field_combination(rho,stages,dt*tableau.A(k,:),g);
        [stages{k},rate,change]=transport_rhs(stage,g,v);
        maximumCfl=max(maximumCfl,dt*rate);maximumActiveSyncChange=max(maximumActiveSyncChange,change);
    end
    rho=field_combination(rho,stages,dt*tableau.b(:)',g);
    [~,endpointRate,change,rho]=transport_rhs(rho,g,v);
    maximumCfl=max(maximumCfl,dt*endpointRate);
    maximumActiveSyncChange=max(maximumActiveSyncChange,change);
    assert(maximumCfl<=.3*(1+128*eps),'ipm_cr:FourthTransportCfl', ...
        'The shared RK4 transport stages exceeded the full physical CFL.');
    if mod(step,500)==0,fprintf('FOURTH SPACE progress %d/%d tau %.6g\n',step,n,T*step/n);end
end
exact=sample_transport(g,T);
row=struct('tau',T,'steps',n,'dt',dt,'seconds',toc(timer), ...
    'maximumStageCfl',maximumCfl,'maximumActiveChangeDuringSynchronization',maximumActiveSyncChange, ...
    'errors',transport_errors(rho,g,queries,exact,false), ...
    'method','same_classic_RK4_tableau_as_production_advance','temporalErrorRemovedByExtrapolation',false);
end

function [F,rate,change,synced]=transport_rhs(input,g,v)
synced=ipm_cr.mesh.fill(input,g);
change=max([max(abs(synced.core(g.core.active)-input.core(g.core.active))), ...
    max(abs(synced.polar(g.polar.active)-input.polar(g.polar.active)))]);
assert(change==0,'ipm_cr:FourthActiveProjection','Synchronization changed an active state.');
value=ipm_cr.op.apply(ipm_cr.op.prepare(synced,g,struct()),v.x,v.y);
F=value.B;rate=value.rate;
end

function result=transport_errors(value,g,q,exact,isDifference)
if isDifference,d=value;else,d=field_subtract(value,exact);end
result=struct('patches',struct());combined=0;bandMaximum=0;wallMaximum=0;
for patch={'core','polar'}
    key=patch{1};p=g.(key);band=interface_band(g,key);
    row=struct('allActive',metric(d.(key),p.active,p.X,p.Y), ...
        'interfaceBand',metric(d.(key),band,p.X,p.Y), ...
        'physicalWall',metric(d.(key),p.active&p.Y==0,p.X,p.Y));
    result.patches.(key)=row;combined=max(combined,row.allActive.absoluteInfinity);
    bandMaximum=max(bandMaximum,row.interfaceBand.absoluteInfinity);
    wallMaximum=max(wallMaximum,row.physicalWall.absoluteInfinity);
end
if isDifference
    sampled=ipm_cr.mesh.sample(value,g,q.x,q.y);difference=sampled;
else
    sampled=ipm_cr.mesh.sample(value,g,q.x,q.y);
    difference=sampled-transport_field(q.x,q.y,manufacture(),q.endpoint);
end
result.allActiveInfinity=combined;result.interfaceBandInfinity=bandMaximum;
result.physicalWallInfinity=wallMaximum;
result.fixedCommonNearInfinity=max(abs(difference));
result.fixedCommonNearRms=sqrt(mean(difference.^2));
[~,index]=max(abs(difference));result.fixedCommonNearMaximumLocation=struct('X',q.x(index),'Y',q.y(index));
end

function f=sample_transport(g,t)
m=manufacture();f=struct('core',transport_field(g.core.X,g.core.Y,m,t), ...
    'polar',transport_field(g.polar.X,g.polar.Y,m,t));
end

function f=transport_field(x,y,m,t)
f=zeros(size(x));
for j=1:numel(m.amplitudes)
    f=f+m.amplitudes(j)*exp(-((x-m.centers(j,1)-m.velocity(1)*t)/m.widths(j,1)).^2- ...
        ((y-m.centers(j,2))/m.widths(j,2)).^2);
end
end

function m=manufacture()
m=struct('velocity',[.6,0],'amplitudes',[1,.55,.4], ...
    'centers',[.62,.35;1,.70;.95,0],'widths',[.45,.50;.50,.45;.55,.50], ...
    'description','Three broad nonradial translated Gaussians; third peak moves along the physical wall. No compact narrow bump or density-origin resetting is used.');
end

function v=velocity(g)
v=struct('x',struct('core',.6*ones(size(g.core.X)),'polar',.6*ones(size(g.polar.X))), ...
    'y',struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X))));
end

function c=crossings(g,T)
c=struct('annularPeakCrosses',hypot(.62,.35)<g.polar.innerActiveRadius&& ...
    hypot(.62+.6*T,.35)>g.polar.innerActiveRadius, ...
    'corePeakCrosses',1<g.core.activeXBounds(end)&&1+.6*T>g.core.activeXBounds(end), ...
    'wallPeakCrosses',.95<g.core.activeXBounds(end)&&.95+.6*T>g.core.activeXBounds(end));
end

function result=poisson_test(o)
plan=grid_plan();rows=struct([]);
result=struct('kind','analytic_decaying_W_to_q_to_u_chain', ...
    'definition','psi=x*y*exp(-r^2); W=2*y*(r^2-2)*exp(-r^2); q=4*x*y*(3-r^2)*exp(-r^2)', ...
    'velocityDefinition','u1=-x*(1-2*y^2)*exp(-r^2); u2=y*(1-2*x^2)*exp(-r^2)', ...
    'samePoissonDegreesForAllNativeLevels',true,'solverOptions',o.poissonSolverOptions, ...
    'fullExteriorDomainIncluded',true,'records',rows,'convergence',struct([]), ...
    'rawNullspaceField','W=exp(-y^2), exact q=0; full native domain retained', ...
    'rawNullspaceIncludedInNewMethodQualification',false, ...
    'analyticNullspaceModel','Actual gaussian_primitive provider with amplitude=1/sqrt(pi), widths=1, skew=0, centerX=-2e10; G=exp(-Y^2), H=expm1(-Y^2), Gx=0 on R<=1e10 to floating-point accuracy.', ...
    'observedAtLeastFourth',false);
result.radialIntervals=[plan.nr]-1;
result.radialIntervalRefinementRatios=result.radialIntervals(2:end)./result.radialIntervals(1:end-1);
result.angularIntervals=[plan.ntheta]-1;
result.angularIntervalRefinementRatios=result.angularIntervals(2:end)./result.angularIntervals(1:end-1);
for j=1:numel(plan)
    timer=tic;g=make_geometry(plan(j));W=struct();exactq=struct();exactu=struct();null=struct();
    for patch={'core','polar'}
        key=patch{1};p=g.(key);r2=p.X.^2+p.Y.^2;e=exp(-r2);
        W.(key)=2*p.Y.*(r2-2).*e;exactq.(key)=4*p.X.*p.Y.*(3-r2).*e;
        exactu.(key)=struct('u1',-p.X.*(1-2*p.Y.^2).*e,'u2',p.Y.*(1-2*p.X.^2).*e);
        null.(key)=exp(-p.Y.^2);
    end
    W=ipm_cr.mesh.fill(W,g);d=ipm_cr.mesh.diff(W,g);
    q=ipm_cr.mesh.fill(struct('core',d.core.x,'polar',d.polar.x),g);
    planPoisson=ipm_cr.poisson.plan(g,o.poissonSolverOptions);
    [flow,diagnostics]=ipm_cr.poisson.apply(q,g,planPoisson,[],true,true);
    [sampledSourceFlow,sampledSourceDiagnostics]=ipm_cr.poisson.apply(exactq,g,planPoisson);
    n=ipm_cr.mesh.diff(ipm_cr.mesh.fill(null,g),g);
    nullq=struct('core',n.core.x,'polar',n.polar.x);
    zero=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
    plateauModel=struct('kind','gaussian_primitive','amplitude',1/sqrt(pi), ...
        'widthX',1,'widthY',1,'skewWeight',0,'centerX',-2e10);
    plateauScale=struct('logC_l',0,'logC_omega',0,'X_shift',0);
    plateauSource=@(X,Y)background_source(plateauModel,plateauScale,X,Y);
    [plateauFlow,plateauInfo]=ipm_cr.poisson.apply(zero,g,planPoisson,plateauSource);
    plateauQ=struct();plateauHError=struct();zeroVelocity=struct();
    for patch={'core','polar'}
        key=patch{1};p=g.(key);bg=ipm_cr.model.background(plateauModel,plateauScale,p.X,p.Y);
        plateauQ.(key)=bg.x;plateauHError.(key)=bg.H-expm1(-p.Y.^2);
        zeroVelocity.(key)=struct('u1',zeros(size(p.X)),'u2',zeros(size(p.X)));
    end
    [analyticFlow,analyticInfo]=ipm_cr.poisson.apply(zero,g,planPoisson,@manufactured_source);
    qError=native_field_metrics(field_subtract(q,exactq),g);
    velocityError=velocity_metrics(flow,exactu,g);
    transferSolveError=velocity_metrics(sampledSourceFlow,exactu,g);
    row=struct('level',plan(j),'qError',qError,'velocityError',velocityError, ...
        'actualSpacing',spacing_summary(g), ...
        'analyticSourceSamplingSolveReturnError',transferSolveError, ...
        'rawDerivativeNullspaceSourceError',native_field_metrics(nullq,g), ...
        'rawNullspaceDiagnosticOnly',true, ...
        'analyticBackgroundNullspaceSourceError',native_field_metrics(plateauQ,g), ...
        'analyticBackgroundNullspaceVelocityError',velocity_metrics(plateauFlow,zeroVelocity,g), ...
        'analyticBackgroundPlateauHError',native_field_metrics(plateauHError,g), ...
        'analyticBackgroundNonzeroSourceVelocityError',velocity_metrics(analyticFlow,exactu,g), ...
        'analyticBackgroundNullspacePoisson',plateauInfo,'analyticBackgroundNonzeroPoisson',analyticInfo, ...
        'poisson',diagnostics,'analyticSampledSourcePoisson',sampledSourceDiagnostics, ...
        'seconds',toc(timer));
    rows=append(rows,row);
    fprintf('FOURTH POISSON h=%.4g q/u/null-global %.3e %.3e %.3e\n', ...
        plan(j).h,qError.allActiveInfinity,velocityError.allActiveInfinity,row.rawDerivativeNullspaceSourceError.allActiveInfinity);
end
result.records=rows;
for name={'qError','velocityError','analyticBackgroundNullspaceSourceError','analyticBackgroundNullspaceVelocityError'}
    key=name{1};errors=arrayfun(@(r)r.(key).allActiveInfinity,rows);
    orders=log2(errors(1:end-1)./errors(2:end));
    floorValue=1024*eps;
    passed=orders(end)>=4||all(errors<=floorValue);
    result.convergence=append(result.convergence,struct('quantity',key,'errors',errors,'orders',orders, ...
        'atLeastFourthObservedOrAtRoundoff',passed,'roundoffThreshold',floorValue));
end
result.observedAtLeastFourth=all([result.convergence.atLeastFourthObservedOrAtRoundoff]);
result.nonzeroAnalyticBackgroundSourcePathPassed=all(arrayfun(@(r) ...
    r.analyticBackgroundNonzeroSourceVelocityError.allActiveInfinity<=1e-8,rows));
result.observedAtLeastFourth=result.observedAtLeastFourth&&result.nonzeroAnalyticBackgroundSourcePathPassed;
end

function q=background_source(model,scale,X,Y)
b=ipm_cr.model.background(model,scale,X,Y);q=b.x;
end

function q=manufactured_source(X,Y)
r2=X.^2+Y.^2;q=4*X.*Y.*(3-r2).*exp(-r2);
end

function result=temporal_test(o,folder,timer,report)
method=o.timeExplicitOptions;method.timeIntegrator='rk4';method.cfl=.3;
[initial,ctx]=ipm_cr.init(method);g=ctx.mesh;wallOps=ctx.wall;method=ctx.opts;
rows=struct([]);states=cell(1,o.timeLevels);
result=struct('kind','complete_IPM_production_RK4_fixed_grid_self_convergence', ...
    'targetTau',o.timeFinalTau,'baseSteps',o.timeBaseSteps,'levels',o.timeLevels, ...
    'sameInitialStateAndWallAuthority',true,'spaceGridHeldFixed',true, ...
    'PoissonGaugeTransportScalarsAndClockAllEnabled',true,'timeIntegrator','rk4', ...
    'explicitOptions',method,'records',rows,'differences',struct([]),'orders',struct([]), ...
    'fieldFourthOrderObserved',false,'coupledScalarFourthOrderObserved',false, ...
    'allScalarFourthOrderObserved',false, ...
    'fourthOrderWindow',[3.5,4.5],'roundoffSeparationFactor',1024, ...
    'gaugeScope','Initial and accepted endpoint residuals; stage scope is additionally recorded when production advance exposes it.');
for level=1:o.timeLevels
    n=o.timeBaseSteps*2^(level-1);dt=o.timeFinalTau/n;state=initial;runTimer=tic;
    method.maxDt=dt;maximumCfl=0;maximumGauge=max(abs(initial.rhs.gaugeResidual));
    productionMethods={};peakStageRecords=struct([]);eventRecords=struct([]);
    actualSteps=0;minimumActualDt=dt;
    for step=1:n
        target=initial.scale.canonicalTime+o.timeFinalTau*step/n;
        while target-state.scale.canonicalTime>32*eps(max(1,abs(target)))
            remaining=target-state.scale.canonicalTime;
            [next,info]=ipm_cr.time.step(state,g,wallOps,method,target,ctx);
            actual=next.scale.canonicalTime-state.scale.canonicalTime;
            shortened=actual<remaining-256*eps(max(1,o.timeFinalTau));
            eventLimited=isfield(info,'peakEvent')&&info.peakEvent.limited;
            assert(actual>0&&actual<=dt+256*eps(max(1,o.timeFinalTau))&& ...
                (~shortened||eventLimited), ...
                'ipm_cr:FourthTimeGrid', ...
                'Only original-peak event splits may subdivide a nominal interval; CFL reduction requires a finer common plan.');
            assert(contains(lower(info.method),'rk4'),'ipm_cr:FourthProductionMethod', ...
                'The actual production advance did not report RK4.');
            actualSteps=actualSteps+1;minimumActualDt=min(minimumActualDt,actual);
            maximumCfl=max(maximumCfl,info.realizedMaximumStageCfl);
            maximumGauge=max(maximumGauge,max(abs(next.rhs.gaugeResidual)));
            if isfield(info,'maximumStageGaugeResidual')
                maximumGauge=max(maximumGauge,info.maximumStageGaugeResidual);
                result.gaugeScope='All production RK4 stages and the accepted endpoint, reported by advance.maximumStageGaugeResidual.';
            end
            if eventLimited
                event=info.peakEvent;event.nominalStep=step;event.actualStep=actualSteps;
                eventRecords=append(eventRecords,event);
            end
            if isfield(info,'allAttemptPeakDiagnostics')
                snapshots=info.allAttemptPeakDiagnostics;
                for k=1:numel(snapshots),snapshots(k).step=actualSteps;end
                if isempty(peakStageRecords),peakStageRecords=snapshots;
                else,peakStageRecords=[peakStageRecords,snapshots];end %#ok<AGROW>
            end
            state=next;if isempty(productionMethods),productionMethods={info.method};end
        end
        if mod(step,100)==0,fprintf('FOURTH TIME level=%d nominal=%d/%d actual=%d tau=%.6g\n',level,step,n,actualSteps,state.scale.canonicalTime);end
    end
    assert(maximumCfl<=.3*(1+128*eps),'ipm_cr:FourthTimeCfl');
    states{level}=state;
    row=struct('steps',actualSteps,'nominalIntervals',n,'dt',dt,'minimumActualDt',minimumActualDt, ...
        'originalPeakEventSplits',eventRecords,'actualTau',state.scale.canonicalTime, ...
        'seconds',toc(runTimer),'maximumStageCfl',maximumCfl,'maximumRecordedGaugeResidual',maximumGauge, ...
        'productionMethods',{productionMethods},'scalars',scalar_values(state), ...
        'finalGaugeResidual',state.rhs.gaugeResidual,'rates',[state.rhs.c_l,state.rhs.c_omega,state.rhs.c_r]);
    row.analyticBackgroundVariableChange=isfield(state,'eta');
    row.peakDiagnostics=peak_summary(peakStageRecords);
    rows=append(rows,row);result.records=rows;
    report.temporal=result;report.elapsedSeconds=toc(timer);persist(folder,report);
    fprintf('FOURTH TIME n=%d dt=%.6g CFL=%.6g gauge=%.3e %.3fs\n',n,dt,maximumCfl,maximumGauge,row.seconds);
end
differences=struct([]);
for j=1:numel(states)-1
    a=states{j};b=states{j+1};fields=native_field_metrics(field_subtract(a.W,b.W),g);
    av=scalar_values(a);bv=scalar_values(b);delta=av.values-bv.values;
    delta(3)=(a.scale.physicalTime-b.scale.physicalTime)+ ...
        (a.scale.physicalTimeLow-b.scale.physicalTimeLow);
    scalarScale=max(1,abs(bv.values));
    differences=append(differences,struct('coarserDt',rows(j).dt,'finerDt',rows(j+1).dt, ...
        'W',fields,'scalarNames',{av.names},'scalarSignedDifferences',delta, ...
        'scalarAbsoluteDifferences',abs(delta),'scalarMixedRelativeDifferences',abs(delta)./scalarScale, ...
        'coupledScalarMixedInfinity',max(abs(delta)./scalarScale), ...
        'gaugeResidualDifferenceInfinity',max(abs(a.rhs.gaugeResidual-b.rhs.gaugeResidual)), ...
        'residualRepresentationPresent',isfield(a,'eta')&&isfield(b,'eta'), ...
        'eta',struct()));
    if differences(end).residualRepresentationPresent
        differences(end).eta=native_field_metrics(field_subtract(a.eta,b.eta),g);
    else,differences(end).eta=struct();end
end
result.differences=differences;orders=struct([]);
fieldErrors=[differences.W];fieldErrors=[fieldErrors.allActiveInfinity];
scalarErrors=[differences.coupledScalarMixedInfinity];
result.WOrders=log2(fieldErrors(1:end-1)./fieldErrors(2:end));
result.coupledScalarOrders=log2(scalarErrors(1:end-1)./scalarErrors(2:end));
Wscale=max([1;abs(states{end}.W.core(g.core.active));abs(states{end}.W.polar(g.polar.active))]);
result.WRoundoffScreen=1024*eps(Wscale);
result.fieldFourthOrderObserved=order_window(result.WOrders,fieldErrors,result.WRoundoffScreen);
result.coupledScalarFourthOrderObserved=order_window(result.coupledScalarOrders,scalarErrors,1024*eps);
for k=1:numel(differences(1).scalarNames)
    e=arrayfun(@(d)d.scalarAbsoluteDifferences(k),differences);
    p=log2(e(1:end-1)./e(2:end));
    scale=max(1,abs(result.records(end).scalars.values(k)));
    orders=append(orders,struct('name',differences(1).scalarNames{k},'absoluteDifferences',e, ...
        'observedOrders',p,'roundoffThreshold',1024*eps(scale), ...
        'fourthOrderObservedAboveRoundoff',order_window(p,e,1024*eps(scale))));
end
result.orders=orders;
result.allScalarFourthOrderObserved=all([orders.fourthOrderObservedAboveRoundoff]);
result.temporalInterpretation='Unextrapolated same-spatial-operator differences at nested exact time steps. A fourth-order window must be visible above the recorded rounding screen; a small or zero difference does not invent an order.';
end

function result=peak_summary(records)
% Stage observations diagnose candidate branch changes; no event is imposed.
result=struct('scope','All attempted stages, reused initial RHS and fresh endpoints', ...
    'records',records,'recordCount',numel(records),'availableCount',0, ...
    'uniqueCellIndices',[],'tieCount',0,'cellChangeCount',0,'xRange',[NaN,NaN], ...
    'stageChangesAreNotLocatedPhysicalEvents',true);
if isempty(records),return;end
valid=[records.available]&[records.rhsReturned];r=records(valid);
result.availableCount=numel(r);if isempty(r),return;end
cells=[r.cellIndex];cells=cells(isfinite(cells));
result.uniqueCellIndices=unique(cells);
result.cellChangeCount=sum(diff(cells)~=0);
result.tieCount=sum([r.tieNonunique]);
x=[r.x];x=x(isfinite(x));
if ~isempty(x),result.xRange=[min(x),max(x)];end
end

function yes=order_window(p,errors,floorValue)
yes=any(p>=3.5&p<=4.5&errors(1:end-1)>floorValue&errors(2:end)>floorValue);
end

function info=scalar_values(state)
s=state.scale;
info=struct('names',{{'logC_l','logC_omega','physicalTimeCompensated','X_shift','densityOrigin'}}, ...
    'values',[s.logC_l,s.logC_omega,s.physicalTime+s.physicalTimeLow,s.X_shift,s.densityOrigin], ...
    'physicalTimeHigh',s.physicalTime,'physicalTimeLow',s.physicalTimeLow);
end

function result=velocity_metrics(flow,exact,g)
field=struct();
for patch={'core','polar'}
    key=patch{1};field.(key)=hypot(flow.(key).u1-exact.(key).u1,flow.(key).u2-exact.(key).u2);
end
result=native_field_metrics(field,g);
end

function result=native_field_metrics(field,g)
result=struct('patches',struct(),'allActiveInfinity',0,'nearInfinity',0,'interfaceBandInfinity',0);
for patch={'core','polar'}
    key=patch{1};p=g.(key);v=field.(key);
    row=struct('allActive',metric(v,p.active,p.X,p.Y), ...
        'near',metric(v,p.active&hypot(p.X,p.Y)<=3,p.X,p.Y), ...
        'interfaceBand',metric(v,interface_band(g,key),p.X,p.Y));
    result.patches.(key)=row;result.allActiveInfinity=max(result.allActiveInfinity,row.allActive.absoluteInfinity);
    result.nearInfinity=max(result.nearInfinity,row.near.absoluteInfinity);
    result.interfaceBandInfinity=max(result.interfaceBandInfinity,row.interfaceBand.absoluteInfinity);
end
end

function mask=interface_band(g,key)
p=g.(key);mask=false(size(p.active));
if strcmp(key,'core')
    rows=p.activeRows;cols=p.activeColumns;
    mask(rows,cols(1:min(3,end)))=true;mask(rows,cols(max(1,end-2):end))=true;
    mask(rows(max(1,end-2):end),cols)=true;
else,mask(p.activeRows(1:3),p.activeColumns)=true;end
mask=mask&p.active;
end

function q=common_queries()
r=3*((0:48)'+.371)/49;theta=pi*((0:96)+.413)/97;
[R,T]=ndgrid(r,theta);x=R(:).*cos(T(:));y=R(:).*sin(T(:));
wall=3*((0:96)'+.617)/97;
q=struct('x',[x;-flipud(wall);0;wall],'y',[y;zeros(2*numel(wall)+1,1)],'endpoint',.6);
end

function g=make_geometry(level)
g=ipm_cr.mesh.build(struct('coreHalfWidth',1.2,'coreHeight',1.2, ...
    'coreSpacing',level.h,'polarStart',.72,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',level.nr,'ntheta',level.ntheta,'firstSpacing',level.h,'rmax',1e10)));
end

function plan=grid_plan()
plan=struct('h',{.06,.03,.015},'nr',{193,385,769},'ntheta',{65,129,257});
end

function info=spacing_summary(g)
p=g.polar;r=p.r(:);n=numel(r);targets=[1,3.5,4,20,100,r(end)];rows=struct([]);
for j=1:numel(targets)
    target=targets(j);lower=find(r<=target,1,'last');lower=max(1,min(n-1,lower));
    [~,node]=min(abs(r-target));first=max(1,min(n-6,node-3));last=first+6;
    metricWidth=p.Jr(node)/(n-1);
    rows=append(rows,struct('requestedRadius',target,'bracketLowerIndex',lower, ...
        'bracketLowerRadius',r(lower),'bracketUpperRadius',r(lower+1), ...
        'localDeltaR',r(lower+1)-r(lower),'nearestNativeIndex',node, ...
        'nearestNativeRadius',r(node),'JOverRadialIntervals',metricWidth, ...
        'sevenPointFirstIndex',first,'sevenPointLastIndex',last, ...
        'sevenPointPhysicalSpan',r(last)-r(first), ...
        'angularSpacing',p.dtheta,'physicalAngularSpacingAtRequestedRadius',target*p.dtheta));
end
info=struct('coreSpacing',[g.core.dx,g.core.dy], ...
    'originalRadialIntervals',p.originalNr-1,'retainedAnnularIntervals',n-1, ...
    'requestedFirstSpacing',g.gridOptions.firstSpacing, ...
    'radialMapFamily',g.gridOptions.radialMap, ...
    'angularIntervals',p.ntheta-1,'angularSpacing',p.dtheta, ...
    'probeRows',rows, ...
    'metricDefinition','Native polar Jr/(retainedAnnularNr-1), equal to the original mesh radial metric divided by originalNr-1.', ...
    'samplingDefinition','Gap brackets the fixed requested radius; metric and seven-node FD span use the nearest native radial node, with the actual one-sided closure at the outer endpoint.');
end

function value=field_combination(base,stages,weights,g)
value=base;
for k=1:numel(weights)
    if weights(k)==0,continue;end
    value.core(g.core.active)=value.core(g.core.active)+weights(k)*stages{k}.core(g.core.active);
    value.polar(g.polar.active)=value.polar(g.polar.active)+weights(k)*stages{k}.polar(g.polar.active);
end
end

function d=field_subtract(a,b)
d=struct('core',a.core-b.core,'polar',a.polar-b.polar);
end

function info=metric(v,mask,x,y)
indices=find(mask);assert(~isempty(indices),'ipm_cr:FourthRegion','A required metric region is empty.');
values=v(indices);assert(all(isfinite(values)),'ipm_cr:FourthFinite');
[maximum,k]=max(abs(values));index=indices(k);
info=struct('sampleCount',numel(indices),'absoluteInfinity',maximum,'rms',sqrt(mean(values.^2)), ...
    'maximumLocation',struct('X',x(index),'Y',y(index),'r',hypot(x(index),y(index))));
end

function inventory=source_inventory()
testPath=mfilename('fullpath');if ~endsWith(testPath,'.m'),testPath=[testPath,'.m'];end
projectRoot=fileparts(fileparts(testPath));files=dir(fullfile(projectRoot,'src','+ipm_cr','**','*.m'));
paths={testPath};
for j=1:numel(files),paths{end+1}=fullfile(files(j).folder,files(j).name);end %#ok<AGROW>
smoothSource=fullfile(projectRoot,'src','+ipm_cr','+op','wenoMex.cpp');
if isfile(smoothSource),paths{end+1}=smoothSource;end
smoothNative=which('ipm_cr.op.wenoMex');
if ~isempty(smoothNative),paths{end+1}=smoothNative;end
paths=sort(unique(paths));
inventory=struct([]);
for j=1:numel(paths)
    assert(~isempty(paths{j})&&isfile(paths{j}),'ipm_cr:FourthSource','A required numerical source is missing.');
    item=dir(paths{j});inventory=append(inventory,struct('path',paths{j},'bytes',item.bytes, ...
        'sha256',ipm_cr.io.sha(paths{j})));
end
end

function persist(folder,report)
temporary=fullfile(folder,'fourth_order.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'fourth_order.mat'),'f');
assert(ok,'ipm_cr:FourthSave','%s',message);
temporary=fullfile(folder,'fourth_order.json.writing');fid=fopen(temporary,'w');
assert(fid>=0,'ipm_cr:FourthSave');cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'fourth_order.json'),'f');
assert(ok,'ipm_cr:FourthSave','%s',message);
end

function a=append(a,b)
if isempty(a),a=b;else,a(end+1)=b;end
end
