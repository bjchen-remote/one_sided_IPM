function report=test_mainline_accuracy(outputDir,controls)
%TEST_MAINLINE_ACCURACY Actual-mesh MMS and complete-RK4 time convergence.
% Completion, finite-pair order gates and complexity qualification are
% separate. Default runs both studies; no production default is changed.
if nargin<2,controls=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','mainline_accuracy');end
assert(isstruct(controls)&&isscalar(controls)&&all(ismember(fieldnames(controls), ...
    {'wallSeconds','runSpatial','runTemporal'})), ...
    'ipm_cr:MainlineAccuracyOptions','Use wallSeconds, runSpatial and runTemporal.');
budget=600;if isfield(controls,'wallSeconds'),budget=controls.wallSeconds;end
validateattributes(budget,{'numeric'},{'scalar','real','finite','positive','<=',600});
runSpatial=true;runTemporal=true;
if isfield(controls,'runSpatial'),runSpatial=controls.runSpatial;end
if isfield(controls,'runTemporal'),runTemporal=controls.runTemporal;end
validateattributes(runSpatial,{'logical'},{'scalar'});validateattributes(runTemporal,{'logical'},{'scalar'});
assert(runSpatial||runTemporal,'ipm_cr:MainlineAccuracyOptions','Select at least one study.');
if ~isfolder(outputDir),mkdir(outputDir);end
path=fullfile(outputDir,'mainline_accuracy.json');
assert(~isfile(path),'ipm_cr:MainlineAccuracyFresh','Choose a fresh result directory.');
oldThreads=maxNumCompThreads(1);threadCleanup=onCleanup(@()maxNumCompThreads(oldThreads)); %#ok<NASGU>
timer=tic;inventory=ipm_cr.io.identity('source');cloud=queries();
report=struct('kind','actual_mesh_spatial_and_complete_RK4_accuracy_v1', ...
    'status','running','passed',false,'timeAdvanced',0,'sourceInventory',inventory, ...
    'requestedStudies',struct('spatial',runSpatial,'temporal',runTemporal), ...
    'sourceDefinition','psi=y^2*exp(-(x^2+y^2)/1.5^2); q=-Delta psi; u=(-psi_y,psi_x)', ...
    'exteriorModel','known complete smooth source; exact Dirichlet trace, no finite-disk corner enrichment', ...
    'setup',struct('degrees',[4,5],'quadOrder',8,'nodalLayout','native','linearSolver','sparse_direct', ...
        'warmupSolves',10,'timedSolves',20,'threads',1,'wallSeconds',budget), ...
    'registeredOrderGates',struct('metric','fixedCommonU, maximum component error at identical physical queries', ...
        'pair',{{'coarse','medium'}},'degree',[4,5],'minimumOrder',[3.5,4.2], ...
        'roundoffResolutionFloor',5e-12,'interpretation', ...
        'finite-pair screening only; unresolved floor-level errors do not pass or fail an order claim'), ...
    'queryRegions',cloud.description,'records',struct([]),'skippedDiagnostics',{{}}, ...
    'finitePair',struct([]),'warmScaling',struct([]),'temporal',struct(),'sourcesUnchanged',false, ...
    'qualification',struct('finitePairOrderGatesPassed',false,'asymptoticOrderCertified',false, ...
        'productionSpatialScreenPassed',false, ...
        'completeRK4TimeScreenPassed',false, ...
        'nearLinearComplexityCertified',false,'velocityClassicallyDifferentiableCertified',false, ...
        'rssMeasured',false,'rssInterpretation','peak RSS must be supplied by the external bounded launcher'));
persist();
try
    if runSpatial
    for level={'quarter','coarse','medium'}
        name=level{1};check_budget(timer,budget);g=[];
        try
            if strcmp(name,'quarter'),g=quarter_geometry();
            else,g=ipm_cr.mesh.build(struct('preset',name));end
        catch problem
            if strcmp(name,'quarter')&&geometry_rejection(problem)
                skip(name,0,problem);continue
            end
            rethrow(problem)
        end
        geometry=geometry_record(g,name);
        native=struct('core',mms(g.core.X,g.core.Y,'q'), ...
            'polar',mms(g.polar.X,g.polar.Y,'q'));
        boundary=struct('kind','green_dirichlet','exteriorStatus','known_complete_source', ...
            'value',@(x,y)mms(x,y,'psi'));
        for degree=[4,5]
            check_budget(timer,budget);obj=[];phase=tic;
            try
                obj=ipm_cr.poisson.fem.solver(g,struct('degree',degree,'quadOrder',8, ...
                    'nodalLayout','native','stride',[],'outerBoundary','green_dirichlet', ...
                    'linearSolver','sparse_direct'));
            catch problem
                if strcmp(name,'quarter')&&geometry_rejection(problem)
                    skip(name,degree,problem);continue
                end
                rethrow(problem)
            end
            setupSeconds=toc(phase);
            for input={'analytic','native'}
                check_budget(timer,budget);mode=input{1};source=native;
                if strcmp(mode,'analytic'),source=@(x,y)mms(x,y,'q');end
                phase=tic;[flow,info]=obj.solve(source,boundary);firstSolveSeconds=toc(phase);
                layout=native_contract(obj,flow,g,info);
                assert(info.solverCompleted&&~info.accuracyAssessed&&info.finiteOutput&& ...
                    ~isfield(info,'accepted')&&info.weakBackwardError<1e-10&&info.traceAbsolute<1e-9&& ...
                    info.essentialBoundaryAbsolute<1e-9,'ipm_cr:MainlineAccuracyOperator', ...
                    'Original weak/trace/boundary gates failed for %s degree%d %s.',name,degree,mode);
                errors=measure(obj,flow,g,cloud);
                timing=warm_timing(obj,source,boundary,timer,budget);
                row=struct('level',name,'degree',degree,'nominalVelocityOrder',degree, ...
                    'input',mode,'geometry',geometry,'setupSeconds',setupSeconds, ...
                    'firstSolveSeconds',firstSolveSeconds,'errors',errors,'warm',timing, ...
                    'information',info,'nativeLayout',layout,'greenResponseCacheBytes',0,'rssMeasured',false, ...
                    'derivativeMeaning','native D7 applied to each returned FE velocity array, not a broken-element Hessian or a classical derivative on faces');
                if isempty(report.records),report.records=row;else,report.records(end+1)=row;end %#ok<AGROW>
                clear flow;persist();
            end
            % Only one FE object (including its source maps and factors)
            % is retained at a time; do not accumulate medium contexts.
            clear obj
        end
        clear g native
    end
    [report.finitePair,report.warmScaling]=comparisons(report.records,report.registeredOrderGates);
    report.qualification.finitePairOrderGatesPassed= ...
        numel(report.finitePair)==4&&all([report.finitePair.gatePassed]);
    defaultMethod=ipm_cr.config();report.productionDegree=defaultMethod.femOptions.degree;
    selected=report.finitePair([report.finitePair.degree]==report.productionDegree);
    report.qualification.productionSpatialScreenPassed=numel(selected)==2&&all([selected.gatePassed]);
    end
    if runTemporal
        report.temporal=temporal_test(outputDir,timer,budget,@publish_temporal);
        report.timeAdvanced=report.temporal.targetTau;
        report.qualification.completeRK4TimeScreenPassed=report.temporal.timeScreenPassed;
    end
    report.sourcesUnchanged=isequaln(inventory,ipm_cr.io.identity('source'));
    assert(report.sourcesUnchanged,'ipm_cr:MainlineAccuracySources');check_budget(timer,budget);
    report.passed=true;report.status='completed_requested_measurements_qualification_recorded_separately';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');persist();rethrow(problem);
end
persist();

    function persist()
        report.elapsedSeconds=toc(timer);ipm_cr.io.json(path,report);
    end
    function skip(name,degree,problem)
        item=struct('level',name,'degree',degree,'status','diagnostic_geometry_not_admissible', ...
            'errorIdentifier',problem.identifier,'reason',problem.message, ...
            'meshChangedToBypassGuard',false,'acceptanceGrid',false);
        report.skippedDiagnostics{end+1}=item;persist();
    end
    function publish_temporal(value)
        report.temporal=value;persist();
    end
end

function g=quarter_geometry()
% This is an extra diagnostic, not the user's coarse acceptance mesh.
c=ipm_cr.mesh.build(struct('preset','coarse'));
a=struct('x',c.core.x(1:2:end),'y',c.core.y(1:2:end), ...
    'r',c.polar.r(1:2:end),'theta',c.polar.theta(1:2:end), ...
    'dxMetric',2*c.core.dxMetric(1:2:end),'dyMetric',2*c.core.dyMetric(1:2:end), ...
    'drMetric',2*c.polar.drMetric(1:2:end));
g=ipm_cr.mesh.build(struct('axes',a));
end
function yes=geometry_rejection(problem)
yes=any(strcmp(problem.identifier,{'ipm_cr:ExplicitOverlap','ipm_cr:FEMOverlap', ...
    'ipm_cr:ExplicitGeometry','ipm_cr:ExplicitAxes'}));
end
function record=native_contract(obj,flow,g,info)
assert(strcmp(info.nodalLayout,'native')&&info.sameNativeGrid&&info.noIndependentPoissonGrid&& ...
    info.coefficientDof==numel(g.core.X)+numel(g.polar.X)&&g.ghostWidth==3&& ...
    isequal(obj.c.x(:),g.core.x(:))&&isequal(obj.c.y(:),g.core.y(:))&& ...
    isequal(obj.p.r(:),g.polar.r(:))&&isequal(obj.p.theta(:),g.polar.theta(:)), ...
    'ipm_cr:MainlineAccuracyNativeGrid','Accuracy/cost comparisons must use exact native unknowns.');
assert(isequal(flow.core.psi,flow.coefficients.core)&&isequal(flow.polar.psi,flow.coefficients.polar), ...
    'ipm_cr:MainlineAccuracyNativeReturn','The raw native streamfunction must not be interpolated back.');
d=[obj.c.xAxis.elementDegrees(:);obj.c.yAxis.elementDegrees(:); ...
    obj.p.tAxis.elementDegrees(:);obj.p.thetaAxis.elementDegrees(:)];
assert(all(d==info.degree|d==info.degree+1),'ipm_cr:MainlineAccuracyNativeDegree');
record=struct('axesExactlyNative',true,'coefficientDof',info.coefficientDof, ...
    'targetDegree',info.degree,'actualDegreeRange',[min(d),max(d)], ...
    'higherDegreeRemainderElements',nnz(d==info.degree+1),'streamfunctionDirectNativeReturn',true);
end
function r=geometry_record(g,name)
r=struct('name',name,'coreShape',size(g.core.X),'polarShape',size(g.polar.X), ...
    'nativeFullDof',numel(g.core.X)+numel(g.polar.X),'nativeActiveDof',g.activeDofCount, ...
    'outerRadius',g.polar.r(end),'coreMinimumMetricSpacing',[min(g.core.dxMetric),min(g.core.dyMetric)], ...
    'activeCoreXBounds',g.core.activeXBounds,'activeCoreYBounds',g.core.activeYBounds, ...
    'innerActiveRadius',g.polar.innerActiveRadius,'ghostWidth',g.ghostWidth);
end
function c=queries()
[x,y]=meshgrid(linspace(-2.5,2.5,43)+.0031,linspace(.013,1.8,31));
c.core=struct('x',x,'y',y);
[theta,r]=meshgrid(linspace(.013,pi-.013,47),logspace(log10(2.6),log10(12),19));
c.polar=struct('x',r.*cos(theta),'y',r.*sin(theta));
x=linspace(-4.5,4.5,121);c.coreWall=struct('x',x,'y',zeros(size(x)));
r=logspace(log10(2.6),10,97);c.polarWall=struct('x',[r,-r],'y',zeros(1,2*numel(r)));
theta=linspace(0,pi,99);c.outer=struct('x',1e10*cos(theta),'y',1e10*sin(theta));
c.outer.y([1,end])=0;
c.description=struct('coreOffgrid',[min(c.core.x,[],'all'),max(c.core.x,[],'all'),.013,1.8], ...
    'polarOffgridRadius',[2.6,12],'polarOffgridTheta',[.013,pi-.013], ...
    'wallCoreX',[-4.5,4.5],'wallPolarRadius',[2.6,1e10], ...
    'outerRadius',1e10,'samePhysicalQueriesAtEveryLevel',true);
end
function e=measure(obj,flow,g,c)
e=struct('nativeFullPsi',0,'nativeActivePsi',0,'nativeFullU',0,'nativeActiveU',0, ...
    'nativeFullJacobian',0,'nativeActiveJacobian',0,'nativeWallPsi',0,'nativeWallU',0, ...
    'fixedCommonPsi',0,'fixedCommonU',0,'offgridCorePsi',0,'offgridCoreU',0, ...
    'offgridPolarPsi',0,'offgridPolarU',0,'fixedWallPsi',0,'fixedWallU',0, ...
    'fixedOuterPsi',0,'fixedOuterU',0);
d1=ipm_cr.mesh.diff(struct('core',flow.core.u1,'polar',flow.polar.u1),g);
d2=ipm_cr.mesh.diff(struct('core',flow.core.u2,'polar',flow.polar.u2),g);
for patch={'core','polar'}
    key=patch{1};p=g.(key);v=flow.(key);
    ep=abs(v.psi-mms(p.X,p.Y,'psi'));
    eu=max(abs(v.u1-mms(p.X,p.Y,'u1')),abs(v.u2-mms(p.X,p.Y,'u2')));
    ej=max(max(abs(d1.(key).x-mms(p.X,p.Y,'u1x')),abs(d1.(key).y-mms(p.X,p.Y,'u1y'))), ...
        max(abs(d2.(key).x-mms(p.X,p.Y,'u2x')),abs(d2.(key).y-mms(p.X,p.Y,'u2y'))));
    e.nativeFullPsi=max(e.nativeFullPsi,max(ep,[],'all'));
    e.nativeActivePsi=max(e.nativeActivePsi,max(ep(p.active)));
    e.nativeFullU=max(e.nativeFullU,max(eu,[],'all'));
    e.nativeActiveU=max(e.nativeActiveU,max(eu(p.active)));
    e.nativeFullJacobian=max(e.nativeFullJacobian,max(ej,[],'all'));
    e.nativeActiveJacobian=max(e.nativeActiveJacobian,max(ej(p.active)));
    wall=p.Y==0;e.nativeWallPsi=max(e.nativeWallPsi,max(ep(wall)));
    e.nativeWallU=max(e.nativeWallU,max(eu(wall)));
end
[e.offgridCorePsi,e.offgridCoreU]=query_error(obj,flow,c.core,'core');
[e.offgridPolarPsi,e.offgridPolarU]=query_error(obj,flow,c.polar,'polar');
e.fixedCommonPsi=max(e.offgridCorePsi,e.offgridPolarPsi);
e.fixedCommonU=max(e.offgridCoreU,e.offgridPolarU);
[cp,cu]=query_error(obj,flow,c.coreWall,'core');
[pp,pu]=query_error(obj,flow,c.polarWall,'polar');
e.fixedWallPsi=max(cp,pp);e.fixedWallU=max(cu,pu);
[e.fixedOuterPsi,e.fixedOuterU]=query_error(obj,flow,c.outer,'polar');
% Compare the actual two element limits in FE coordinates. A physical
% atan2/log roundtrip can miss an exact face by one ULP and hide a jump.
cc=flow.coefficients.core;pa=flow.coefficients.polar;
e.coreXFaceVelocityJump=face_jump(obj.c.xAxis,cc,false,[]);
e.coreYFaceVelocityJump=face_jump(obj.c.yAxis,cc,true,[]);
e.polarRadialFaceVelocityJump=face_jump(obj.p.tAxis,pa,true,'radial');
e.polarAngularFaceVelocityJump=face_jump(obj.p.thetaAxis,pa,false,obj.p.r(:));
assert(all(isfinite(cell2mat(struct2cell(e)))),'ipm_cr:MainlineAccuracyFinite');
end
function jump=face_jump(a,coefficients,rows,scale)
faces=a.breaks(2:end-1);jump=0;if isempty(faces),return;end
[~,left]=ipm_cr.poisson.fem.assembly('evaluate',a,faces,'left');
[~,right]=ipm_cr.poisson.fem.assembly('evaluate',a,faces,'right');
if rows,value=(left-right)*coefficients;else,value=coefficients*(left-right)';end
if ischar(scale),value=value./exp(faces(:));
elseif ~isempty(scale),value=value./scale;end
jump=max(abs(value),[],'all');
end
function [ep,eu]=query_error(obj,flow,q,chart)
v=obj.evaluate(flow,q.x,q.y,chart);
ep=max(abs(v.psi-mms(q.x,q.y,'psi')),[],'all');
eu=max([abs(v.u1(:)-mms(q.x(:),q.y(:),'u1'));abs(v.u2(:)-mms(q.x(:),q.y(:),'u2'))]);
end
function r=warm_timing(obj,source,boundary,timer,budget)
for j=1:10,check_budget(timer,budget);obj.solve(source,boundary);end
samples=zeros(1,20);
for j=1:20
    check_budget(timer,budget);phase=tic;obj.solve(source,boundary);samples(j)=toc(phase);
end
middle=median(samples);ordered=sort(samples);
r=struct('samplesSeconds',samples,'medianSeconds',middle, ...
    'medianAbsoluteDeviationSeconds',median(abs(samples-middle)), ...
    'interquartileWidthSeconds',ordered(15)-ordered(5), ...
    'minimumSeconds',min(samples),'maximumSeconds',max(samples), ...
    'includes','source/load, cached sparse solves, full native psi/u return, original row and overlap diagnostics', ...
    'excludes','constructor, first native-source map build, offgrid MMS error audit, JSON output');
end
function [pair,scaling]=comparisons(records,gates)
pair=struct([]);scaling=struct([]);
for degree=[4,5]
    for input={'analytic','native'}
        mode=input{1};select=[records.degree]==degree&strcmp({records.input},mode);
        selected=records(select);a=selected(strcmp({selected.level},'coarse'));
        b=selected(strcmp({selected.level},'medium'));
        assert(numel(a)==1&&numel(b)==1,'ipm_cr:MainlineAccuracyRequiredMeshes');
        names=fieldnames(a.errors);orders=struct();
        for k=1:numel(names)
            key=names{k};orders.(key)=slope(a.errors.(key),b.errors.(key));
        end
        resolved=a.errors.fixedCommonU>gates.roundoffResolutionFloor&& ...
            b.errors.fixedCommonU>gates.roundoffResolutionFloor;
        minimum=gates.minimumOrder(gates.degree==degree);
        row=struct('degree',degree,'input',mode,'levels',{{'coarse','medium'}}, ...
            'orders',orders,'minimumFixedCommonUOrder',minimum,'resolvedAboveFloor',resolved, ...
            'gatePassed',resolved&&orders.fixedCommonU>=minimum,'asymptoticOrderCertified',false);
        if isempty(pair),pair=row;else,pair(end+1)=row;end %#ok<AGROW>
        na=a.information.coefficientDof;nb=b.information.coefficientDof;
        fa=a.information.linearSolverStorage;fb=b.information.linearSolverStorage;
        row=struct('degree',degree,'input',mode,'coefficientDof',[na,nb], ...
            'warmMedianSeconds',[a.warm.medianSeconds,b.warm.medianSeconds], ...
            'warmEmpiricalExponent',exponent(a.warm.medianSeconds,b.warm.medianSeconds,na,nb), ...
            'patchFactorNnz',[fa.patchFactorNnz,fb.patchFactorNnz], ...
            'patchFactorNnzEmpiricalExponent',exponent(fa.patchFactorNnz,fb.patchFactorNnz,na,nb), ...
            'patchFactorBytes',[fa.patchFactorBytes,fb.patchFactorBytes], ...
            'assembledOperatorNnz',[a.information.sparseOperatorNnz,b.information.sparseOperatorNnz], ...
            'setupSeconds',[a.setupSeconds,b.setupSeconds],'theoreticalComplexityCertified',false, ...
            'scope','two-scale empirical ratio; factor bytes exclude source/return maps, Schur and assembled operators');
        if isempty(scaling),scaling=row;else,scaling(end+1)=row;end %#ok<AGROW>
    end
end
end
function p=slope(a,b)
p=NaN;if a>0&&b>0,p=log(a/b)/log(2);end
end
function p=exponent(a,b,na,nb)
p=NaN;if a>0&&b>0&&nb>na,p=log(b/a)/log(nb/na);end
end
function value=mms(x,y,quantity)
b=1/1.5^2;e=exp(-b*(x.^2+y.^2));
switch quantity
    case 'psi',value=y.^2.*e;
    case 'q',value=(-2+12*b*y.^2-4*b^2*y.^2.*(x.^2+y.^2)).*e;
    case 'u1',value=-2*y.*(1-b*y.^2).*e;
    case 'u2',value=-2*b*x.*y.^2.*e;
    case 'u1x',value=4*b*x.*y.*(1-b*y.^2).*e;
    case 'u1y',value=(-2+10*b*y.^2-4*b^2*y.^4).*e;
    case 'u2x',value=(-2*b*y.^2+4*b^2*x.^2.*y.^2).*e;
    case 'u2y',value=-4*b*x.*y.*(1-b*y.^2).*e;
end
end
function result=temporal_test(folder,timer,budget,publish)
targetTau=.4;nominalDt=[.01,.005,.0025];rows=struct([]);paths=cell(1,3);
defaults=ipm_cr.config();degree=defaults.femOptions.degree;transport=defaults.transport;
result=struct('kind',['complete_default_FEM_',transport.reconstruction,'_RK4_fixed_coarse_self_convergence_v2'], ...
    'targetTau',targetTau,'nominalMaxDt',nominalDt,'preset','coarse','degree',degree, ...
    'transport',transport,'transportSelectedFromCurrentDefaults',true, ...
    'remeshEnabled',false,'spaceOperatorHeldFixed',true,'originalIPCRetained',true, ...
    'directTableauStep',true,'stageCflMeasured',true,'cflSelectsInitialInterval',true, ...
    'timePartition','ordinary original-spline-knot crossings; fixed four secant iterates of knot slope', ...
    'sspGuaranteeClaimed',false, ...
    'registeredScreen',struct('minimumOrder',3.6,'relativeResolutionThreshold',1e-12, ...
        'requiredField','W','primaryScalars',{{'logC_l','logC_omega','X_shift','densityOrigin','physicalTimeCompensated'}}, ...
        'unresolvedPolicy','no observed order is assigned below threshold; unresolved scalars are excluded from the resolved-scalar gate'), ...
    'records',rows,'differences',struct([]),'fieldOrders',struct([]),'scalarOrders',struct([]), ...
    'sameInitialStateAndRhsIdentity',false,'WTimeScreenPassed',false, ...
    'resolvedPrimaryScalarScreenPassed',false,'timeScreenPassed',false, ...
    'globalSpatialOrderCertified',false,'longTrajectoryQualified',false);
publish(result);referenceKey=[];initialReference=[];
for level=1:3
    check_budget(timer,budget);runTimer=tic;state=[];ctx=[];
    method=struct('geometry',struct('preset','coarse'),'poissonBackend','fem', ...
        'femOptions',struct('degree',degree,'quadOrder',8,'nodalLayout','native','stride',[], ...
            'linearSolver','sparse_direct'), ...
        'transport',transport, ...
        'remesh',struct('enabled',false),'timeIntegrator','rk4','cfl',.3,'maxDt',nominalDt(level));
    [state,ctx]=ipm_cr.init(method);g=ctx.mesh;
    assert(strcmp(ctx.opts.poissonBackend,'fem')&&isequaln(ctx.opts.transport,transport)&& ...
        ctx.opts.femOptions.degree==degree&&~ctx.opts.remesh.enabled&& ...
        state.rhs.poisson.sameNativeGrid&&state.rhs.poisson.noIndependentPoissonGrid&& ...
        state.rhs.poisson.coefficientDof==numel(g.core.X)+numel(g.polar.X), ...
        'ipm_cr:MainlineTimeMethod');
    temporal_poisson_gate(state.rhs.poisson);
    initial=struct('W',state.W,'eta',[],'scale',state.scale,'originalAuthority',state.originalAuthority);
    if isfield(state,'eta'),initial.eta=state.eta;end
    if level==1,referenceKey=ctx.key;initialReference=initial;
    else
        assert(isequaln(ctx.key,referenceKey)&&isequaln(initial,initialReference), ...
            'ipm_cr:MainlineTimeIdentity','Time refinement changed the initial data, IPC or spatial method.');
    end
    steps=struct([]);maximumCfl=0;maximumGauge=0;rhsCalls=0;partitionCalls=0;partitions=struct([]);
    n=round(targetTau/nominalDt(level));
    try
        for interval=1:n
            goal=targetTau*interval/n;
            while state.scale.canonicalTime<goal
                check_budget(timer,budget);before=state.scale.canonicalTime;
                expectedDt=min([ctx.opts.maxDt,.95*ctx.opts.cfl/state.rhs.rate,goal-before]);
                [state,info]=ipm_cr.step(state,ctx,goal);
                assert(state.scale.canonicalTime>before&&strcmp(info.timeIntegrator,'rk4')&& ...
                    info.initialProposedDt==expectedDt&&info.dt==info.timePartition.segmentDt&& ...
                    info.dt>0&&info.dt<=expectedDt&&info.cflStepFormulaApplied&&info.terminalRhsEvaluated&& ...
                    isequal(info.stageCfl,info.dt*info.stageRates)&& ...
                    info.realizedMaximumStageCfl<=.3*(1+128*eps), ...
                    'ipm_cr:MainlineTimeStep','Direct RK4 interval, endpoint or independently measured stage CFL failed.');
                temporal_poisson_gate(state.rhs.poisson);
                assert(info.physicalClockUpdate.wordCount==4&& ...
                    isequal(clock_words(state.scale),info.physicalClockUpdate.afterWords), ...
                    'ipm_cr:MainlineTimeClockWords');
                maximumCfl=max(maximumCfl,info.realizedMaximumStageCfl);
                maximumGauge=max(maximumGauge,info.maximumStageGaugeResidual);
                rhsCalls=rhsCalls+info.rhsEvaluations;
                partitionCalls=partitionCalls+info.timePartition.extraRhsEvaluations;
                if info.timePartition.fitted
                    fitted=info.timePartition;fitted.endpointTau=state.scale.canonicalTime;
                    if isempty(partitions),partitions=fitted;else,partitions(end+1)=fitted;end %#ok<AGROW>
                end
                entry=struct('step',state.step,'nominalInterval',interval,'tauBefore',before, ...
                    'tauAfter',state.scale.canonicalTime,'dt',info.dt, ...
                    'maximumStageCfl',info.realizedMaximumStageCfl,'stageCfl',info.stageCfl, ...
                    'stageRates',info.stageRates,'stageAbscissae',info.stageAbscissae, ...
                    'stageGaugeResidualInfinity',info.stageGaugeResidualInfinity, ...
                    'rhsEvaluations',info.rhsEvaluations,'stagePeakDiagnostics',info.stagePeakDiagnostics, ...
                    'timePartition',info.timePartition, ...
                    'scale',state.scale,'physicalClockWords',clock_words(state.scale));
                if isempty(steps),steps=entry;else,steps(end+1)=entry;end %#ok<AGROW>
            end
            assert(state.scale.canonicalTime==goal,'ipm_cr:MainlineTimeEndpoint');
        end
    catch problem
        % Preserve the actual last completed state if this independent test
        % stops before the requested endpoint.
        failedPath=fullfile(folder,sprintf('temporal_level%d_last_completed.mat',level));
        save_time_state(failedPath,state,ctx,steps);
        result.failure=struct('level',level,'completedTau',state.scale.canonicalTime, ...
            'checkpoint',failedPath,'sha256',ipm_cr.io.sha(failedPath), ...
            'errorIdentifier',problem.identifier,'errorMessage',problem.message);
        publish(result);rethrow(problem)
    end
    assert(state.scale.canonicalTime==targetTau&&maximumGauge<=1e-5, ...
        'ipm_cr:MainlineTimeFinalState');
    paths{level}=fullfile(folder,sprintf('temporal_level%d_final.mat',level));
    save_time_state(paths{level},state,ctx,steps);
    row=struct('level',level,'nominalDt',nominalDt(level),'nominalIntervals',n, ...
        'completedSteps',numel(steps),'actualTau',state.scale.canonicalTime, ...
        'actualDt',[steps.dt],'minimumActualDt',min([steps.dt]),'maximumActualDt',max([steps.dt]), ...
        'maximumStageCfl',maximumCfl,'maximumStageGaugeResidual',maximumGauge, ...
        'rhsEvaluations',rhsCalls,'partitionExtraRhsEvaluations',partitionCalls, ...
        'fittedKnotCount',numel(partitions),'fittedKnots',partitions, ...
        'scalars',time_scalars(state),'scale',state.scale,'physicalClockWords',clock_words(state.scale), ...
        'finalGaugeResidual',state.rhs.gaugeResidual,'elapsedSeconds',toc(runTimer), ...
        'finalStateFile',paths{level},'finalStateSHA256',ipm_cr.io.sha(paths{level}), ...
        'checkpointContainsNativeStateAndGeometry',true,'coldReplayVerified',false, ...
        'sourceAndInputIdentity',ctx.identity,'stepRecords',steps);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    result.records=rows;result.sameInitialStateAndRhsIdentity=true;publish(result);
    clear ctx state g initial
end
differences=struct([]);
for pair=1:2
    a=load(paths{pair},'state','geometry');b=load(paths{pair+1},'state');
    row=struct('coarserDt',nominalDt(pair),'finerDt',nominalDt(pair+1), ...
        'fields',struct(),'scalarNames',{{}},'scalarSignedDifferences',[],'scalarAbsoluteDifferences',[]);
    for name={'W','q','F','Q'}
        key=name{1};av=time_field(a.state,key);bv=time_field(b.state,key);
        row.fields.(key)=field_difference(av,bv,a.geometry);
    end
    av=time_scalars(a.state);bv=time_scalars(b.state);delta=av.values-bv.values;
    delta(5)=compensated_sum([clock_words(a.state.scale),-clock_words(b.state.scale)]);
    row.scalarNames=av.names;row.scalarSignedDifferences=delta;row.scalarAbsoluteDifferences=abs(delta);
    if isempty(differences),differences=row;else,differences(end+1)=row;end %#ok<AGROW>
    clear a b
end
result.differences=differences;finest=load(paths{3},'state','geometry');
orders=struct([]);
for name={'W','q','F','Q'}
    key=name{1};e=[differences(1).fields.(key).allActiveInfinity,differences(2).fields.(key).allActiveInfinity];
    scale=field_magnitude(time_field(finest.state,key),finest.geometry);
    entry=time_order(key,e,scale,strcmp(key,'W'));
    if isempty(orders),orders=entry;else,orders(end+1)=entry;end %#ok<AGROW>
end
result.fieldOrders=orders;result.WTimeScreenPassed=orders(1).gatePassed;
scalars=time_scalars(finest.state);orders=struct([]);
for k=1:numel(scalars.names)
    e=[differences(1).scalarAbsoluteDifferences(k),differences(2).scalarAbsoluteDifferences(k)];
    entry=time_order(scalars.names{k},e,max(1,abs(scalars.values(k))),k<=5);
    if isempty(orders),orders=entry;else,orders(end+1)=entry;end %#ok<AGROW>
end
result.scalarOrders=orders;primary=[orders.required]&[orders.resolvedAboveThreshold];
result.resolvedPrimaryScalarScreenPassed=any(primary)&&all([orders(primary).gatePassed]);
result.timeScreenPassed=result.WTimeScreenPassed&&result.resolvedPrimaryScalarScreenPassed;
result.interpretation='Actual unextrapolated full coupled solutions with identical fixed spatial discretization; below-threshold differences remain unresolved and receive no invented order.';
publish(result);
end
function temporal_poisson_gate(info)
assert(info.solverCompleted&&~info.accuracyAssessed&&info.finiteOutput&& ...
    ~isfield(info,'accepted'),'ipm_cr:MainlineTimePoissonCompletion');
assert(info.weakBackwardError<1e-10&&info.traceAbsolute<1e-9&& ...
    info.essentialBoundaryAbsolute<1e-9,'ipm_cr:MainlineTimePoissonAccuracy');
end
function save_time_state(path,state,ctx,stepRecords)
geometry=ctx.mesh;method=ctx.opts;wallOps=ctx.wall;sourceAndInputIdentity=ctx.identity; %#ok<NASGU>
temporary=[path,'.writing.mat'];
save(temporary,'state','geometry','method','wallOps','sourceAndInputIdentity','stepRecords','-v7.3');
[ok,message]=movefile(temporary,path,'f');assert(ok,'ipm_cr:MainlineAccuracySave','%s',message);
end
function f=time_field(state,name)
if strcmp(name,'W'),f=state.W;else,f=state.rhs.(name);end
assert(isstruct(f)&&all(isfield(f,{'core','polar'}))&& ...
    all(isfinite(f.core),'all')&&all(isfinite(f.polar),'all'),'ipm_cr:MainlineTimeFullField');
end
function r=field_difference(a,b,g)
r=struct('allFullInfinity',0,'allActiveInfinity',0,'nearActiveInfinity',0);
for name={'core','polar'}
    key=name{1};p=g.(key);d=abs(a.(key)-b.(key));near=p.active&hypot(p.X,p.Y)<=3;
    r.allFullInfinity=max(r.allFullInfinity,max(d,[],'all'));
    r.allActiveInfinity=max(r.allActiveInfinity,max(d(p.active)));
    if any(near,'all'),r.nearActiveInfinity=max(r.nearActiveInfinity,max(d(near)));end
end
end
function scale=field_magnitude(f,g)
scale=max([1;abs(f.core(g.core.active));abs(f.polar(g.polar.active))]);
end
function row=time_order(name,errors,scale,required)
threshold=1e-12*scale;resolved=errors(2)>threshold;p=NaN;
if resolved,p=slope(errors(1),errors(2));end
row=struct('quantity',name,'absoluteDifferences',errors,'scale',scale, ...
    'resolutionThreshold',threshold,'resolvedAboveThreshold',resolved,'observedOrder',p, ...
    'minimumOrder',3.6,'required',required,'gateApplicable',resolved, ...
    'gatePassed',resolved&&isfinite(p)&&p>=3.6);
end
function s=time_scalars(state)
p=state.scale;words=clock_words(p);
s=struct('names',{{'logC_l','logC_omega','X_shift','densityOrigin','physicalTimeCompensated', ...
    'c_l','c_omega','c_r','delta','Uref'}}, ...
    'values',[p.logC_l,p.logC_omega,p.X_shift,p.densityOrigin,compensated_sum(words), ...
        state.rhs.c_l,state.rhs.c_omega,state.rhs.c_r,state.rhs.delta,state.rhs.Uref], ...
    'physicalClockWords',words,'clockDifferenceUsesAllWords',true);
end
function words=clock_words(scale)
tail=[0,0];if isfield(scale,'physicalTimeTail'),tail=scale.physicalTimeTail;end
words=[scale.physicalTime,scale.physicalTimeLow,reshape(tail,1,2)];
end
function value=compensated_sum(words)
% Grow an error-free expansion before taking the final binary64 difference;
% never combine the two large clocks before subtracting their lower words.
e=[];
for x=words
    q=x;next=[];
    for y=e
        s=q+y;b=s-q;r=(q-(s-b))+(y-b);q=s;
        if r~=0,next(end+1)=r;end %#ok<AGROW>
    end
    if q~=0||isempty(next),next(end+1)=q;end %#ok<AGROW>
    e=next;
end
value=sum(e);
end
function check_budget(timer,budget)
assert(toc(timer)<budget,'ipm_cr:MainlineAccuracyBudget','Static accuracy budget exhausted.');
end
