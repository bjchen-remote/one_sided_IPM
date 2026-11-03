function report=test_fem_poisson(outputDir,controls)
%TEST_FEM_POISSON Bounded public FEM operator and complete-RHS regression.
% No evolution, private checkpoint, long job, or publication is performed.
% A smooth manufactured Dirichlet problem exercises nonzero wall source;
% it does not certify the finite-cutoff Green problem's corner accuracy.
if nargin<2,controls=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','fem_poisson');end
assert(isstruct(controls)&&isscalar(controls)&&all(ismember(fieldnames(controls),{'wallSeconds'})), ...
    'ipm_cr:FEMTestOptions','Only wallSeconds may be overridden.');
budget=120;if isfield(controls,'wallSeconds'),budget=controls.wallSeconds;end
validateattributes(budget,{'numeric'},{'real','finite','scalar','positive'});
outputDir=char(java.io.File(char(outputDir)).getCanonicalPath());
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'fem_poisson.json')),'ipm_cr:FEMTestFresh','Choose a fresh result directory.');
timer=tic;inventory=ipm_cr.io.identity('source');
report=struct('kind','portable_fem_operator_rhs_regression_v1','status','running', ...
    'passed',false,'timeAdvanced',0,'records',struct([]),'convergence',struct(), ...
    'options',struct(),'rhs',struct(),'sourceInventory',inventory, ...
    'sourcesUnchanged',false,'continuumAccuracyCertified',false, ...
    'nonzeroWallFiniteCutoffFourthOrderCertified',false);
persist();
try
    report.options=option_tests();
    radius=8;
    for L=[16,32]
        assert(toc(timer)<budget,'ipm_cr:FEMTestBudget');
        g=ipm_cr.mesh.build(geometry(L,radius));
        obj=ipm_cr.poisson.fem.solver(g,struct('outerBoundary','green_dirichlet', ...
            'nodalLayout','native','stride',[]));
        boundary=struct('kind','green_dirichlet','exteriorStatus','manufactured_boundary', ...
            'value',@(x,y)manufacture(x,y,radius,'psi'));
        source=@(x,y)manufacture(x,y,radius,'q');
        [flow,info]=obj.solve(source,boundary,[]);operator_gate(info);
        assert(info.sameNativeGrid&&info.noIndependentPoissonGrid&& ...
            info.coefficientDof==numel(g.core.X)+numel(g.polar.X)&& ...
            isequal(obj.c.x(:),g.core.x(:))&&isequal(obj.c.y(:),g.core.y(:))&& ...
            isequal(obj.p.r(:),g.polar.r(:))&&isequal(obj.p.theta(:),g.polar.theta(:)), ...
            'ipm_cr:FEMTestNativeGrid','All FEM unknown coordinates must be the exact transport nodes.');
        assert(isequal(flow.core.psi,flow.coefficients.core)&& ...
            isequal(flow.polar.psi,flow.coefficients.polar),'ipm_cr:FEMTestNativeReturn');
        row=struct('L',L,'coreShape',size(g.core.X),'polarShape',size(g.polar.X), ...
            'sourceAtWallOrigin',source(0,0),'errors',errors(flow,g,radius), ...
            'information',info,'splitDifferenceInfinity',NaN,'ghostPoisonDifferenceInfinity',NaN);
        % A constant native residual is reproduced exactly. The remaining
        % continuous source is evaluated at Gauss points, never resampled
        % as a complete native source and then counted a second time.
        native=struct('core',.25*ones(size(g.core.X)),'polar',.25*ones(size(g.polar.X)));
        offset=@(x,y)source(x,y)-.25;
        [split,splitInfo]=obj.solve(native,boundary,offset);operator_gate(splitInfo);
        row.splitDifferenceInfinity=distance(flow,split);
        native.core(g.core.ghost)=1e12;native.polar(g.polar.ghost)=-1e12;
        [poisoned,poisonedInfo]=obj.solve(native,boundary,offset);operator_gate(poisonedInfo);
        row.ghostPoisonDifferenceInfinity=distance(split,poisoned);
        if isempty(report.records),report.records=row;else,report.records(end+1)=row;end
        persist();
        scale=max(1,max(abs([flow.core.psi(:);flow.polar.psi(:); ...
            flow.core.u1(:);flow.core.u2(:);flow.polar.u1(:);flow.polar.u2(:)])));
        assert(source(0,0)==-2,'ipm_cr:FEMTestWallSource','Nonzero wall source was changed.');
        assert(row.splitDifferenceInfinity<2e-10*scale&&row.ghostPoisonDifferenceInfinity==0, ...
            'ipm_cr:FEMTestSource','The analytic source split or active-only ownership changed.');
        assert(~info.strongResidualAvailable&&~info.outerSourceTailCertified&&~info.massLumped, ...
            'ipm_cr:FEMTestScope');
        assert(norm(obj.c.xAxis.M-diag(diag(obj.c.xAxis.M)),inf)>0, ...
            'ipm_cr:FEMTestMass','The assembled Galerkin mass must not be lumped.');
        clear obj
    end
    coarse=report.records(1).errors;fine=report.records(2).errors;
    report.convergence=struct('problem','smooth Dirichlet manufacture with q(0,0)=-2', ...
        'levels',[16,32],'activePsiOrder',log(coarse.activePsi/fine.activePsi)/log(2), ...
        'activeUOrder',log(coarse.activeU/fine.activeU)/log(2), ...
        'scope','bounded regression; does not qualify a finite-cutoff Green corner');
    persist();
    assert(fine.activePsi<coarse.activePsi/8&&fine.activeU<coarse.activeU/4, ...
        'ipm_cr:FEMTestRefinement','Smooth manufactured errors did not decrease under refinement.');
    assert(toc(timer)<budget,'ipm_cr:FEMTestBudget');
    report.rhs=rhs_tests();
    report.sourcesUnchanged=isequaln(inventory,ipm_cr.io.identity('source'));
    assert(report.sourcesUnchanged,'ipm_cr:FEMTestSources');
    assert(toc(timer)<budget,'ipm_cr:FEMTestBudget');
    report.passed=true;report.status='passed_bounded_operator_and_complete_rhs_checks';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    persist();rethrow(problem);
end
persist();

    function persist()
        report.elapsedSeconds=toc(timer);
        ipm_cr.io.json(fullfile(outputDir,'fem_poisson.json'),report);
    end
end

function r=option_tests()
defaults=ipm_cr.config();
assert(strcmp(defaults.poissonBackend,'fem')&& ...
    strcmp(defaults.femOptions.linearSolver,'sparse_direct')&& ...
    defaults.femOptions.degree==5&&isempty(defaults.femOptions.stride)&& ...
    strcmp(defaults.femOptions.nodalLayout,'native')&& ...
    strcmp(defaults.transport.reconstruction,'linear5'),'ipm_cr:FEMTestDefault');
resolved=ipm_cr.config(struct('poissonBackend','fem','femOptions',struct('quadOrder',10)));
assert(strcmp(resolved.poissonBackend,'fem')&&resolved.femOptions.degree==5&& ...
    resolved.femOptions.quadOrder==10&&isempty(resolved.femOptions.stride)&& ...
    strcmp(resolved.femOptions.nodalLayout,'native')&& ...
    isequaln(resolved,ipm_cr.config(resolved)),'ipm_cr:FEMTestOptions');
reject(@()ipm_cr.config(struct('femOptions',struct('quadorder',8))),'ipm_cr:FEMOptions');
reject(@()ipm_cr.config(struct('femOptions',struct('exteriorModel','unknown_tail'))),'ipm_cr:FEMExterior');
reject(@()ipm_cr.config(struct('femOptions',struct('outerBoundary','green_dirichlet'))),'ipm_cr:FEMExterior');
reject(@()ipm_cr.config(struct('femOptions',struct('stride',2))),'ipm_cr:FEMNativeGrid');
reject(@()ipm_cr.config(struct('femOptions',struct('nodalLayout','gll_reference'))),'ipm_cr:FEMNativeGrid');
compatible=ipm_cr.config(struct('femOptions',struct('stride',5)));
assert(strcmp(compatible.femOptions.nodalLayout,'native'),'ipm_cr:FEMTestNativeCompatibility');
r=struct('defaultIsSparseFem',true,'defaultTransportIsLinear5',true, ...
    'femPartialOptionsRoundtrip',true, ...
    'unknownOptionRejected',true,'undeclaredExteriorRejected',true, ...
    'unsupportedRuntimeBoundaryRejected',true,'independentStrideAndGllRejected',true, ...
    'degreeCompatibleStrideAccepted',true);
end

function r=rhs_tests()
overrides=struct('poissonBackend','fem','geometry',geometry(16,1e10), ...
    'solverOptions',struct('rOuter',1e10));
[s,ctx]=ipm_cr.init(overrides);g=ctx.mesh;
assert(s.step==0&&s.scale.canonicalTime==0&&isfield(s,'eta')&& ...
    s.rhs.analyticBackgroundSource,'ipm_cr:FEMTestInitial');
assert(strcmp(ctx.opts.transport.reconstruction,'linear5')&& ...
    strcmp(ctx.opts.transport.wenoAlphaMode,'local')&& ...
    ctx.opts.transport.smoothSpeedCoefficient==.1,'ipm_cr:FEMTestTransport');
p=s.rhs.poisson;
operator_gate(p);
assert(strcmp(p.poissonBackend,'fem')&& ...
    ~p.strongResidualAvailable&&isnan(p.strongResidual)&&~p.outerSourceTailCertified&& ...
    strcmp(p.exteriorModel,'finite_source_disk_zero_exterior')&&p.sameNativeGrid&& ...
    p.noIndependentPoissonGrid&&p.coefficientDof==numel(g.core.X)+numel(g.polar.X), ...
    'ipm_cr:FEMTestMetadata');
assert(max(abs(s.rhs.gaugeResidual))<=1e-5,'ipm_cr:FEMTestGauge');
assert(p.wallNormalVelocityInfinity<1e-10,'ipm_cr:FEMTestWall');
zero=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
[full,fi]=ctx.velocity(zero,[],true);[light,li]=ctx.velocity(zero,[],false);
operator_gate(fi);operator_gate(li);
assert(distance(full,light)==0&&isempty(light.core.psi)&&isempty(light.polar.psi)&& ...
    fi.solverCompleted&&li.solverCompleted,'ipm_cr:FEMTestZeroExterior');
assert(max(abs([full.core.psi(:);full.polar.psi(:);full.core.u1(:);full.core.u2(:); ...
    full.polar.u1(:);full.polar.u2(:)]))==0,'ipm_cr:FEMTestZero');
direct=ipm_cr.poisson.fem.solver(g,struct('outerBoundary','forced_dtn','nodalLayout','native','stride',[]));
[~,boundary]=ipm_cr.poisson.fem.options();
missing=rmfield(boundary,'exteriorStatus');unknown=boundary;unknown.exteriorStatus='unknown_tail';
reject(@()direct.solve(zero),'ipm_cr:FEMExteriorData');
reject(@()direct.solve(zero,missing,[]),'ipm_cr:FEMExteriorData');
reject(@()direct.solve(zero,unknown,[]),'ipm_cr:FEMExteriorData');
r=struct('tau',s.scale.canonicalTime,'step',s.step,'gaugeResidualInfinity',max(abs(s.rhs.gaugeResidual)), ...
    'wallNormalVelocityInfinity',p.wallNormalVelocityInfinity,'poisson',p,'stats',ctx.stats(), ...
    'analyticBackgroundSplitRetained',true,'zeroExteriorFieldExact',true, ...
    'missingAndUnknownExteriorRejected',true,'lightVelocityUnchanged',true);
end

function operator_gate(info)
% Accuracy qualification belongs to this independent test, not solve().
assert(info.solverCompleted&&~info.accuracyAssessed&&info.finiteOutput&& ...
    ~isfield(info,'accepted'),'ipm_cr:FEMTestCompletion');
assert(info.weakBackwardError<1e-10&&info.traceAbsolute<1e-9&& ...
    info.essentialBoundaryAbsolute<1e-9,'ipm_cr:FEMTestOperator', ...
    'Original weak-row, trace and essential-boundary accuracy gates failed.');
end

function geo=geometry(L,radius)
a=1.2;geo=struct('coreHalfWidth',a,'coreHeight',a,'coreSpacing',a/L, ...
    'polarStart',.65,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',4*L+1,'ntheta',2*L+1,'rmax',radius,'firstSpacing',a/(2*L), ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01,'stencilWidth',7), ...
    'coreMap',struct('enabled',true,'intervals',[2*L,L],'minSpacing',[a/L,a/L]/sinh(1)));
end

function value=manufacture(x,y,R,field)
switch field
    case 'psi',value=y.^2.*(1-(x.^2+y.^2)/R^2);
    case 'q',value=-2+(2*x.^2+14*y.^2)/R^2;
    case 'u1',value=-2*y+(2*x.^2.*y+4*y.^3)/R^2;
    case 'u2',value=-2*x.*y.^2/R^2;
end
end

function r=errors(flow,g,R)
r=struct('activePsi',0,'activeU',0,'fullPsi',0,'fullU',0);
for patch={'core','polar'}
    key=patch{1};p=g.(key);f=flow.(key);
    ep=abs(f.psi-manufacture(p.X,p.Y,R,'psi'));
    eu=max(abs(f.u1-manufacture(p.X,p.Y,R,'u1')),abs(f.u2-manufacture(p.X,p.Y,R,'u2')));
    r.activePsi=max(r.activePsi,max(ep(p.active)));r.activeU=max(r.activeU,max(eu(p.active)));
    r.fullPsi=max(r.fullPsi,max(ep,[],'all'));r.fullU=max(r.fullU,max(eu,[],'all'));
end
end

function d=distance(a,b)
d=0;
for patch={'core','polar'}
    key=patch{1};
    for field={'u1','u2'}
        name=field{1};d=max(d,max(abs(a.(key).(name)-b.(key).(name)),[],'all'));
    end
    if ~isempty(a.(key).psi)&&~isempty(b.(key).psi)
        d=max(d,max(abs(a.(key).psi-b.(key).psi),[],'all'));
    end
end
end

function reject(callback,identifier)
try
    callback();error('ipm_cr:FEMTestExpectedFailure','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:FEMTestWrongFailure', ...
        'Expected %s, received %s.',identifier,problem.identifier);
end
end
