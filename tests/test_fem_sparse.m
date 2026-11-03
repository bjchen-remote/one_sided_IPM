function report=test_fem_sparse(outputDir)
%TEST_FEM_SPARSE Same weak operator, sparse warm solves versus modal oracle.
% This does not advance time or certify continuum order/complexity. Both
% exterior operators, native/analytic splitting and Green enrichment are
% compared, including coefficients, full velocity and original row gates.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','fem_sparse');end
if ~isfolder(outputDir),mkdir(outputDir);end
path=fullfile(outputDir,'fem_sparse.json');
assert(~isfile(path),'ipm_cr:FEMSparseTestFresh','Choose a fresh result directory.');
timer=tic;inventory=ipm_cr.io.identity('source');
report=struct('kind','same_operator_sparse_modal_fem_regression_v1', ...
    'passed',false,'status','running','timeAdvanced',0,'records',struct([]), ...
    'enriched',struct(),'sourcesUnchanged',false,'sourceInventory',inventory, ...
    'continuumAccuracyCertified',false,'nearLinearComplexityCertified',false);
persist();
try
    defaults=ipm_cr.poisson.fem.options();
    assert(strcmp(defaults.linearSolver,'sparse_direct')&&defaults.degree==5&& ...
        strcmp(defaults.nodalLayout,'native')&&isempty(defaults.stride),'ipm_cr:FEMSparseTestDefault');
    for L=[16,32]
        g=ipm_cr.mesh.build(geometry(L));
        for kind={'green_dirichlet','forced_dtn'}
            outer=kind{1};options=struct('outerBoundary',outer,'linearSolver','sparse_direct', ...
                'nodalLayout','native','stride',[]);
            a=ipm_cr.poisson.fem.solver(g,options);options.linearSolver='modal';
            b=ipm_cr.poisson.fem.solver(g,options);
            boundary=struct('kind',outer,'exteriorStatus','manufactured_boundary', ...
                'value',@(x,y)y.^2.*(1-(x.^2+y.^2)/64));
            if strcmp(outer,'forced_dtn')
                boundary.exteriorStatus='declared_compact_source';boundary.value=@(x,y)zeros(size(x));
            end
            source=@(x,y)-2+(2*x.^2+14*y.^2)/64;
            [fa,ia]=a.solve(source,boundary);[fb,ib]=b.solve(source,boundary);
            operator_gate(ia);operator_gate(ib);
            layout=native_contract(a,fa,g,ia);native_contract(b,fb,g,ib);
            difference=distance(fa,fb);scale=magnitude(fb);
            assert(difference<2e-10*scale,'ipm_cr:FEMSparseModalMismatch', ...
                'L%d %s sparse/modal full-field difference %.17g, scale %.17g.',L,outer,difference,scale);
            assert(~isfield(a.coupling,'coreResponse')&&~isfield(a.c,'Vx')&& ...
                ~isfield(a.p,'harmonic')&&issparse(a.c.factor.L)&&issparse(a.p.factor.L)&& ...
                ia.linearSolverStorage.coreResponseBytes==0&&~ia.denseModalWarmTransforms, ...
                'ipm_cr:FEMSparseStorage','Sparse mode retained dense warm responses/transforms.');
            native=struct('core',.25*ones(size(g.core.X)),'polar',.25*ones(size(g.polar.X)));
            offset=@(x,y)source(x,y)-.25;
            [split,splitInfo]=a.solve(native,boundary,offset);operator_gate(splitInfo);
            splitDifference=distance(fa,split);
            assert(splitDifference<2e-10*scale,'ipm_cr:FEMSparseSplit');
            [repeated,repeatedInfo]=a.solve(native,boundary,offset);operator_gate(repeatedInfo);
            assert(distance(split,repeated)==0,'ipm_cr:FEMSparseRepeat');
            row=struct('L',L,'outerBoundary',outer,'fullFieldAndCoefficientDifference',difference, ...
                'referenceMagnitude',scale,'nativeAnalyticSplitDifference',splitDifference, ...
                'coefficientDof',ia.coefficientDof,'sparseOperatorNnz',ia.sparseOperatorNnz, ...
                'nativeLayout',layout, ...
                'storage',ia.linearSolverStorage,'sparseInformation',ia,'modalInformation',ib);
            if isempty(report.records),report.records=row;else,report.records(end+1)=row;end %#ok<AGROW>
            persist();clear a b
        end
    end
    report.enriched=enrichment_test();
    report.smoothExterior=smooth_exterior_test();
    report.largeDiskConstantSplit=large_disk_split();
    report.sourcesUnchanged=isequaln(inventory,ipm_cr.io.identity('source'));
    assert(report.sourcesUnchanged,'ipm_cr:FEMSparseSources');
    report.passed=true;report.status='passed_same_operator_sparse_modal_equivalence';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');persist();rethrow(problem);
end
persist();

    function persist()
        report.elapsedSeconds=toc(timer);ipm_cr.io.json(path,report);
    end
end

function r=large_disk_split()
geo=geometry(16);geo.grid.rmax=1e10;g=ipm_cr.mesh.build(geo);
a=ipm_cr.poisson.fem.solver(g,struct('outerBoundary','forced_dtn','nodalLayout','native','stride',[]));
boundary=struct('kind','forced_dtn','exteriorStatus','declared_compact_source', ...
    'value',@(x,y)zeros(size(x)));
source=@(x,y)exp(-x.^2-y.^2);
native=struct('core',.125*ones(size(g.core.X)),'polar',.125*ones(size(g.polar.X)));
[value,map]=ipm_cr.mesh.sample(native,g,a.p.qX,a.p.qY);
assert(all(value==.125,'all')&&isequal(value,ipm_cr.mesh.sample(native,map)), ...
    'ipm_cr:FEMSourceConstant','Cold and cached sampling must preserve constants exactly.');
[reference,referenceInfo]=a.solve(source,boundary);operator_gate(referenceInfo);
native_contract(a,reference,g,referenceInfo);
[split,splitInfo]=a.solve(native,boundary,@(x,y)source(x,y)-.125);operator_gate(splitInfo);
[repeat,repeatInfo]=a.solve(native,boundary,@(x,y)source(x,y)-.125);operator_gate(repeatInfo);
defect=distance(reference,split);scale=magnitude(reference);
assert(defect<1e-10*scale&&distance(split,repeat)==0, ...
    'ipm_cr:FEMSourceSplitRange','r^2 amplified a mathematically cancelling source split: %.6g.',defect);
r=struct('radius',1e10,'constant',.125,'constantSamplingExactlyPreserved',true, ...
    'fullFieldAndCoefficientDefect',defect,'referenceMagnitude',scale,'cachedRepeatExact',true);
end

function r=enrichment_test()
g=ipm_cr.mesh.build(geometry(16));
options=struct('outerBoundary','forced_dtn','linearSolver','sparse_direct','nodalLayout','native','stride',[]);
a=ipm_cr.poisson.fem.enriched(g,options);options.linearSolver='modal';
b=ipm_cr.poisson.fem.enriched(g,options);
boundary=struct('kind','forced_dtn','exteriorStatus','declared_compact_source', ...
    'value',@(x,y)zeros(size(x)));
source=@(x,y)1+.1*(1-(x.^2+y.^2)/64).^4;
[fa,ia]=a.solve(source,boundary);[fb,ib]=b.solve(source,boundary);
operator_gate(ia);operator_gate(ib);
layout=native_contract(a,fa,g,ia);native_contract(b,fb,g,ib);
difference=distance(fa,fb);scale=magnitude(fb);
assert(difference<2e-10*scale&&ia.cornerCorrection&&ib.cornerCorrection&& ...
    isequal(fa.enrichment.coefficients,fb.enrichment.coefficients)&& ...
    ia.solverCompleted&&ib.solverCompleted&&source(8,0)==1,'ipm_cr:FEMSparseGreen');
r=struct('nonzeroOuterWallSource',1,'fullFieldAndCoefficientDifference',difference, ...
    'referenceMagnitude',scale,'nativeLayout',layout,'sparseInformation',ia,'modalInformation',ib);
end

function r=smooth_exterior_test()
% A smooth, noncompact source has no finite-cutoff corner singularity.
% The optional enrichment wrapper must preserve the raw Galerkin problem
% exactly for all noncompact/unspecified-cutoff exterior statuses.
g=ipm_cr.mesh.build(geometry(16));options=struct('outerBoundary','green_dirichlet','nodalLayout','native','stride',[]);
a=ipm_cr.poisson.fem.enriched(g,options);b=ipm_cr.poisson.fem.solver(g,options);
source=@(x,y)smooth_field(x,y,'q');r=struct([]);
for status={'known_complete_source','modeled_tail_with_bound','manufactured_boundary'}
    boundary=struct('kind','green_dirichlet','exteriorStatus',status{1}, ...
        'value',@(x,y)smooth_field(x,y,'psi'));
    if strcmp(status{1},'modeled_tail_with_bound'),boundary.errorBound=0;end
    [fa,ia]=a.solve(source,boundary);[fb,ib]=b.solve(source,boundary);
    operator_gate(ia);operator_gate(ib);
    native_contract(a,fa,g,ia);native_contract(b,fb,g,ib);
    difference=distance(fa,fb);
    assert(difference==0&&ia.solverCompleted&&ib.solverCompleted&&~ia.cornerCorrection&& ...
        ia.cornerCorrectionSkipped&&~isfield(fa,'enrichment')&&a.Counters.cornerFits==0&& ...
        source(8,0)~=0,'ipm_cr:FEMKnownExteriorEnrichment', ...
        'A known smooth exterior problem must not acquire finite-disk corner enrichment.');
    x=[-.77,.13,.81];y=[.21,.38,.67];
    va=a.evaluate(fa,x,y,'core');vb=b.evaluate(fb,x,y,'core');
    [ac,ap]=a.overlap(fa);[bc,bp]=b.overlap(fb);
    assert(isequal(va,vb)&&isequal(ac,bc)&&isequal(ap,bp), ...
        'ipm_cr:FEMKnownExteriorEvaluation');
    row=struct('exteriorStatus',status{1},'fullFieldAndCoefficientDifference',difference, ...
        'outerWallSource',source(8,0),'cornerFits',a.Counters.cornerFits, ...
        'reason',ia.cornerCorrectionReason,'evaluateAndOverlapIdentical',true);
    if isempty(r),r=row;else,r(end+1)=row;end %#ok<AGROW>
end
end
function value=smooth_field(x,y,quantity)
e=exp(-(x.^2+y.^2)/25);f=y.^2+.3*y;
if strcmp(quantity,'psi'),value=f.*e;
else,value=-e.*(2-4*y.*(2*y+.3)/25+(4*(x.^2+y.^2)/625-4/25).*f);end
end

function g=geometry(L)
a=1.2;g=struct('coreHalfWidth',a,'coreHeight',a,'coreSpacing',a/L, ...
    'polarStart',.65,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',4*L+1,'ntheta',2*L+1,'rmax',8,'firstSpacing',a/(2*L), ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01,'stencilWidth',7), ...
    'coreMap',struct('enabled',true,'intervals',[2*L,L],'minSpacing',[a/L,a/L]/sinh(1)));
end
function value=magnitude(f)
value=1;
for name={'core','polar'}
    k=name{1};value=max([value;abs(f.(k).psi(:));abs(f.(k).u1(:));abs(f.(k).u2(:));abs(f.coefficients.(k)(:))]);
end
end
function value=distance(a,b)
value=0;
for name={'core','polar'}
    k=name{1};
    for field={'psi','u1','u2'}
        q=field{1};value=max(value,max(abs(a.(k).(q)-b.(k).(q)),[],'all'));
    end
    value=max(value,max(abs(a.coefficients.(k)-b.coefficients.(k)),[],'all'));
end
end

function operator_gate(info)
assert(info.solverCompleted&&~info.accuracyAssessed&&info.finiteOutput&& ...
    ~isfield(info,'accepted'),'ipm_cr:FEMSparseCompletion');
assert(info.weakBackwardError<1e-10&&info.traceAbsolute<1e-9&& ...
    info.essentialBoundaryAbsolute<1e-9,'ipm_cr:FEMSparseOriginalRows', ...
    'Independent original weak/trace/boundary checks failed for the direct scheme.');
end

function record=native_contract(obj,flow,g,info)
assert(strcmp(info.nodalLayout,'native')&&info.sameNativeGrid&&info.noIndependentPoissonGrid&& ...
    info.coefficientDof==numel(g.core.X)+numel(g.polar.X)&&g.ghostWidth==3, ...
    'ipm_cr:FEMSparseNativeDof','Native FEM must retain every physical node, including ghost collars.');
assert(isequal(obj.c.x(:),g.core.x(:))&&isequal(obj.c.y(:),g.core.y(:))&& ...
    isequal(obj.p.r(:),g.polar.r(:))&&isequal(obj.p.theta(:),g.polar.theta(:)), ...
    'ipm_cr:FEMSparseNativeAxes','Coordinate equality must be exact, not merely equal shapes.');
for patch={'core','polar'}
    key=patch{1};assert(isequal(size(flow.coefficients.(key)),size(g.(key).X)));
    if ~isfield(flow,'enrichment')
        assert(isequal(flow.(key).psi,flow.coefficients.(key)),'ipm_cr:FEMSparseNativeReturn');
    end
end
d=[obj.c.xAxis.elementDegrees(:);obj.c.yAxis.elementDegrees(:); ...
    obj.p.tAxis.elementDegrees(:);obj.p.thetaAxis.elementDegrees(:)];
assert(all(d==info.degree|d==info.degree+1),'ipm_cr:FEMSparseNativeDegree');
axes={obj.c.xAxis,obj.c.yAxis,obj.p.tAxis,obj.p.thetaAxis};
for j=1:numel(axes)
    assert(axes{j}.symmetryRoundingBeforeProjection<1e-12,'ipm_cr:FEMSparseAssemblySymmetry');
end
assert(rcond(obj.coupling.S)>1e-12,'ipm_cr:FEMSparseCouplingCondition');
if isfield(obj.c,'Vx')
    assert(all(obj.c.lx>0)&&all(obj.c.ly>0)&&obj.c.eigenDefectX<1e-9&& ...
        obj.c.eigenDefectY<1e-9,'ipm_cr:FEMSparseCoreEigensystem');
end
if isfield(obj.p,'V')
    assert(all(obj.p.mu>0)&&obj.p.eigenDefect<1e-9,'ipm_cr:FEMSparsePolarEigensystem');
end
record=struct('axesExactlyNative',true,'coefficientDof',info.coefficientDof, ...
    'targetDegree',info.degree,'actualDegreeRange',[min(d),max(d)], ...
    'higherDegreeRemainderElements',nnz(d==info.degree+1),'fullGhostLayersRetained',true);
end
