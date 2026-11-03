function report=test_fem_native_grid(outputDir)
%TEST_FEM_NATIVE_GRID Exact transport/Poisson nodes and real Galerkin forms.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'result','fem_native_grid');end
if ~isfolder(outputDir),mkdir(outputDir);end
report=struct('passed',false,'records',struct([]),'lastSolveDiagnostic',struct(), ...
    'sources',ipm_cr.io.identity('source'));
for name={'coarse','medium'}
    g=ipm_cr.mesh.build(struct('preset',name{1}));
    obj=ipm_cr.poisson.fem.solver(g,struct('degree',5,'quadOrder',8, ...
        'nodalLayout','native','outerBoundary','green_dirichlet'));
    assert(isequal(obj.c.x(:),g.core.x(:))&&isequal(obj.c.y(:),g.core.y(:))&& ...
        isequal(obj.p.r(:),g.polar.r(:))&&isequal(obj.p.theta(:),g.polar.theta(:))&& ...
        isequal(obj.c.X,g.core.X)&&isequal(obj.c.Y,g.core.Y)&& ...
        isequal(obj.p.X,g.polar.X)&&isequal(obj.p.Y,g.polar.Y), ...
        'ipm_cr:FEMNativeTestCoordinates','Same counts alone do not prove the same grid.');
    axes={obj.c.xAxis,obj.c.yAxis,obj.p.tAxis,obj.p.thetaAxis};basisError=zeros(1,4);
    for j=1:4
        a=axes{j};span=a.nodes(end)-a.nodes(1);z=(a.nodes-a.nodes(1))/span;
        [E,D]=ipm_cr.poisson.fem.assembly('evaluate',a,a.qx);
        truth=(a.qx-a.nodes(1))/span;
        for p=0:5
            v=z.^p;expected=zeros(size(truth));
            if p>0,expected=p*truth.^(p-1)/span;end
            basisError(j)=max([basisError(j);abs(E*v-truth.^p);abs(D*v-expected)]);
        end
        assert(basisError(j)<1e-9&&a.quadratureCount>a.dofCount&&~a.massLumped, ...
            'ipm_cr:FEMNativeTestBasis','Native nodal polynomials or independent mass quadrature failed.');
        assert(a.symmetryRoundingBeforeProjection<1e-12,'ipm_cr:FEMNativeTestSymmetry');
        [~,flag]=chol(a.M);assert(flag==0,'ipm_cr:FEMNativeTestMass');
    end
    psi=@(x,y)y.^2.*exp(-(x.^2+y.^2)/2.25);
    q=@(x,y)(-2+12*y.^2/2.25-4*y.^2.*(x.^2+y.^2)/2.25^2).*exp(-(x.^2+y.^2)/2.25);
    boundary=struct('kind','green_dirichlet','exteriorStatus','known_complete_source','value',psi);
    [f,info]=obj.solve(q,boundary);
    report.lastSolveDiagnostic=struct('preset',name{1},'input','analytic','information',info);
    ipm_cr.io.json(fullfile(outputDir,'fem_native_grid.json'),report);
    operator_gate(info,[name{1},' analytic']);
    assert(info.solverCompleted&&info.coefficientDof==numel(g.core.X)+numel(g.polar.X)&& ...
        isequal(f.core.psi,f.coefficients.core)&&isequal(f.polar.psi,f.coefficients.polar)&& ...
        isequal(size(f.core.psi),size(g.core.X))&&isequal(size(f.polar.psi),size(g.polar.X)), ...
        'ipm_cr:FEMNativeTestUnknowns','Streamfunction must be the native coefficient array, without regridding.');
    assert(max(abs(f.core.psi),[],'all')>.1&& ...
        max(abs(f.core.u2(1,:)))<1e-12&&max(abs(f.polar.u2(:,[1,end])),[],'all')<1e-12);
    native=struct('core',q(g.core.X,g.core.Y),'polar',q(g.polar.X,g.polar.Y));
    [n,ni]=obj.solve(native,boundary);
    report.lastSolveDiagnostic=struct('preset',name{1},'input','native','information',ni);
    ipm_cr.io.json(fullfile(outputDir,'fem_native_grid.json'),report);
    operator_gate(ni,[name{1},' native']);defect=0;
    for chart={'core','polar'}
        k=chart{1};for field={'psi','u1','u2'}
            v=field{1};defect=max(defect,max(abs(n.(k).(v)-f.(k).(v)),[],'all'));
        end
    end
    assert(ni.solverCompleted&&defect<1e-5,'ipm_cr:FEMNativeTestQuadratureSource');
    row=struct('preset',name{1},'exactFourAxes',true,'nativePsiExactlyCoefficients',true, ...
        'nativeDof',numel(g.core.X)+numel(g.polar.X),'basisPolynomialError',basisError, ...
        'nativeVersusAnalyticSourceDefect',defect,'information',info);
    if isempty(report.records),report.records=row;else,report.records(end+1)=row;end %#ok<AGROW>
    clear obj f n g
end
expect(@()ipm_cr.config(struct('femOptions',struct('nodalLayout','gll_reference'))), ...
    'ipm_cr:FEMNativeGrid');
expect(@()ipm_cr.config(struct('femOptions',struct('stride',2))), 'ipm_cr:FEMNativeGrid');
report.sourcesUnchanged=isequaln(report.sources,ipm_cr.io.identity('source'));
assert(report.sourcesUnchanged);report.passed=true;
test_report(outputDir,'fem_native_grid',report);
end
function operator_gate(info,label)
assert(info.solverCompleted&&~info.accuracyAssessed&&info.finiteOutput&& ...
    ~isfield(info,'accepted'),'ipm_cr:FEMNativeTestCompletion', ...
    '%s: completed=%d, accuracyAssessed=%d, finiteOutput=%d.', ...
    label,info.solverCompleted,info.accuracyAssessed,info.finiteOutput);
assert(info.weakBackwardError<1e-10&&info.traceAbsolute<1e-9&& ...
    info.essentialBoundaryAbsolute<1e-9,'ipm_cr:FEMNativeTestWeakRows', ...
    '%s: weak=%.17g (core=%.17g, polar=%.17g; limit 1e-10), trace=%.17g, essential=%.17g (limits 1e-9).', ...
    label,info.weakBackwardError,info.coreWeakBackwardError,info.polarWeakBackwardError, ...
    info.traceAbsolute,info.essentialBoundaryAbsolute);
end
function expect(call,identifier)
try,call();catch problem,assert(strcmp(problem.identifier,identifier),'Unexpected %s',problem.identifier);return;end
error('ipm_cr:FEMNativeExpectedRejection','Expected %s.',identifier);
end
