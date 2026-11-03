function report=test_fem_green(outputDir,controls)
%TEST_FEM_GREEN Portable numerical-fixture regression of exact Green kernels.
% Independent image-Green quadrature and SciPy Li2 fixtures exercise all 32
% polynomial basis sources, including nonzero wall source and both corners.
% This finite quadrature regression is not a universal accuracy certificate
% or a convergence-order test of the FEM discretization. No time is advanced.
if nargin<2,controls=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','fem_green');end
assert(isstruct(controls)&&isscalar(controls)&&all(ismember(fieldnames(controls),{'wallSeconds'})), ...
    'ipm_cr:FEMGreenTestOptions','Only wallSeconds may be overridden.');
budget=30;if isfield(controls,'wallSeconds'),budget=controls.wallSeconds;end
validateattributes(budget,{'numeric'},{'real','finite','scalar','positive'});
outputDir=char(java.io.File(char(outputDir)).getCanonicalPath());
if ~isfolder(outputDir),mkdir(outputDir);end
path=fullfile(outputDir,'fem_green.json');
assert(~isfile(path),'ipm_cr:FEMGreenTestFresh','Choose a fresh result directory.');
timer=tic;inventory=ipm_cr.io.identity('source');
fixturePath=fullfile(root,'tests','reference','fem_green_reference.json');
fixture=jsondecode(fileread(fixturePath));
assert(strcmp(fixture.format,'fem_independent_green_reference_v1'), ...
    'ipm_cr:FEMGreenFixture','Unsupported independent reference fixture.');
report=struct('kind','portable_fem_analytic_green_regression_v1','status','running','passed',false, ...
    'timeAdvanced',0,'fixture',struct('path','tests/reference/fem_green_reference.json', ...
    'sha256',ipm_cr.io.sha(fixturePath),'provenance',fixture.provenance,'scope',fixture.scope), ...
    'testSHA256',ipm_cr.io.sha([mfilename('fullpath'),'.m']),'sourceInventory',inventory, ...
    'records',struct([]),'maximumNormalizedErrors',struct('psi',0,'u1',0,'u2',0), ...
    'endpoints',struct(),'dilog',struct(),'sourcesUnchanged',false, ...
    'continuumAccuracyCertified',false,'femFourthOrderCertified',false);
persist();
try
    data=fixture.green;powers=data.powers;x=reshape(data.x,1,[]);y=reshape(data.y,1,[]);
    assert(data.radius==1&&isequal(size(powers),[32,2])&&numel(x)==8&&numel(y)==8&& ...
        isequal(size(data.psi),[32,8])&&isequal(size(data.u1),[32,8])&&isequal(size(data.u2),[32,8]), ...
        'ipm_cr:FEMGreenFixture','Independent Green fixture dimensions changed.');
    tolerance=1e-13;
    for R=[1,3.7]
        for j=1:size(powers,1)
            assert(toc(timer)<budget,'ipm_cr:FEMGreenTestBudget');
            a=powers(j,1);b=powers(j,2);coeff=zeros(8,4);coeff(a+1,b+1)=1;
            actual=ipm_cr.poisson.fem.greenPolynomial(coeff,R*x,R*y,R);
            ref=struct('psi',R^2*data.psi(j,:),'u1',R*data.u1(j,:),'u2',R*data.u2(j,:));
            assert(all(isfinite([actual.psi(:);actual.u1(:);actual.u2(:);actual.q(:)])), ...
                'ipm_cr:FEMGreenFinite','Green response contains nonfinite values.');
            errors=struct('psi',0,'u1',0,'u2',0);
            for field={'psi','u1','u2'}
                key=field{1};scale=R;if strcmp(key,'psi'),scale=R^2;end
                errors.(key)=max(abs(actual.(key)-ref.(key)),[],'all')/scale;
                report.maximumNormalizedErrors.(key)=max(report.maximumNormalizedErrors.(key),errors.(key));
                assert(errors.(key)<=tolerance,'ipm_cr:FEMGreenReference', ...
                    'Independent Green fixture mismatch for x^%d*y^%d, R=%.17g, %s: %.17g.',a,b,R,key,errors.(key));
            end
            expectedQ=x.^a.*y.^b;
            assert(max(abs(actual.q-expectedQ),[],'all')<=tolerance,'ipm_cr:FEMGreenSource');
            row=struct('a',a,'b',b,'radius',R,'normalizedErrors',errors, ...
                'sourceAtPositiveWallCorner',actual.q(5),'sourceAtNegativeWallCorner',actual.q(6));
            if isempty(report.records),report.records=row;else,report.records(end+1)=row;end %#ok<AGROW>
        end
    end
    report.endpoints=endpoint_tests();
    report.dilog=dilog_tests(fixture.dilog,tolerance);
    report.sourcesUnchanged=isequaln(inventory,ipm_cr.io.identity('source'));
    assert(report.sourcesUnchanged,'ipm_cr:FEMGreenTestSources');
    assert(toc(timer)<budget,'ipm_cr:FEMGreenTestBudget');
    report.passed=true;report.status='passed_independent_green_and_dilog_fixtures';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    persist();rethrow(problem);
end
persist();

    function persist()
        report.elapsedSeconds=toc(timer);ipm_cr.io.json(path,report);
    end
end

function r=endpoint_tests()
% Exact source moments give these finite wall-corner velocities independently
% of the helper's modal/Laurent evaluation. q=1 and q=x remain nonzero there.
powers=[0,0;1,0;0,1];expected=[-1/pi,-1/pi;-1/(2*pi),1/(2*pi);-1/8,-1/8];
maximum=0;nearMaximum=0;wallMaximum=0;
for j=1:3
    c=zeros(8,4);c(powers(j,1)+1,powers(j,2)+1)=1;
    value=ipm_cr.poisson.fem.greenPolynomial(c,[1,-1],[0,0],1);
    maximum=max(maximum,max(abs(value.u1-expected(j,:))));
    wallMaximum=max(wallMaximum,max(abs([value.psi,value.u2])));
    close=ipm_cr.poisson.fem.greenPolynomial(c,[1-1e-12,-1+1e-12],[1e-13,1e-13],1);
    assert(all(isfinite([close.psi,close.u1,close.u2])),'ipm_cr:FEMGreenCornerFinite');
    nearMaximum=max(nearMaximum,max(abs(close.u1-expected(j,:))));
end
assert(maximum<=1e-14&&wallMaximum<=1e-14&&nearMaximum<=1e-9, ...
    'ipm_cr:FEMGreenCornerLimit','Finite nonzero-wall source corner limit failed.');
r=struct('powers',powers,'expectedCornerU1',expected,'maximumCornerU1Error',maximum, ...
    'maximumWallPsiOrU2',wallMaximum,'maximumNearCornerU1Difference',nearMaximum, ...
    'nonzeroWallSourcesIncluded',true,'nearCornerCoordinates',[1-1e-12,1e-13;-1+1e-12,1e-13]);
end

function r=dilog_tests(data,tolerance)
z=complex(data.zReal,data.zImag);reference=complex(data.valueReal,data.valueImag);
assert(numel(z)==876&&isequal(size(z),size(reference)), ...
    'ipm_cr:FEMGreenDilogFixture','SciPy Li2 fixture dimensions changed.');
actual=ipm_cr.poisson.fem.complexDilog(z);
error=max(abs(actual-reference),[],'all');
assert(all(isfinite(actual),'all')&&error<=tolerance,'ipm_cr:FEMGreenDilogReference', ...
    'Principal Li2 disagrees with independent SciPy fixture: %.17g.',error);
assert(max(abs(ipm_cr.poisson.fem.complexDilog([1,-1])-[pi^2/6,-pi^2/12]))<=1e-15, ...
    'ipm_cr:FEMGreenDilogEndpoints');
r=struct('sampleCount',numel(z),'maximumAbsoluteError',error,'referenceDefinition',data.definition, ...
    'closedUnitDisk',true,'positiveAndNegativeImaginaryAxesIncluded',true,'exactEndpointsChecked',true);
end
