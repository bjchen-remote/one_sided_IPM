function report=test_ipm_cr(outputDir,options)
%TEST_IPM_CR Independent Cartesian/annular manufactured-solution checks.
%   test_ipm_cr(outputDir) writes report.json only after every assertion has
%   passed. Optional options.degrees changes the p-refinement levels;
%   options.includeFarTail=false skips the separately labelled R=1e10 case.
%   options.includeIsolatedRefinement=true also varies each patch axis alone.
%   Exact source, velocity and Laplacian formulas do not use either solver.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir)
    outputDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
end
if nargin<2,options=struct();end
degrees=option(options,'degrees',[16,28,40]);
assert(isnumeric(degrees)&&numel(degrees)>=2&&all(degrees>=8)&& ...
    all(degrees==fix(degrees))&&all(diff(degrees)>0), ...
    'ipm_cr:test:Degrees','Use increasing integer polynomial degrees.');
timer=tic;records=struct([]);convergence=struct([]);
names={'mixed_gaussian','mixed_rational','exterior_gaussian'};
for caseIndex=1:numel(names)
    name=names{caseIndex};rows=struct([]);
    for degree=degrees
        cfg=settings(degree);
        started=tic;solver=ipm_cr.psolver(cfg);
        [solution,info]=solver.solve(@(x,y)source(name,x,y), ...
            @(x,y)streamfunction(name,x,y));
        errors=patch_errors(solver,solution,name);
        assert(all(isfinite([errors.core.psi,errors.core.velocity, ...
            errors.core.laplacian,errors.polar.psi,errors.polar.velocity, ...
            errors.polar.laplacian])), 'ipm_cr:test:Finite', ...
            'Every patch and derivative must be finite off grid.');
        assert(info.linearBackwardError<1e-10&&info.traceResidual<1e-10, ...
            'ipm_cr:test:Algebraic','The direct coupled system must solve accurately.');
        assert(all(isfinite([info.overlapValueMismatch,info.overlapVelocityMismatch, ...
            info.normalFluxMismatch,info.strongResidual,info.sourceScale])), ...
            'ipm_cr:test:Audit','Independent interface and PDE audits must be finite.');
        row=struct('name',name,'coreDegree',degree, ...
            'polarDegree',cfg.polarDegree,'angularDegree',cfg.angularDegree, ...
            'rOuter',cfg.rOuter,'seconds',toc(started),'errors',errors, ...
            'diagnostics',info,'counters',solver.Counters);
        rows=append(rows,row);records=append(records,row);
        fprintf('ipm_cr %-18s n=%d psi=%.3e u=%.3e lap=%.3e %.3fs\n', ...
            name,degree,max(errors.core.psi,errors.polar.psi), ...
            max(errors.core.velocity,errors.polar.velocity), ...
            max(errors.core.laplacian,errors.polar.laplacian),row.seconds);
    end
    psiErrors=arrayfun(@(r)max(r.errors.core.psi,r.errors.polar.psi),rows);
    velocityErrors=arrayfun(@(r)max(r.errors.core.velocity,r.errors.polar.velocity),rows);
    finest=rows(end).errors;
    assert(max(finest.core.psi,finest.polar.psi)<1e-8&& ...
        max(finest.core.velocity,finest.polar.velocity)<1e-7, ...
        'ipm_cr:test:ManufacturedAccuracy', ...
        '%s finest off-grid psi/u errors are %.6g/%.6g.',name,psiErrors(end),velocityErrors(end));
    assert(improves(psiErrors)&&improves(velocityErrors), ...
        'ipm_cr:test:Refinement','%s does not improve under independent p refinement.',name);
    convergence=append(convergence,struct('name',name,'degrees',degrees, ...
        'psiInfinityErrors',psiErrors,'velocityInfinityErrors',velocityErrors, ...
        'psiErrorRatios',psiErrors(2:end)./max(psiErrors(1:end-1),realmin), ...
        'velocityErrorRatios',velocityErrors(2:end)./max(velocityErrors(1:end-1),realmin)));
end

% The last case has its Gaussian centre strictly outside the core. Retain
% the independently known nonzero induced core field in the evidence.
probe=manufacture('exterior_gaussian',0,.4);
assert(abs(probe.psi)>1e-4,'ipm_cr:test:ExteriorForcing', ...
    'The exterior-source manufacture must exert a measurable core field.');
queryX=[-.91,-.37,.13,.79;-.77,-.19,.31,.98];
queryY=[.07,.23,.49,.81;.13,.37,.71,1.03];
query=solver.evaluate(solution,queryX,queryY);
exact=manufacture('exterior_gaussian',queryX,queryY);
assert(isequal(size(query.psi),size(queryX))&& ...
    max(abs(query.psi-exact.psi),[],'all')<1e-8, ...
    'ipm_cr:test:QueryShape','Public evaluation must retain same-shape query arrays.');

strictCfg=settings(degrees(end));strictCfg.requireAccuracy=true;
strictSolver=ipm_cr.psolver(strictCfg);
[strictSolution,strictInfo]=strictSolver.solve(@(x,y)source('mixed_rational',x,y), ...
    @(x,y)streamfunction('mixed_rational',x,y));
strictErrors=patch_errors(strictSolver,strictSolution,'mixed_rational');
assert(strictInfo.accepted,'ipm_cr:test:StrictAcceptance', ...
    'A finest manufactured solve must pass the default accuracy gates.');
coarseCfg=settings(8);coarseCfg.requireAccuracy=true;
coarseSolver=ipm_cr.psolver(coarseCfg);
expect_error(@()coarseSolver.solve(@(x,y)source('mixed_rational',x,y), ...
    @(x,y)streamfunction('mixed_rational',x,y)),'ipm_cr:AccuracyGate');
badSource=@(x,y)invalid_off_grid(x,y,size(strictSolver.Core.X),size(strictSolver.Polar.X));
expect_error(@()strictSolver.solve(badSource,zeros(1,strictSolver.Polar.nt)),'ipm_cr:SourceAudit');

% Cache ownership is explicit. A repeated source must reuse the geometry;
% invalid rebuilds must preserve the accepted solver and old solutions.
before=solver.Counters;
[warm,warmInfo]=solver.solve(@(x,y)source('exterior_gaussian',x,y), ...
    @(x,y)streamfunction('exterior_gaussian',x,y));
after=solver.Counters;
assert(after.builds==before.builds&&after.solves==before.solves+1&& ...
    after.factorReuses>before.factorReuses&& ...
    max(abs(warm.core.psi-solution.core.psi),[],'all')<1e-13&& ...
    max(abs(warm.polar.psi-solution.polar.psi),[],'all')<1e-13, ...
    'ipm_cr:test:FactorReuse','A warm solve must reuse the prepared factors.');
beforeInvalid=solver.evaluate(warm,queryX,queryY);builds=after.builds;
expect_error(@()solver.reconfigure(struct('rInner',2)),'ipm_cr:Overlap');
afterInvalid=solver.evaluate(warm,queryX,queryY);
assert(solver.Counters.builds==builds&&isequaln(beforeInvalid,afterInvalid), ...
    'ipm_cr:test:ReconfigureRollback','A failed configuration must preserve the active solver.');
solver.reconfigure(struct('rInner',.9));
expect_error(@()solver.evaluate(warm,queryX,queryY),'ipm_cr:StaleSolution');
[differentOverlap,overlapInfo]=solver.solve(@(x,y)source('exterior_gaussian',x,y), ...
    @(x,y)streamfunction('exterior_gaussian',x,y));
overlapErrors=patch_errors(solver,differentOverlap,'exterior_gaussian');
assert(max(overlapErrors.core.psi,overlapErrors.polar.psi)<1e-8&& ...
    max(overlapErrors.core.velocity,overlapErrors.polar.velocity)<1e-7, ...
    'ipm_cr:test:OverlapIndependence','Changing the overlap must retain manufactured accuracy.');

% A constant density has exactly zero continuous horizontal derivative.
% Exercise both the function and native-patch velocity entry points.
zeroBoundary=zeros(1,solver.Polar.nt);
[constant,constantInfo]=solver.velocity(@(x,y)7+zeros(size(x)),zeroBoundary);
constantError=max(abs([constant.core.psi(:);constant.polar.psi(:)]));
zeroInput=struct('core',zeros(size(solver.Core.X)), ...
    'polar',zeros(size(solver.Polar.X)));
[zero,zeroInfo]=solver.velocity(zeroInput,zeroBoundary);
zeroError=max(abs([zero.core.psi(:);zero.polar.psi(:)]));
zeroValue=solver.evaluate(zero,queryX,queryY);
assert(constantError<1e-10&&zeroError==0&& ...
    all(zeroValue.u1==0,'all')&&all(zeroValue.u2==0,'all'), ...
    'ipm_cr:test:ConstantDensity','Constant/zero density may not generate Darcy flow.');
% Nonconstant density-to-source-to-velocity path, with a known primitive of
% q=-Delta(y exp(-x^2-y^2)); both patch derivatives are exercised directly.
density=@(x,y)y.*exp(-y.^2).*((6-4*y.^2)*sqrt(pi)/2.*erf(x)+2*x.*exp(-x.^2));
psiExact=@(x,y)y.*exp(-x.^2-y.^2);
[native,nativeInfo]=solver.velocity(density,psiExact);
[shifted,~]=solver.velocity(@(x,y)density(x,y)+7,psiExact);
[xx,yy]=meshgrid(linspace(-1.17,1.17,29),linspace(0,2.4,31));
nv=solver.evaluate(native,xx,yy);shiftedValue=solver.evaluate(shifted,xx,yy);
E=exp(-xx.^2-yy.^2);
nativeVelocityError=max([abs(nv.u1+(1-2*yy.^2).*E),abs(nv.u2+2*xx.*yy.*E)],[],'all');
backgroundError=max([abs(nv.u1-shiftedValue.u1),abs(nv.u2-shiftedValue.u2)],[],'all');
assert(nativeVelocityError<1e-6&&backgroundError<1e-10,'ipm_cr:test:DensityDerivative', ...
    'Native density-source velocity errors %.3g; additive background %.3g.',nativeVelocityError,backgroundError);
management=struct('cacheBefore',before,'cacheAfter',after, ...
    'warmDiagnostics',warmInfo,'invalidReconfigureRolledBack',true, ...
    'staleSolutionRejected',true,'differentOverlapErrors',overlapErrors, ...
    'differentOverlapDiagnostics',overlapInfo,'constantPsiError',constantError, ...
    'zeroPsiError',zeroError,'constantDiagnostics',constantInfo,'zeroDiagnostics',zeroInfo, ...
    'coarseAccuracyGateRejected',true,'invalidOffGridSourceRejected',true, ...
    'nativeDensityVelocityError',nativeVelocityError,'additiveBackgroundVelocityDifference',backgroundError, ...
    'nativeDensityDiagnostics',nativeInfo, ...
    'strictDefaultGatesPassed',true,'strictErrors',strictErrors,'strictDiagnostics',strictInfo);

isolated=struct('tested',false);
if option(options,'includeIsolatedRefinement',false)
    axisDegrees=option(options,'isolatedDegrees',[8,14,20]);
    assert(isnumeric(axisDegrees)&&numel(axisDegrees)==3&& ...
        all(axisDegrees>=6)&&all(axisDegrees==fix(axisDegrees))&&all(diff(axisDegrees)>0), ...
        'ipm_cr:test:IsolatedDegrees','Use three increasing isolated polynomial degrees.');
    axisNames={'core','radial','angular'};axisRecords=struct([]);
    for axisIndex=1:numel(axisNames)
        axisName=axisNames{axisIndex};axisRows=struct([]);
        for degree=axisDegrees
            cfg=settings(degrees(end));
            switch axisName
                case 'core',cfg.coreDegree=[degree,degree];
                case 'radial',cfg.polarDegree=degree;
                case 'angular',cfg.angularDegree=degree;
            end
            started=tic;axisSolver=ipm_cr.psolver(cfg);
            [axisSolution,axisInfo]=axisSolver.solve(@(x,y)source('mixed_rational',x,y), ...
                @(x,y)streamfunction('mixed_rational',x,y));
            axisErrors=patch_errors(axisSolver,axisSolution,'mixed_rational');
            axisRows=append(axisRows,struct('degree',degree,'settings',cfg, ...
                'errors',axisErrors,'diagnostics',axisInfo,'seconds',toc(started)));
            fprintf('ipm_cr isolated %-7s n=%d psi=%.3e u=%.3e\n',axisName,degree, ...
                max(axisErrors.core.psi,axisErrors.polar.psi), ...
                max(axisErrors.core.velocity,axisErrors.polar.velocity));
        end
        ep=arrayfun(@(r)max(r.errors.core.psi,r.errors.polar.psi),axisRows);
        eu=arrayfun(@(r)max(r.errors.core.velocity,r.errors.polar.velocity),axisRows);
        assert((ep(end)<.5*ep(1)&&ep(end)<.8*ep(end-1))||ep(end)<1e-9, ...
            'ipm_cr:test:IsolatedPsiRefinement','Isolated %s psi refinement failed.',axisName);
        assert((eu(end)<.5*eu(1)&&eu(end)<.8*eu(end-1))||eu(end)<1e-8, ...
            'ipm_cr:test:IsolatedVelocityRefinement','Isolated %s velocity refinement failed.',axisName);
        axisRecords=append(axisRecords,struct('axis',axisName,'degrees',axisDegrees, ...
            'psiInfinityErrors',ep,'velocityInfinityErrors',eu,'records',axisRows));
    end
    isolated=struct('tested',true,'manufacture','mixed_rational','records',axisRecords);
end

farTail=struct('tested',false);
if option(options,'includeFarTail',true)
    cfg=settings(degrees(end));cfg.outerBoundary='decay';cfg.rOuter=1e10;
    cfg.polarDegree=max(12,ceil(degrees(end)/3));
    started=tic;farSolver=ipm_cr.psolver(cfg);
    [far,farInfo]=farSolver.solve(@(x,y)source('mixed_rational',x,y));
    nearErrors=patch_errors(farSolver,far,'mixed_rational');
    radii=logspace(1,10,31);angles=[.19,.61,1.17,1.91,2.73]';
    x=cos(angles)*radii;y=sin(angles)*radii;
    value=farSolver.evaluate(far,x,y);exact=manufacture('mixed_rational',x,y);
    assert(all(isfinite(value.psi),'all')&&all(isfinite(value.u1),'all')&& ...
        all(isfinite(value.u2),'all')&& ...
        max(nearErrors.core.psi,nearErrors.polar.psi)<1e-8&& ...
        max(nearErrors.core.velocity,nearErrors.polar.velocity)<1e-7, ...
        'ipm_cr:test:LargeDomain','The complete R=1e10 solve must retain core accuracy.');
    farTail=struct('tested',true,'rOuter',cfg.rOuter,'seconds',toc(started), ...
        'nearFieldErrors',nearErrors,'diagnostics',farInfo, ...
        'sampleRadii',radii,'sampleAngles',angles', ...
        'tailPsiAbsoluteError',max(abs(value.psi-exact.psi),[],'all'), ...
        'tailVelocityAbsoluteError',max([abs(value.u1-exact.u1),abs(value.u2-exact.u2)],[],'all'), ...
        'tailPsiRelativeErrors',abs(value.psi-exact.psi)./max(abs(exact.psi),realmin), ...
        'tailRelativeAccuracyAsserted',false, ...
        'boundaryQualification','Homogeneous modal decay at finite R; algebraic source continues beyond R.');
end
report=struct('kind','ipm_cr_independent_elliptic_manufactured_suite_v1', ...
    'passed',true,'seconds',toc(timer),'refinementDegrees',degrees, ...
    'records',records,'convergence',convergence,'management',management, ...
    'farTail',farTail,'isolatedRefinement',isolated,'exteriorInducedCorePsi',probe.psi, ...
    'manufacturedSourcesIndependentOfDiscreteOperators',true, ...
    'offGridBothPatchesTested',true,'timeEvolutionTested',false, ...
    'trajectoryAccuracyClaimed',false);
if ~exist(outputDir,'dir'),mkdir(outputDir);end
path=fullfile(outputDir,'report.json');
fid=fopen(path,'w');assert(fid>=0,'ipm_cr:test:Report','Cannot write %s.',path);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');
fprintf('ipm_cr independent elliptic tests passed in %.3f s; report: %s\n',report.seconds,path);
end

function cfg=settings(degree)
cfg=struct('coreDegree',[degree,degree],'polarDegree',max(8,ceil(degree/2)), ...
    'angularDegree',degree+8,'rOuter',8,'outerBoundary','dirichlet', ...
    'requireAccuracy',false);
end

function q=invalid_off_grid(x,~,coreSize,polarSize)
q=zeros(size(x));
if ~isequal(size(x),coreSize)&&~isequal(size(x),polarSize),q(:)=NaN;end
end

function errors=patch_errors(solver,solution,name)
% Fixed physical queries avoid collocation nodes and interface-only norms.
x=linspace(-1.173,1.169,23);y=linspace(.019,1.181,21);
[X,Y]=meshgrid(x,y);X=X(:);Y=Y(:);
value=ipm_cr.poisson.core.evaluate(solver.Core,solution.core,X,Y);
exact=manufacture(name,X,Y);errors.core=error_norms(value,exact);
r=exp(linspace(log(.913),log(7.891),31));theta=linspace(.017,pi-.013,37)';
X=cos(theta)*r;Y=sin(theta)*r;X=X(:);Y=Y(:);
value=ipm_cr.poisson.polar.evaluate(solver.Polar,solution.polar,X,Y);
exact=manufacture(name,X,Y);errors.polar=error_norms(value,exact);
end

function errors=error_norms(value,exact)
errors=struct('psi',max(abs(value.psi(:)-exact.psi(:))), ...
    'velocity',max([abs(value.u1(:)-exact.u1(:));abs(value.u2(:)-exact.u2(:))]), ...
    'laplacian',max(abs(value.laplacian(:)+exact.source(:))));
end

function result=manufacture(name,x,y)
% Independent continuous formulas; no finite differences or solver calls.
if strcmp(name,'exterior_gaussian')
    a=2.5;b=1;dx=x-a;dy=y-b;E=exp(-dx.^2-dy.^2);
    result=struct('psi',y.*E,'source',(8*y-4*b-4*y.*(dx.^2+dy.^2)).*E, ...
        'u1',-(1-2*y.*dy).*E,'u2',-2*y.*dx.*E);return;
end
beta=1.7;
P=y+.2*x.*y+.3*y.^2+.1*x.^2.*y;
Px=.2*y+.2*x.*y;Py=1+.2*x+.6*y+.1*x.^2;
lapP=.2*y+.6;
if strcmp(name,'mixed_gaussian')
    E=exp(-x.^2-beta*y.^2);
    q=-E.*(lapP-4*(x.*Px+beta*y.*Py)+ ...
        (4*x.^2+4*beta^2*y.^2-2-2*beta).*P);
    result=struct('psi',P.*E,'source',q, ...
        'u1',-(Py-2*beta*y.*P).*E,'u2',(Px-2*x.*P).*E);
elseif strcmp(name,'mixed_rational')
    D=1+x.^2+beta*y.^2;E=D.^(-3);
    q=E.*(-lapP+(12*(x.*Px+beta*y.*Py)+6*(1+beta)*P)./D- ...
        48*(x.^2+beta^2*y.^2).*P./D.^2);
    result=struct('psi',P.*E,'source',q, ...
        'u1',-(Py-6*beta*y.*P./D).*E,'u2',(Px-6*x.*P./D).*E);
else
    error('ipm_cr:test:Manufacture','Unknown manufacture %s.',name);
end
end

function q=source(name,x,y)
e=manufacture(name,x,y);q=e.source;
end

function psi=streamfunction(name,x,y)
e=manufacture(name,x,y);psi=e.psi;
end

function yes=improves(errors)
% Once errors reach roundoff, further monotonic decay is not meaningful.
yes=errors(end)<.25*errors(1)||errors(end)<5e-11;
yes=yes&&(errors(end)<.75*errors(end-1)||errors(end)<5e-11);
end

function expect_error(call,identifier)
try
    call();
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:test:WrongError', ...
        'Expected %s, received %s.',identifier,problem.identifier);return;
end
error('ipm_cr:test:MissingError','Expected error %s.',identifier);
end

function value=option(options,name,fallback)
value=fallback;if isfield(options,name),value=options.(name);end
end

function values=append(values,row)
if isempty(values),values=row;else,values(end+1)=row;end
end
