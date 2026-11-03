function report=test_fem_wall(outputDir,controls)
%TEST_FEM_WALL Static full-domain fourth-order finite-cutoff wall-source test.
% q=1+.1*(1-r^2/R^2)^4 inside R=8, zero outside; q=1 at both wall corners.
% The radial Green oracle is independent of greenPolynomial. complexDilog is
% used only as a low-level special function, checked by frozen SciPy/mpmath
% full-field sentinels and independently integrated source-Green wall values.
if nargin<2,controls=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','fem_wall');end
assert(isstruct(controls)&&isscalar(controls)&&all(ismember(fieldnames(controls),{'wallSeconds'})), ...
    'ipm_cr:FEMWallOptions','Only wallSeconds may be overridden.');
budget=120;if isfield(controls,'wallSeconds'),budget=controls.wallSeconds;end
validateattributes(budget,{'numeric'},{'real','finite','scalar','positive'});
outputDir=char(java.io.File(char(outputDir)).getCanonicalPath());if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'fem_wall.json')),'ipm_cr:FEMWallFresh','Choose a fresh output directory.');
threads=maxNumCompThreads(1);threadCleanup=onCleanup(@()maxNumCompThreads(threads)); %#ok<NASGU>
timer=tic;inventory=ipm_cr.io.identity('source');
defaults=ipm_cr.poisson.fem.options();degree=defaults.degree;
report=struct('kind','portable_fem_nonzero_wall_cutoff_accuracy_v1','status','running','passed',false, ...
    'productionDegree',degree, ...
    'timeAdvanced',0,'source','q=1+.1*(1-r^2/64)^4 for r<=8; zero exterior', ...
    'sourceAtOuterWallCorners',1,'reference',struct(),'negativeControl',struct(), ...
    'records',struct([]),'convergence',struct([]),'sourceInventory',inventory,'sourcesUnchanged',false, ...
    'orderGate',struct('minimum',3.8,'nominal',4,'explanation', ...
    'Fine-pair full velocity order has a 0.2 finite-resolution tolerance; measured order is recorded without rounding.'), ...
    'scope','static source family and stated mesh refinement; not a certificate for arbitrary nonsmooth sources');
persist();
try
    report.reference=reference_checks();persist();
    g=ipm_cr.mesh.build(geometry(16));raw=ipm_cr.poisson.fem.solver(g,struct('degree',4, ...
        'nodalLayout','native','stride',[],'outerBoundary','forced_dtn'));
    bd=boundary('forced_dtn');[f,info]=raw.solve(@source,bd,[]);gate(info);
    native_contract(raw,f,g,info);
    report.negativeControl=struct('level',16,'backend','unenriched native degree4/5', ...
        'errors',metrics(raw,f,g),'information',info);
    assert(report.negativeControl.errors.fullU>1e-4,'ipm_cr:FEMWallNegativeControl', ...
        'The raw finite-cutoff nonzero-wall corner error was not exercised.');
    clear raw f g;persist();
    for L=[64,128]
        check_budget();g=ipm_cr.mesh.build(geometry(L));
        native=ipm_cr.mesh.fill(struct('core',source(g.core.X,g.core.Y),'polar',source(g.polar.X,g.polar.Y)),g);
        for kinds={'forced_dtn','green_dirichlet'}
            kind=kinds{1};check_budget();obj=ipm_cr.poisson.fem.enriched(g,struct('degree',degree, ...
                'quadOrder',8,'nodalLayout','native','stride',[],'outerBoundary',kind));bd=boundary(kind);
            for inputs={'analytic','native'}
                inputKind=inputs{1};q=@source;if strcmp(inputKind,'native'),q=native;end
                check_budget();[flow,info]=obj.solve(q,bd,[]);gate(info);
                layout=native_contract(obj,flow,g,info);
                row=struct('L',L,'boundary',kind,'input',inputKind,'errors',metrics(obj,flow,g), ...
                    'information',info,'nativeLayout',layout,'lastPhysicalElementWidth', ...
                    exp(obj.p.tAxis.breaks(end))-exp(obj.p.tAxis.breaks(end-1)));
                assert(row.errors.wallNormal<1e-10,'ipm_cr:FEMWallNormal');
                if isempty(report.records),report.records=row;else,report.records(end+1)=row;end %#ok<AGROW>
                persist();
            end
            clear obj flow
        end
    end
    for kinds={'forced_dtn','green_dirichlet'}
        for inputs={'analytic','native'}
            kind=kinds{1};inputKind=inputs{1};
            select=strcmp({report.records.boundary},kind)&strcmp({report.records.input},inputKind);
            pair=report.records(select);assert(isequal([pair.L],[64,128]),'ipm_cr:FEMWallPair');
            names=fieldnames(pair(1).errors);orders=struct();
            for j=1:numel(names)
                name=names{j};a=pair(1).errors.(name);b=pair(2).errors.(name);
                orders.(name)=NaN;if a>0&&b>0,orders.(name)=log(a/b)/log(2);end
            end
            value=struct('boundary',kind,'input',inputKind,'levels',[64,128],'orders',orders);
            if isempty(report.convergence),report.convergence=value;else,report.convergence(end+1)=value;end %#ok<AGROW>
            assert(orders.fullU>=3.8&&orders.activeU>=3.8&&orders.cornerU>=3.8&& ...
                pair(2).errors.fullU<1e-6,'ipm_cr:FEMWallFourthOrder', ...
                'Full-domain or corner velocity did not meet the explicit fourth-order regression gate.');
        end
    end
    report.sourcesUnchanged=isequaln(inventory,ipm_cr.io.identity('source'));
    assert(report.sourcesUnchanged,'ipm_cr:FEMWallSources');check_budget();
    report.passed=true;report.status='passed_static_nonzero_wall_full_velocity_fourth_order';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');persist();rethrow(problem);
end
persist();
    function check_budget()
        assert(toc(timer)<budget,'ipm_cr:FEMWallBudget','Static test wall-time budget exceeded.');
    end
    function persist()
        report.elapsedSeconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'fem_wall.json'),report);
    end
end

function g=geometry(L)
a=1.2;g=struct('coreHalfWidth',a,'coreHeight',a,'coreSpacing',a/L, ...
    'polarStart',.65,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',4*L+1,'ntheta',2*L+1,'rmax',8,'firstSpacing',a/(2*L), ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01,'stencilWidth',7));
end
function q=source(x,y)
q=1+.1*(1-(x.^2+y.^2)/64).^4;
end
function bd=boundary(kind)
% The independent oracle includes this source only inside the finite disk.
% Declare that support explicitly so the corner enrichment is applicable.
bd=struct('kind',kind,'exteriorStatus','declared_compact_source','value',@(x,y)zeros(size(x)));
if strcmp(kind,'green_dirichlet'),bd.value=@(x,y)reference_value(x,y,'psi');end
end
function gate(info)
assert(info.solverCompleted&&~info.accuracyAssessed&&info.finiteOutput&& ...
    ~isfield(info,'accepted'),'ipm_cr:FEMWallCompletion');
assert(info.weakBackwardError<1e-10&&info.traceAbsolute<1e-9&& ...
    info.essentialBoundaryAbsolute<1e-9,'ipm_cr:FEMWallOriginalGates', ...
    'Original weak rows, traces and essential boundary gates remain mandatory.');
end
function record=native_contract(obj,flow,g,info)
assert(info.sameNativeGrid&&info.noIndependentPoissonGrid&&strcmp(info.nodalLayout,'native')&& ...
    info.coefficientDof==numel(g.core.X)+numel(g.polar.X)&&g.ghostWidth==3&& ...
    isequal(obj.c.x(:),g.core.x(:))&&isequal(obj.c.y(:),g.core.y(:))&& ...
    isequal(obj.p.r(:),g.polar.r(:))&&isequal(obj.p.theta(:),g.polar.theta(:)), ...
    'ipm_cr:FEMWallNativeGrid','Wall accuracy must use every original transport coordinate.');
assert(isequal(size(flow.core.psi),size(g.core.X))&&isequal(size(flow.polar.psi),size(g.polar.X)));
if ~isfield(flow,'enrichment')
    assert(isequal(flow.core.psi,flow.coefficients.core)&&isequal(flow.polar.psi,flow.coefficients.polar));
end
d=[obj.c.xAxis.elementDegrees(:);obj.c.yAxis.elementDegrees(:); ...
    obj.p.tAxis.elementDegrees(:);obj.p.thetaAxis.elementDegrees(:)];
assert(all(d==info.degree|d==info.degree+1),'ipm_cr:FEMWallNativeDegree');
record=struct('axesExactlyNative',true,'fullGhostLayersRetained',true, ...
    'coefficientDof',info.coefficientDof,'targetDegree',info.degree, ...
    'actualDegreeRange',[min(d),max(d)],'higherDegreeRemainderElements',nnz(d==info.degree+1));
end
function r=metrics(obj,flow,g)
r=struct('activePsi',0,'activeU',0,'fullPsi',0,'fullU',0,'wallPsi',0,'wallU',0, ...
    'wallNormal',0,'outerPsi',0,'outerU',0,'offgridCorePsi',0,'offgridCoreU',0, ...
    'offgridPolarPsi',0,'offgridPolarU',0,'cornerPsi',0,'cornerU',0);
for keys={'core','polar'}
    key=keys{1};p=g.(key);v=flow.(key);truth=radial_reference(p.X,p.Y);
    eu=max(abs(v.u1-truth.u1),abs(v.u2-truth.u2));ep=abs(v.psi-truth.psi);wall=p.Y==0;
    if strcmp(key,'polar'),wall(:,[1,end])=true;end
    r.activePsi=max(r.activePsi,max(ep(p.active)));r.activeU=max(r.activeU,max(eu(p.active)));
    r.fullPsi=max(r.fullPsi,max(ep,[],'all'));r.fullU=max(r.fullU,max(eu,[],'all'));
    r.wallPsi=max(r.wallPsi,max(ep(wall)));r.wallU=max(r.wallU,max(eu(wall)));
    r.wallNormal=max(r.wallNormal,max(abs(v.u2(wall))));
    if strcmp(key,'polar'),r.outerPsi=max(ep(end,:));r.outerU=max(eu(end,:));end
end
[x,y]=meshgrid(linspace(-1.09,1.09,24),linspace(.007,1.09,17));
[r.offgridCorePsi,r.offgridCoreU]=query_errors(obj,flow,x,y,'core');
[t,a]=meshgrid(linspace(.001,pi-.001,61),[.8,1.1,1.6,3,5,7,7.9,8]);
[r.offgridPolarPsi,r.offgridPolarU]=query_errors(obj,flow,a.*cos(t),a.*sin(t),'polar');
[t,a]=meshgrid([0,1e-5,1e-4,1e-3,.01,.025,.05,.1],8*[.96,.98,.99,.995,.999,1]);
[r.cornerPsi,r.cornerU]=query_errors(obj,flow,[a.*cos(t),-a.*cos(t)],[a.*sin(t),a.*sin(t)],'polar');
end
function [p,u]=query_errors(obj,flow,x,y,chart)
v=obj.evaluate(flow,x,y,chart);e=radial_reference(x,y);
p=max(abs(v.psi(:)-e.psi(:)));u=max([abs(v.u1(:)-e.u1(:));abs(v.u2(:)-e.u2(:))]);
end
function value=reference_value(x,y,name)
f=radial_reference(x,y);value=f.(name);
end
function f=radial_reference(x,y)
% For q=rho^(m-2), -Delta[(Re(z^m)-|z|^m)/m^2]=rho^(m-2).
% Odd sine-mode radial moments sum to F_m below. No Green basis is used.
R=8;z=complex(x,y)/R;rho2=abs(z).^2;corner=z==1|z==-1;zz=z;zz(corner)=0;
az=atanh(zz);azOverZ=ones(size(z));nonzero=zz~=0;azOverZ(nonzero)=az(nonzero)./zz(nonzero);
T=.5*(ipm_cr.poisson.fem.complexDilog(z)-ipm_cr.poisson.fem.complexDilog(-z));
f=struct('psi',zeros(size(x)),'u1',zeros(size(x)),'u2',zeros(size(x)));
coefficient=[1.1,-.4,.6,-.4,.1];
for k=0:4
    m=2*k+2;c=coefficient(k+1);S=zeros(size(z));dS=S;ratio=S;
    for j=1:2:m-1,S=S+z.^(m-j)/j;dS=dS+(m-j)*z.^(m-j-1)/j;end
    for j=0:m/2-1,ratio=ratio-z.^(2*j);end
    F=((z.^m-1).*az-S)/m^2-T/m;
    dF=(z.^m-1).*azOverZ/m+(ratio-dS)/m^2;
    dF(corner)=-sum(1./(1:2:m-1))/m;
    f.psi=f.psi+c*R^2*(real(z.^m)-rho2.^(m/2))/m^2-c*2*R^2/pi*imag(F);
    f.u1=f.u1+c*R/m*(imag(z.^(m-1))+rho2.^k.*imag(z))+c*2*R/pi*real(dF);
    f.u2=f.u2+c*R/m*(real(z.^(m-1))-rho2.^k.*real(z))-c*2*R/pi*imag(dF);
end
f.psi(y==0)=0;f.u2(y==0)=0;
end
function r=reference_checks()
% 80-digit mpmath full-field values, independently cross-checked by SciPy
% spence and direct angular-integrated source Green quadrature on the wall.
fixture=[0,0,0,-5.2999101303388745,0; ...
    .3,.2,1.0374378146666243,-5.078209650284359,-.0042709187299703327; ...
    -1,.7,3.4207342285909101,-4.5155540712960232,.047584574474442885; ...
    3,2,7.9461273392820937,-3.0101866494765028,-.36685429492867383; ...
    -4,1,4.2378154403172008,-3.7458526991847672,.26192517762646222; ...
    7.2,.15,.49359774075158669,-3.2182573454132402,-.10507909190033066; ...
    7.999999,1e-6,2.5671791194240267e-6,-2.5671787102688972,-4.925659326517838e-6; ...
    8,0,0,-2.567174284610148,0;-8,0,0,-2.567174284610148,0; ...
    0,8,12.992852153173859,1.4711140070005126,0];
f=radial_reference(fixture(:,1),fixture(:,2));values=[f.psi,f.u1,f.u2];
error=max(abs(values-fixture(:,3:5)),[],'all');
wall=[0,-5.299910130338874;1,-5.26380983114502;4,-4.741825668007871; ...
    7.9,-2.7301394287101814;8,-2.5671742846101475;-8,-2.5671742846101475];
f=radial_reference(wall(:,1),zeros(size(wall,1),1));wallError=max(abs(f.u1-wall(:,2)));
assert(error<1e-12&&wallError<1e-12,'ipm_cr:FEMWallOracle', ...
    'Independent high-precision and source-Green reference fixtures failed.');
r=struct('fullFieldSentinelCount',size(fixture,1),'fullFieldSentinelMaxError',error, ...
    'independentWallGreenCount',size(wall,1),'independentWallGreenMaxError',wallError, ...
    'fixtureTolerance',1e-12,'fullOracleUsesGreenPolynomial',false, ...
    'specialFunctionPrimitive','complexDilog for principal Li2 only');
end
