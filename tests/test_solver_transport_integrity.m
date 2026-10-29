function report=test_solver_transport_integrity(outputDir,options)
%TEST_SOLVER_TRANSPORT_INTEGRITY Outer inflow gates and measured cache proposal.
% The derivative-plan candidate is local to this test. A microbenchmark alone
% cannot qualify a production optimization, and never changes mesh.diff.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','solver_transport_integrity');end
if nargin<2,options=struct();end
assert(isstruct(options)&&isscalar(options)&& ...
    all(ismember(fieldnames(options),{'benchmark','iterations','rounds'})), ...
    'ipm_cr:TransportIntegrityOptions','Use benchmark, iterations and rounds.');
benchmark=option(options,'benchmark',true);
validateattributes(benchmark,{'logical'},{'scalar'});
iterations=option(options,'iterations',30);rounds=option(options,'rounds',5);
validateattributes(iterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(rounds,{'numeric'},{'scalar','integer','>=',3});
if ~isfolder(outputDir),mkdir(outputDir);end
timer=tic;threads=maxNumCompThreads(1);
restoreThreads=onCleanup(@()maxNumCompThreads(threads)); %#ok<NASGU>
g=ipm_cr.mesh.build(geometry(.06,193,65));
rho=field(g);z=zero(g);
tr=struct('reconstruction','linear5','smoothSpeedCoefficient',.1);
standalone=ipm_cr.op.prepare(rho,g,tr);
strict=tr;strict.outerInflowPolicy='error';strict.outerInflowTolerance=0;
guarded=ipm_cr.op.prepare(rho,g,strict);
[outX,outY]=radial_velocity(g,.2);
a=standalone(outX,outY);b=guarded(outX,outY);
assert(isequaln(a.B,b.B)&&isequaln(a.rate,b.rate)&& ...
    ~a.outerBoundary.unresolvedInflow&&~b.outerBoundary.unresolvedInflow&& ...
    b.outerBoundary.enforced&&b.outerBoundary.incomingCount==0, ...
    'ipm_cr:TransportIntegrityOutflow','Outflow policy changed the spatial operator or CFL.');
assert(isequaln(a.outerBoundary,a.position.outerBoundary)&& ...
    isequaln(b.outerBoundary,b.position.outerBoundary), ...
    'ipm_cr:TransportIntegrityReport','Complete CFL diagnostics must retain the exterior closure report.');
[inX,inY]=radial_velocity(g,-.2);
incoming=standalone(inX,inY);
assert(incoming.outerBoundary.unresolvedInflow&& ...
    incoming.outerBoundary.incomingCount==g.polar.nt&& ...
    strcmp(incoming.outerBoundary.status,'unprovided_physical_inflow')&& ...
    ~incoming.outerBoundary.physicalInflowDataSupplied&& ...
    ~incoming.outerBoundary.exteriorTailCertified, ...
    'ipm_cr:TransportIntegrityInflow','Reporting mode must expose the missing physical inflow data.');
expect(@()guarded(inX,inY),'ipm_cr:UnprovidedOuterInflow');
probe=guarded(inX,inY,struct('cfl',false,'partials',true));
assert(isequaln(probe.B,incoming.B)&&isempty(probe.rate)&&isempty(probe.position)&& ...
    probe.outerBoundary.unresolvedInflow&&~probe.outerBoundary.enforced&& ...
    ~probe.outerBoundary.completeEvaluation&&~isempty(probe.partials), ...
    'ipm_cr:TransportIntegrityProbe','A private gauge probe must not reject an intermediate incoming speed.');
neutral=guarded(z,z);
assert(neutral.outerBoundary.incomingCount==0&&neutral.rate>0&& ...
    neutral.position.smoothPhysicalSpeedFloor>0, ...
    'ipm_cr:TransportIntegrityFloor','The exterior guard must preserve the LF floor in the complete CFL.');
permissive=strict;permissive.outerInflowTolerance=1e-8;
permitted=ipm_cr.op.prepare(rho,g,permissive);
[smallX,smallY]=radial_velocity(g,-5e-9/g.polar.r(end));
small=permitted(smallX,smallY);
assert(small.outerBoundary.incomingCount>0&& ...
    small.outerBoundary.significantIncomingCount==0&& ...
    ~small.outerBoundary.unresolvedInflow&& ...
    strcmp(small.outerBoundary.status,'inward_speed_within_explicit_tolerance'), ...
    'ipm_cr:TransportIntegrityTolerance','An explicit tolerance must remain visible rather than hide incoming speeds.');
invalid=tr;invalid.outerInflowPolicy='analytic_tail';
expect(@()ipm_cr.op.prepare(rho,g,invalid),'ipm_cr:ExplicitOuterBoundary');
cfg=ipm_cr.config();
assert(strcmp(cfg.transport.outerInflowPolicy,'error')&& ...
    cfg.transport.outerInflowTolerance==0,'ipm_cr:TransportIntegrityDefault', ...
    'Production configuration must reject unprovided physical inflow by default.');
overridden=ipm_cr.config(struct('transport',struct('reconstruction','linear5')));
assert(strcmp(overridden.transport.outerInflowPolicy,'error')&& ...
    overridden.transport.outerInflowTolerance==0, ...
    'ipm_cr:TransportIntegrityDefault','A reconstruction override must not silently remove the production inflow gate.');
report=struct('kind','ipm_cr_transport_integrity_v1','passed',true, ...
    'outflowSpatialOperatorAndCflBitwiseUnchanged',true, ...
    'standaloneDefaultReportsMissingInflow',true, ...
    'productionDefaultRejectsMissingInflow',true,'privateGaugeProbeAllowed',true, ...
    'explicitToleranceAudited',true,'lfFloorAndActiveCflRetained',true, ...
    'outflow',b.outerBoundary,'incoming',incoming.outerBoundary, ...
    'toleratedIncoming',small.outerBoundary,'performance',struct('performed',false), ...
    'derivativeCacheAdopted',false,'trajectoryAccuracyCertified',false);
if benchmark
    report.performance=measure_performance(iterations,rounds);
end
report.seconds=toc(timer);
ipm_cr.io.json(fullfile(outputDir,'solver_transport_integrity.json'),report);
fprintf('CR TRANSPORT: exterior inflow audited/gated, private probes retained, cache adopted=0 (%.3fs)\n',report.seconds);
end

function result=measure_performance(iterations,rounds)
% ABBA blocks retain the same local-reference subtraction and sparse FIND
% order. Operator equality guards preserve diagnostics replacing a matrix.
plan=[.06,193,65;.03,385,129];rows=struct([]);
for level=1:size(plan,1)
    o=geometry(plan(level,1),plan(level,2),plan(level,3));
    g=ipm_cr.mesh.build(o);q=field(g);p=plans(g);
    baseline=ipm_cr.mesh.diff(q,g);candidate=cached_diff(q,g,p);
    assert(isequaln(baseline,candidate),'ipm_cr:TransportIntegrityCacheArithmetic', ...
        'A cache candidate changed the physical difference arithmetic.');
    modified=g;modified.core.Dx=2*g.core.Dx;
    assert(isequaln(ipm_cr.mesh.diff(q,modified),cached_diff(q,modified,p)), ...
        'ipm_cr:TransportIntegrityCacheIdentity','Changed derivative matrices must invalidate a cached plan.');
    reference=@()ipm_cr.mesh.diff(q,g);proposed=@()cached_diff(q,g,p);
    reference();proposed(); % Warm both implementations before measurement.
    samples=zeros(rounds,4);
    for k=1:rounds
        samples(k,1)=block(reference,iterations);
        samples(k,2)=block(proposed,iterations);
        samples(k,3)=block(proposed,iterations);
        samples(k,4)=block(reference,iterations);
    end
    old=median(reshape(samples(:,[1,4]),[],1));
    fresh=median(reshape(samples(:,[2,3]),[],1));
    % This warmed full RHS is a cost reference only; no candidate production
    % dispatch is installed, so its timing is not presented as an A/B gain.
    method=struct('geometry',o,'solverOptions',struct('coreDegree',[48,48], ...
        'angularDegree',64,'polarDegree',16,'rOuter',1e10));
    [state,ctx]=ipm_cr.init(method);ipm_cr.rhs(state,ctx,true);
    warmRhs=zeros(rounds,1);
    for k=1:rounds
        rhsTimer=tic;ipm_cr.rhs(state,ctx,true);warmRhs(k)=toc(rhsTimer);
    end
    row=struct('coreSpacing',plan(level,1),'coreShape',size(g.core.X), ...
        'polarShape',size(g.polar.X),'iterationsPerBlock',iterations, ...
        'abbaSecondsPerCall',samples,'referenceMedianSeconds',old, ...
        'candidateMedianSeconds',fresh,'candidateOverReference',fresh/old, ...
        'microbenchmarkCandidateFaster',fresh<old, ...
        'warmFullRhsSeconds',warmRhs,'warmFullRhsMedianSeconds',median(warmRhs), ...
        'candidateFullRhsMeasured',false,'arithmeticBitwiseIdentical',true, ...
        'staleDerivativePlanRebuilt',true,'productionOptimizationAdopted',false);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
result=struct('performed',true,'threads',1,'rounds',rounds, ...
    'ordering','ABBA reference/candidate/candidate/reference, warmed, same-state', ...
    'records',rows,'candidatePresentOnlyInsideTest',true, ...
    'netProductionSpeedupEstablished',false, ...
    'decision','Keep current mesh.diff; only microbenchmarks compare the cache candidate. A production change additionally requires warmed full-RHS A/B evidence.');
end

function elapsed=block(callback,n)
timer=tic;for j=1:n,callback();end;elapsed=toc(timer)/n;
end
function p=plans(g)
p=struct('coreX',difference_plan(g.core.Dx),'coreY',difference_plan(g.core.Dy), ...
    'polarR',difference_plan(g.polar.Dr),'polarTheta',difference_plan(g.polar.Dt));
end
function p=difference_plan(D)
[row,column,weight]=find(D);
p=struct('operator',D,'row',row,'column',column, ...
    'sum',sparse(row,(1:numel(row))',weight,size(D,1),numel(row)));
end
function d=cached_diff(q,g,p)
d=struct('core',struct('x',apply(g.core.Dx,q.core',p.coreX)', ...
    'y',apply(g.core.Dy,q.core,p.coreY)));
fr=apply(g.polar.Dr,q.polar,p.polarR);ft=apply(g.polar.Dt,q.polar',p.polarTheta)';
d.polar=struct('x',fr.*g.polar.cosTheta-ft.*g.polar.sinTheta./g.polar.r, ...
    'y',fr.*g.polar.sinTheta+ft.*g.polar.cosTheta./g.polar.r);
end
function v=apply(D,q,p)
if ~isequaln(p.operator,D),p=difference_plan(D);end
v=p.sum*(q(p.column,:)-q(p.row,:));
end
function o=geometry(h,nr,nt)
o=struct('coreSpacing',h,'grid',struct('nr',nr,'ntheta',nt,'firstSpacing',h));
end
function q=field(g)
q=struct();
for patch={'core','polar'}
    name=patch{1};p=g.(name);q.(name)=tanh(p.X).*exp(-p.Y);
end
q=ipm_cr.mesh.fill(q,g);
end
function q=zero(g)
q=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
end
function [vx,vy]=radial_velocity(g,coefficient)
vx=struct('core',coefficient*g.core.X,'polar',coefficient*g.polar.X);
vy=struct('core',coefficient*g.core.Y,'polar',coefficient*g.polar.Y);
end
function v=option(s,name,fallback)
if isfield(s,name),v=s.(name);else,v=fallback;end
end
function expect(callback,identifier)
try
    callback();error('ipm_cr:TransportIntegrityExpectedFailure','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:TransportIntegrityWrongFailure', ...
        'Expected %s, got %s: %s',identifier,problem.identifier,problem.message);
end
end
