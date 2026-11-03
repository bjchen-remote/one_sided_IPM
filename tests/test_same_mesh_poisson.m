function report=test_same_mesh_poisson(outputDir,options)
%TEST_SAME_MESH_POISSON Portable bounded same-node Poisson regression tests.
% No private checkpoint/runtime, remeshing, or external data path is used.
% Two tiny public RK steps compare direct/recovered states; the checkpoint,
% pack and context tests use the repository's public t0 authority.
% Four small geometries: L=16/32, uniform/sinh Cartesian core, and the actual
% mapped annular radius (deliberately NOT uniform in log r). Manufactured
% fields have nonzero wall source and a nonzero exact harmonic exterior.
% Active, ghost, full-array, and independent overlap errors remain distinct.
if nargin<2,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','same_mesh_poisson');end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),{'wallSeconds'})), ...
    'ipm_cr:SameMeshTestOptions','Only wallSeconds may be overridden.');
budget=120;if isfield(options,'wallSeconds'),budget=options.wallSeconds;end
validateattributes(budget,{'numeric'},{'real','finite','scalar','positive'});
outputDir=char(java.io.File(char(outputDir)).getCanonicalPath());
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'same_mesh_poisson.json')), ...
    'ipm_cr:SameMeshTestFresh','Choose a fresh result directory.');
timer=tic;inventory=ipm_cr.io.identity('source');
report=struct('kind','portable_same_physical_nodes_poisson_regression_v1', ...
    'status','running','passed',false,'timeAdvanced',0,'records',struct([]), ...
    'convergence',struct([]),'optionRejections',struct(),'context',struct(), ...
    'clock',struct(),'restart',struct(),'sourceInventory',inventory,'sourcesUnchanged',false, ...
    'fourthOrderDesignTarget',true,'requiredTwoLevelOrder',3, ...
    'continuumAccuracyCertified',false,'trajectoryAccuracyCertified',false, ...
    'limitations',{{'The >=3 two-level regression floor detects serious loss of the fourth-order design; it is not a continuum-order certificate.', ...
    'Ghost/full-array errors are recorded separately, without claiming that the active error applies there.', ...
    'The exterior source is exactly zero, but the solver uses its own semidiscrete angular DtN.', ...
    'The native split equivalence uses an exactly reproducible constant nodal source, not a second interpolation of the complete analytic source.', ...
    'Tiny nonzero clock tails are a synthetic persistence fixture; two short direct/recovered steps are not a physical-trajectory accuracy study.'}});
persist(outputDir,report);
try
    levels=[16,32];variants={'uniform','sinh'};modes=[1,2];
    for iv=1:numel(variants)
        for L=levels
            assert(toc(timer)<budget,'ipm_cr:SameMeshTestBudget','Component-test wall budget exceeded.');
            g=ipm_cr.mesh.build(geometry(L,variants{iv},8));
            assert(all(diff(g.polar.r)>0)&&g.polar.r(1)>0,'ipm_cr:SameMeshTestRadius');
            logGap=diff(log(g.polar.r));
            assert(max(logGap)-min(logGap)>1e-6,'ipm_cr:SameMeshTestRadius', ...
                'This regression must exercise a nonuniform log-radius axis.');
            if strcmp(variants{iv},'sinh')
                assert(max(diff(g.core.x))-min(diff(g.core.x))>1e-6, ...
                    'ipm_cr:SameMeshTestMap','The mapped-core case must be nonuniform.');
            end
            obj=ipm_cr.poisson.fd.solver(g,struct());
            for mode=modes
                q=@(x,y)manufacture_component(mode,x,y,'q');
                [flow,info]=obj.solve(q,true,[],[]);
                row=struct('variant',variants{iv},'L',L,'mode',mode, ...
                    'geometry',struct('coreSize',size(g.core.X),'polarSize',size(g.polar.X), ...
                    'minimumLogGap',min(logGap),'maximumLogGap',max(logGap)), ...
                    'analytic',errors(flow,g,mode),'native',struct(),'independentOverlap',struct(), ...
                    'sameStreamfunctionVelocityInfinity',same_psi_velocity(obj,flow,g), ...
                    'splitEquivalenceInfinity',NaN,'lightEquivalenceInfinity',NaN,'information',info);
                ec=manufacture(mode,g.core.X,g.core.Y);ea=manufacture(mode,g.polar.X,g.polar.Y);
                nodal=ipm_cr.mesh.fill(struct('core',ec.q,'polar',ea.q),g);
                [native,~]=obj.solve(nodal,true,[],[]);row.native=errors(native,g,mode);
                % Every local interpolation exactly reproduces a constant.
                % The analytic offset is evaluated on actual Poisson nodes;
                % no reinterpolated analytic full source defines the reference.
                constant=.25;
                residual=struct('core',constant*ones(size(g.core.X)), ...
                    'polar',constant*ones(size(g.polar.X)));
                offset=@(x,y)manufacture_component(mode,x,y,'q')-constant;
                [split,~]=obj.solve(residual,true,[],offset);
                row.splitEquivalenceInfinity=flow_distance(flow,split);
                [light,lightInfo]=obj.solve(q,false,[],[]);
                row.lightEquivalenceInfinity=flow_distance(flow,light);
                row.lightAccepted=lightInfo.accepted;
                x=obj.auditMaps.x;y=obj.auditMaps.y;
                c=obj.evaluate(flow,x,y,'core');a=obj.evaluate(flow,x,y,'polar');exact=manufacture(mode,x,y);
                row.independentOverlap=struct( ...
                    'psiErrorInfinity',maximum([c.psi(:)-exact.psi(:);a.psi(:)-exact.psi(:)]), ...
                    'uErrorInfinity',maximum([c.u1(:)-exact.u1(:);c.u2(:)-exact.u2(:); ...
                    a.u1(:)-exact.u1(:);a.u2(:)-exact.u2(:)]), ...
                    'patchPsiDifferenceInfinity',maximum(c.psi(:)-a.psi(:)), ...
                    'patchUDifferenceInfinity',maximum([c.u1(:)-a.u1(:);c.u2(:)-a.u2(:)]), ...
                    'sampleCount',numel(x));
                f0=manufacture(mode,0,0);far=manufacture(mode,4,1);
                row.manufactureChecks=struct('wallPsi',f0.psi,'nonzeroWallSource',f0.q, ...
                    'exteriorSource',far.q,'exteriorPsi',far.psi, ...
                    'exteriorDtnDefect',4*far.psi_x+far.psi_y+mode*far.psi, ...
                    'laplacianIdentity',maximum([ec.q(:)+ec.psi_xx(:)+ec.psi_yy(:); ...
                    ea.q(:)+ea.psi_xx(:)+ea.psi_yy(:)]));
                report.records=append(report.records,row);report.elapsedSeconds=toc(timer);persist(outputDir,report);
                assert(info.accepted&&info.linearSolveBackwardError<=1e-10, ...
                    'ipm_cr:SameMeshTestAlgebraic','Audited valid source did not pass original-row/trace checks.');
                assert(isnan(lightInfo.accepted),'ipm_cr:SameMeshTestLight', ...
                    'An unaudited solve must not advertise complete acceptance.');
                scale=max([1;abs(flow.core.psi(:));abs(flow.polar.psi(:)); ...
                    abs(flow.core.u1(:));abs(flow.core.u2(:));abs(flow.polar.u1(:));abs(flow.polar.u2(:))]);
                assert(row.sameStreamfunctionVelocityInfinity<=2e-9*scale&& ...
                    row.splitEquivalenceInfinity<=2e-9*scale&&row.lightEquivalenceInfinity<=2e-9*scale, ...
                    'ipm_cr:SameMeshTestRepresentation','Velocity, analytic offset or light path changed the represented streamfunction.');
                assert(f0.psi==0&&abs(f0.q)>.1&&far.q==0&&abs(far.psi)>1e-5&& ...
                    abs(row.manufactureChecks.exteriorDtnDefect)<1e-13, ...
                    'ipm_cr:SameMeshTestManufacture','The nonzero wall source or harmonic exterior test was lost.');
            end
            clear obj
        end
    end
    for variant=variants
        for mode=modes
            rr=report.records(strcmp({report.records.variant},variant{1})&[report.records.mode]==mode);
            [~,order]=sort([rr.L]);rr=rr(order);
            p=log([rr(1).analytic.activePsi/rr(2).analytic.activePsi, ...
                rr(1).analytic.activeU/rr(2).analytic.activeU])/log(2);
            row=struct('variant',variant{1},'mode',mode,'levels',levels, ...
                'analyticActivePsi',[rr(1).analytic.activePsi,rr(2).analytic.activePsi], ...
                'observedPsiOrder',p(1),'observedUOrder',p(2));
            row.analyticActiveU=[rr(1).analytic.activeU,rr(2).analytic.activeU];
            row.regressionFloorPassed=all(p>=report.requiredTwoLevelOrder);
            report.convergence=append(report.convergence,row);
        end
    end
    report.optionRejections=option_tests(g);
    assert(toc(timer)<budget,'ipm_cr:SameMeshTestBudget');
    report.context=context_tests();
    report.clock=clock_tests();
    assert(toc(timer)<budget,'ipm_cr:SameMeshTestBudget');
    report.restart=restart_tests(outputDir,budget-toc(timer));
    report.timeAdvanced=report.restart.endpointTau-report.restart.startTau;
    assert(toc(timer)<budget,'ipm_cr:SameMeshTestBudget');
    report.sourcesUnchanged=isequaln(inventory,ipm_cr.io.identity('source'));
    persist(outputDir,report);
    assert(all([report.convergence.regressionFloorPassed]),'ipm_cr:SameMeshTestConvergence', ...
        'At least one analytic active psi/U pair lost the >=3 two-level regression floor; inspect the preserved raw errors.');
    assert(report.sourcesUnchanged,'ipm_cr:SameMeshTestSources');
    report.passed=true;report.status='passed_bounded_component_and_identity_checks';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function o=geometry(L,variant,R)
a=1.2;o=struct('coreHalfWidth',a,'coreHeight',a,'coreSpacing',a/L, ...
    'polarStart',.65,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',4*L+1,'ntheta',2*L+1,'rmax',R,'firstSpacing',a/(2*L), ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01,'stencilWidth',7));
if strcmp(variant,'sinh')
    o.coreMap=struct('enabled',true,'intervals',[2*L,L],'minSpacing',[a/L,a/L]/sinh(1));
end
end

function m=errors(flow,g,mode)
names={'activePsi','activeU','ghostPsi','ghostU','fullPsi','fullU'};
for k=1:numel(names),m.(names{k})=0;end
for patch={'core','polar'}
    key=patch{1};p=g.(key);f=manufacture(mode,p.X,p.Y);v=flow.(key);
    ep=abs(v.psi-f.psi);eu=max(abs(v.u1-f.u1),abs(v.u2-f.u2));
    m.activePsi=max(m.activePsi,max(ep(p.active)));m.activeU=max(m.activeU,max(eu(p.active)));
    m.ghostPsi=max(m.ghostPsi,max(ep(p.ghost)));m.ghostU=max(m.ghostU,max(eu(p.ghost)));
    m.fullPsi=max(m.fullPsi,max(ep,[],'all'));m.fullU=max(m.fullU,max(eu,[],'all'));
end
end

function d=same_psi_velocity(obj,flow,g)
% Reconstruct derivatives of stored psi independently at its actual nodes.
d=0;
for patch={'core','polar'}
    key=patch{1};p=g.(key);v=obj.evaluate(flow,p.X,p.Y,key);
    d=max(d,maximum([v.u1(:)-flow.(key).u1(:);v.u2(:)-flow.(key).u2(:)]));
end
end

function r=option_tests(g)
r=struct();
r.unknownFDOption=reject(@()ipm_cr.poisson.fd.solver(g,struct('widht',7)),'ipm_cr:FDOptions');
r.unsupportedWidth=reject(@()ipm_cr.poisson.fd.solver(g,struct('width',5)),'ipm_cr:FDOptions');
r.unsupportedOuterBoundary=reject(@()ipm_cr.poisson.fd.solver(g,struct('outerBoundary','dirichlet')),'ipm_cr:FDOptions');
r.unknownBackend=reject(@()ipm_cr.config(struct('poissonBackend','same_mesh_typo')),'');
end

function r=context_tests()
% Geometry reaches the original wall authority; no PDE/RK stage is evolved.
overrides=struct('geometry',geometry(16,'uniform',1e10), ...
    'solverOptions',struct('coreDegree',[12,12],'angularDegree',16, ...
    'polarDegree',8,'rOuter',1e10,'requireAccuracy',false));
[defaults,cfg]=ipm_cr.config(overrides);
assert(strcmp(defaults.poissonBackend,'fem'),'ipm_cr:SameMeshTestDefault');
roundtrip=ipm_cr.config(defaults);
assert(isequaln(defaults,roundtrip),'ipm_cr:SameMeshTestDefault','Resolved defaults must round-trip.');
explicit=defaults;explicit.poissonBackend='fem';
fd=defaults;fd.poissonBackend='same_mesh_fd';fd.fdOptions=struct();
reference=defaults;reference.poissonBackend='spectral';
g=ipm_cr.mesh.build(defaults.geometry);wall=ipm_cr.model.reference(cfg);
c0=ipm_cr.context(g,wall,defaults);c1=ipm_cr.context(g,wall,explicit);cf=ipm_cr.context(g,wall,fd);
cs=ipm_cr.context(g,wall,reference);
assert(isequaln(c0.key,c1.key)&&~isequaln(c0.key,cf.key)&& ...
    ~isequaln(c0.key,cs.key)&&~isequaln(cf.key,cs.key), ...
    'ipm_cr:SameMeshTestBackendKey','The backend must participate in mathematical context identity.');
source=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
[v0,~]=c0.velocity(source,[],true);[v1,~]=c1.velocity(source,[],true);
assert(flow_distance(v0,v1)==0,'ipm_cr:SameMeshTestDefault','Explicit FEM changed the default field.');
[vf,fi]=cf.velocity(source,[],true);assert(flow_distance(v0,vf)<1e-13,'ipm_cr:SameMeshTestZero');
[vs,~]=cs.velocity(source,[],true);
assert(flow_distance(vs,vf)<1e-13,'ipm_cr:SameMeshTestZero', ...
    'The explicitly selected spectral reference and FD must agree on zero source.');
r=struct('defaultIsFem',true,'resolvedRoundtripExact',true, ...
    'explicitFemKeyUnchanged',true,'backendKeyDistinct',true, ...
    'defaultFemFieldExact',true,'spectralReferenceExplicit',true, ...
    'spectralAndFdZeroAgree',true,'backendIdentityComparedIndependently',true, ...
    'fdZeroSourceInformation',fi);
end

function r=restart_tests(outputDir,budget)
% Exercise durable public interfaces, including actual FD/mapped RHS refresh.
% The tiny tail deliberately changes no high/low clock word. This synthetic
% persistence fixture makes no claim to be a separately qualified initial IVP.
timer=tic;
overrides=struct('geometry',geometry(16,'sinh',1e10), ...
    'poissonBackend','same_mesh_fd','fdOptions',struct(), ...
    'solverOptions',struct('rOuter',1e10),'maxDt',1e-4,'timeIntegrator','rk4');
[seed,ctx]=ipm_cr.init(overrides);g=ctx.mesh;method=ctx.opts;
assert(g.core.coreMap.enabled&&strcmp(method.poissonBackend,'same_mesh_fd'), ...
    'ipm_cr:SameMeshTestRestartSetup');
beforeHighLow=[seed.scale.physicalTime,seed.scale.physicalTimeLow];
seed.scale.physicalTimeTail=[2^-110,2^-165];
[seed.rhs,seed.W]=ipm_cr.rhs(seed,ctx,true);
assert(isequal([seed.scale.physicalTime,seed.scale.physicalTimeLow],beforeHighLow), ...
    'ipm_cr:SameMeshTestRestartClock');
folder=fullfile(outputDir,'fd_mapped_tail_restart');mkdir(folder);
file=java.io.RandomAccessFile(fullfile(folder,'owner.lock'),'rw');
channel=file.getChannel();lock=channel.tryLock();
assert(~isempty(lock),'ipm_cr:SameMeshTestRestartOwner');
ownerCleanup=onCleanup(@()release_restart_owner(lock,channel,file)); %#ok<NASGU>
names=ipm_cr.io.history('columns');history=ipm_cr.io.history('row',seed,struct(),0);
checkpoint=struct('state',seed,'geometry',g,'wallOps',ctx.wall,'explicitOptions',method, ...
    'sourceInventory',ctx.identity.source,'authority',ctx.identity.authority, ...
    'history',history,'historyColumns',{names},'auditHistory',struct([]),'accumulatedWallSeconds',0);
[pointer,~]=ipm_cr.io.checkpoint('save',folder,checkpoint,[]);
[loaded,readPointer]=ipm_cr.io.checkpoint('load',folder);
initialSha=strcmp(ipm_cr.io.sha(fullfile(folder,readPointer.path)),readPointer.sha256);
assert(initialSha&&isequaln(loaded.state,seed)&& ...
    isequal(loaded.state.scale.physicalTimeTail,seed.scale.physicalTimeTail), ...
    'ipm_cr:SameMeshTestRestartReadback','Initial FD state/tail or accepted SHA changed in readback.');
recoveredContext=ipm_cr.context(loaded.geometry,loaded.wallOps,loaded.explicitOptions);
assert(isequaln(ctx.key,recoveredContext.key)&& ...
    strcmp(recoveredContext.opts.poissonBackend,'same_mesh_fd'), ...
    'ipm_cr:SameMeshTestRestartContext','Recovered context changed the backend or mapped method.');
direct=ipm_cr.refresh(seed,ctx);recovered=ipm_cr.refresh(loaded.state,recoveredContext);
refreshDifference=restart_difference(direct,recovered);
assert(refreshDifference<=2e-10,'ipm_cr:SameMeshTestRestartRefresh');
assert(toc(timer)<budget,'ipm_cr:SameMeshTestBudget');
target=seed.scale.canonicalTime+1e-4;
[direct,directInfo]=ipm_cr.step(direct,ctx,target);
[recovered,recoveredInfo]=ipm_cr.step(recovered,recoveredContext,target);
assert(direct.step==seed.step+1&&recovered.step==seed.step+1&& ...
    direct.scale.canonicalTime==recovered.scale.canonicalTime, ...
    'ipm_cr:SameMeshTestRestartStep','The two public one-step branches must reach the same accepted endpoint.');
stepDifference=restart_difference(direct,recovered);
assert(stepDifference<=2e-10&&any(recovered.scale.physicalTimeTail~=0), ...
    'ipm_cr:SameMeshTestRestartStep','Recovery changed the short FD step or discarded its nonzero clock tail.');
history(end+1,:)=ipm_cr.io.history('row',recovered,recoveredInfo,0);
checkpoint=loaded;checkpoint.state=recovered;checkpoint.history=history;
[pointer,~]=ipm_cr.io.checkpoint('save',folder,checkpoint,pointer);
[endpoint,readPointer]=ipm_cr.io.checkpoint('load',folder);
endpointSha=strcmp(ipm_cr.io.sha(fullfile(folder,readPointer.path)),readPointer.sha256);
assert(endpointSha&&isequaln(endpoint.state,recovered)&& ...
    isequal(endpoint.state.scale.physicalTimeTail,recovered.scale.physicalTimeTail), ...
    'ipm_cr:SameMeshTestRestartReadback','Endpoint FD state/tail or accepted SHA changed in readback.');
package=ipm_cr.pack(recovered,recoveredContext,folder, ...
    struct('kind','synthetic_nonzero_tail_restart_regression'),history,names);
fieldPath=fullfile(package.folder,package.fieldFile);
metadataPath=fullfile(package.folder,package.metadataFile);
assert(strcmp(ipm_cr.io.sha(fieldPath),package.fieldSHA256)&& ...
    strcmp(ipm_cr.io.sha(metadataPath),package.metadataSHA256), ...
    'ipm_cr:SameMeshTestRestartPack','Published science-package hashes changed.');
packed=load(fieldPath,'data');metadata=jsondecode(fileread(metadataPath));
assert(isequal(packed.data.scale.physicalTimeTail,recovered.scale.physicalTimeTail)&& ...
    isequal(reshape(metadata.scale.physicalTimeTail,1,[]),recovered.scale.physicalTimeTail)&& ...
    strcmp(packed.data.method.poissonBackend,'same_mesh_fd'), ...
    'ipm_cr:SameMeshTestRestartPack','Native MAT/metadata lost clock words or backend identity.');
sample=ipm_cr.diag.fields(recovered,g,method,[0,.2],[0,.1]);
assert(isequal(sample.scalars.clockTail,recovered.scale.physicalTimeTail), ...
    'ipm_cr:SameMeshTestRestartDiagnostics','Public sampled diagnostics lost the extension words.');
summary=ipm_cr.io.history('summary',struct('targetTau',recovered.scale.canonicalTime), ...
    'target_reached',recovered,g,method,history,names,pointer,0,0);
assert(isequal(summary.physicalTimeTail,recovered.scale.physicalTimeTail), ...
    'ipm_cr:SameMeshTestRestartDiagnostics','Public run summary lost the extension words.');
assert(toc(timer)<budget,'ipm_cr:SameMeshTestBudget');
r=struct('syntheticClockTailFixture',true,'backend','same_mesh_fd','sinhCore',true, ...
    'integrationCalls',2,'startTau',seed.scale.canonicalTime,'endpointTau',recovered.scale.canonicalTime, ...
    'initialClockTail',seed.scale.physicalTimeTail,'endpointClockTail',recovered.scale.physicalTimeTail, ...
    'initialCheckpointShaValid',initialSha,'endpointCheckpointShaValid',endpointSha, ...
    'checkpointStateBitwiseExact',true,'contextKeyPreserved',true, ...
    'refreshMixedDifference',refreshDifference,'oneStepMixedDifference',stepDifference, ...
    'maximumStageCfl',max(directInfo.realizedMaximumStageCfl,recoveredInfo.realizedMaximumStageCfl), ...
    'sciencePackageShaValid',true,'matAndJsonClockTailExact',true, ...
    'sampledClockTailExact',true,'summaryClockTailExact',true,'elapsedSeconds',toc(timer));
end

function difference=restart_difference(a,b)
% Compare saved/public physical state, not solver cache/counter metadata.
difference=0;
for name={'W','eta'}
    field=name{1};
    for patch={'core','polar'}
        key=patch{1};difference=max(difference,mixed_difference(a.(field).(key),b.(field).(key)));
    end
end
for name={'F','q','Q'}
    field=name{1};
    for patch={'core','polar'}
        key=patch{1};difference=max(difference,mixed_difference(a.rhs.(field).(key),b.rhs.(field).(key)));
    end
end
for patch={'core','polar'}
    for name={'u1','u2'}
        key=patch{1};field=name{1};
        difference=max(difference,mixed_difference(a.rhs.physical.(key).(field),b.rhs.physical.(key).(field)));
    end
end
for name={'c_l','c_omega','c_r','Uref','delta'}
    key=name{1};difference=max(difference,mixed_difference(a.rhs.(key),b.rhs.(key)));
end
assert(isequaln(a.scale,b.scale),'ipm_cr:SameMeshTestRestartClock', ...
    'The recovered short-step scales and full four-word clock must agree bitwise.');
end

function difference=mixed_difference(a,b)
difference=maximum(a-b)/max(1,maximum(b));
end

function release_restart_owner(lock,channel,file)
lock.release();channel.close();file.close();
end

function r=clock_tests()
old=struct('logC_l',0,'logC_omega',0,'X_shift',0,'canonicalTime',1,'densityOrigin',0, ...
    'physicalTime',1,'physicalTimeLow',2^-54,'physicalTimeTail',[2^-108,2^-162]);
z=[0;0;old.physicalTime;0;1.01;0];increment=2^-140;
[next,info]=ipm_cr.time.stageScale(old,z,increment);
before=[old.physicalTime,old.physicalTimeLow,old.physicalTimeTail];
after=[next.physicalTime,next.physicalTimeLow,next.physicalTimeTail];
% BigDecimal(double) captures each binary64 word exactly. This comparison is
% independent of the solver's TwoSum/compression and needs no extra toolbox.
newDecimal=decimal_sum(after);oldDecimal=decimal_sum(before);
delta=newDecimal.subtract(oldDecimal);
assert(info.strictlyIncreased&&delta.signum()>0&& ...
    abs(delta.doubleValue()-increment)<=1e-10*increment, ...
    'ipm_cr:SameMeshTestClock','The complete clock must retain a positive tiny increment.');
assert(isequal(after(1:2),before(1:2))&&~isequal(after(3:4),before(3:4)), ...
    'ipm_cr:SameMeshTestClock','This regression must actually require extension words.');
[zero,~]=ipm_cr.time.stageScale(old,z,0);
assert(isequal([zero.physicalTime,zero.physicalTimeLow,zero.physicalTimeTail],before), ...
    'ipm_cr:SameMeshTestClock','Zero offset must preserve incoming clock words.');
opts=ipm_cr.config();model=opts.background;model.tailCorrection=true;
% Deliberately amplified, non-normalized words isolate both consumer paths.
% This is an algebraic field-consumer test, not a physical trajectory clock.
s=old;s.physicalTime=0;s.physicalTimeLow=0;s.physicalTimeTail=[.125,-.03125];
folded=s;folded.physicalTime=sum(s.physicalTimeTail);folded.physicalTimeTail=[0,0];
x=[-.4,.2,1.3,3];y=[.2,.7,.4,1];
b=ipm_cr.model.background(model,s,x,y);reference=ipm_cr.model.background(model,folded,x,y);
defect=0;
for name={'H','x','y','xx','xy'}
    key=name{1};defect=max(defect,maximum(b.(key)-reference.(key)));
end
q=ipm_cr.model.source(model,s,x,y);
defect=max(defect,maximum(q-reference.x));
assert(defect<2e-13,'ipm_cr:SameMeshTestClockConsumer', ...
    'Analytic background/source must consume both physicalTimeTail words.');
r=struct('beforeWords',before,'afterWords',after,'requestedIncrement',increment, ...
    'exactWordDifferenceAsDouble',delta.doubleValue(),'strictlyIncreased',true, ...
    'twoWordClockUnchanged',true,'consumerInfinity',defect,'update',info);
end

function value=decimal_sum(words)
value=javaObject('java.math.BigDecimal',0);
for word=words,value=value.add(javaObject('java.math.BigDecimal',word));end
end

function f=manufacture(mode,x,y)
% A C10 radial splice: pure sin(m theta)/r^m outside r=1.8, finite at r=0.
s=x.^2+y.^2;a=.6^2;b=1.8^2;ell=.75;wall=.15;
c=zeros(size(s));e=ones(size(s));cp=c;cpp=c;outside=s>=b;inside=s<=a;
c(outside)=1;e(outside)=0;middle=~(inside|outside);
t=(s(middle)-a)/(b-a);z=1-t;C=zeros(size(t));E=C;
for k=11:21
    coefficient=nchoosek(21,k);
    C=C+coefficient*t.^k.*z.^(21-k);E=E+coefficient*z.^k.*t.^(21-k);
end
c(middle)=C;e(middle)=E;normalizer=21*nchoosek(20,10);
cp(middle)=normalizer*t.^10.*z.^10/(b-a);
cpp(middle)=10*normalizer*t.^9.*z.^9.*(1-2*t)/(b-a)^2;
g=zeros(size(s));gp=g;gpp=g;
v=s(inside)+ell^2;g(inside)=v.^(-mode);gp(inside)=-mode*v.^(-mode-1);gpp(inside)=mode*(mode+1)*v.^(-mode-2);
v=s(outside);g(outside)=v.^(-mode);gp(outside)=-mode*v.^(-mode-1);gpp(outside)=mode*(mode+1)*v.^(-mode-2);
v=s(middle);vi=v+ell^2;gi=vi.^(-mode);go=v.^(-mode);
gip=-mode*vi.^(-mode-1);gop=-mode*v.^(-mode-1);
gipp=mode*(mode+1)*vi.^(-mode-2);gopp=mode*(mode+1)*v.^(-mode-2);
g(middle)=e(middle).*gi+c(middle).*go;
gp(middle)=e(middle).*gip+c(middle).*gop+cp(middle).*(go-gi);
gpp(middle)=e(middle).*gipp+c(middle).*gopp+2*cp(middle).*(gop-gip)+cpp(middle).*(go-gi);
if mode==1,P=y;Px=zeros(size(x));Py=ones(size(x));else,P=2*x.*y;Px=2*y;Py=2*x;end
h=e;hp=-cp;hpp=-cpp;
f.psi=P.*g+wall*y.^2.*h;
f.psi_x=Px.*g+2*x.*P.*gp+2*wall*x.*y.^2.*hp;
f.psi_y=Py.*g+2*y.*P.*gp+wall*(2*y.*h+2*y.^3.*hp);
f.psi_xx=4*x.*Px.*gp+2*P.*gp+4*x.^2.*P.*gpp+wall*y.^2.*(2*hp+4*x.^2.*hpp);
f.psi_yy=4*y.*Py.*gp+2*P.*gp+4*y.^2.*P.*gpp+wall*(2*h+10*y.^2.*hp+4*y.^4.*hpp);
f.q=-4*P.*(s.*gpp+(mode+1)*gp)-wall*(2*h+12*y.^2.*hp+4*s.*y.^2.*hpp);
f.q(outside)=0;f.u1=-f.psi_y;f.u2=f.psi_x;
end

function value=manufacture_component(mode,x,y,name)
f=manufacture(mode,x,y);value=f.(name);
end

function d=flow_distance(a,b)
d=0;for patch={'core','polar'},for field={'psi','u1','u2'}
    key=patch{1};name=field{1};d=max(d,maximum(a.(key).(name)-b.(key).(name)));
end,end
end

function value=maximum(v)
value=max(abs(v),[],'all');
end

function result=reject(call,identifier)
result=struct('rejected',false,'identifier','');
try,call();catch problem,result.rejected=true;result.identifier=problem.identifier;end
assert(result.rejected&&(isempty(identifier)||strcmp(result.identifier,identifier)), ...
    'ipm_cr:SameMeshTestRejection','Expected rejection %s; observed %s.',identifier,result.identifier);
end

function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end

function persist(folder,report)
test_report(folder,'same_mesh_poisson',report);
end
