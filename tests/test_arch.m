function report=test_arch(outputDir)
%TEST_ARCH Public-state identity, context isolation and generic field readout.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','architecture');end
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'architecture.json')),'ipm_cr:ArchOutput','Use a fresh output directory.');
method=struct('mesh',struct('coreSpacing',.06, ...
    'grid',struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
    'poisson',struct('coreDegree',[48,48],'angularDegree',64,'polarDegree',16,'rOuter',1e10));
timer=tic;[s,a]=ipm_cr.init(method);
assert(isstruct(a.mesh)&&isa(a,'ipm_cr.context'),'ipm_cr:ArchMesh');
assert(all(isfield(s,{'W','eta','rhs','scale','step','ctxId','ctxKey'}))&& ...
    ~any(isfield(s,{'rho','residual','flow'})),'ipm_cr:ArchState');
canonical=ipm_cr.config(struct('geometry',method.mesh,'solverOptions',method.poisson));
assert(isequaln(a.opts,canonical),'ipm_cr:ArchConfig');
assert(strcmp(s.ctxId,a.id)&&isequaln(s.ctxKey,a.key),'ipm_cr:ArchInputIdentity');
% The original wall authority needs no transport-grid or state carrier.
cfg=ipm_cr.problem();[wall,scale,refs,guess]=ipm_cr.model.reference(cfg);
assert(isequaln(wall,a.wall)&&isequaln(scale,rmfield(s.scale,'densityOrigin'))&& ...
    isequaln(refs,s.originalAuthority)&&isequal(size(guess),[1,2])&&all(isfinite(guess)), ...
    'ipm_cr:ArchAuthority');
bad=cfg;bad.authoritySha=repmat('0',1,64);
expect(@()ipm_cr.model.reference(bad),'ipm_cr:ReferenceIdentity');
bad=cfg;bad.physics.gaussianWidthX=.95*bad.physics.gaussianWidthX;
expect(@()ipm_cr.model.reference(bad),'polar:T0AnalyticProvenance');
bad=cfg;bad.scaling.referenceWallIntegral=1.1*refs.I0;
expect(@()ipm_cr.model.reference(bad),'polar:T0ReferenceMismatch');
% A valid file SHA cannot turn a late frame into a t0 reference.
loaded=load(cfg.authorityFile,'frame');frame=loaded.frame;frame.tau=1;
path=[tempname(outputDir),'.mat'];save(path,'frame','-v7');
bad=cfg;bad.authorityFile=path;bad.authoritySha=ipm_cr.io.sha(path);
expect(@()ipm_cr.model.reference(bad),'polar:T0AuthorityClock');
x=[0,.3;-.4,.7];y=[0,.1;.3,.4];v=ipm_cr.sample(s,a,x,y);
assert(isequal(size(v.W),size(x))&&all(isfinite(v.Q),'all'),'ipm_cr:ArchSample');
assert(isequaln(v,ipm_cr.diag.fields(s,a.mesh,a.opts,x,y)),'ipm_cr:ArchReadout');
b=ipm_cr.context(a.mesh,a.wall,a.opts);
assert(~strcmp(a.id,b.id)&&isequaln(a.key,b.key),'ipm_cr:ArchContextIdentity');
recovered=ipm_cr.refresh(s,b);assert(isequaln(s.W,recovered.W),'ipm_cr:ArchRefresh');
assert(all(isfield(recovered,{'W','eta','rhs'}))&& ...
    ~any(isfield(recovered,{'rho','residual','flow'}))&& ...
    strcmp(recovered.ctxId,b.id)&&isequaln(recovered.ctxKey,b.key),'ipm_cr:ArchState');
before=b.stats();ownBefore=a.stats();[fa,Wa]=ipm_cr.rhs(s,a);
after=b.stats();ownAfter=a.stats();
assert(isequaln(before,after)&&isequaln(Wa,s.W),'ipm_cr:ArchIsolation', ...
    'Reevaluating one context must preserve the other context and the represented W exactly.');
freshF=field_difference(fa.F,s.rhs.F);
% Direct Newton starts at the stored rates and can change their last bits.
% Match the independent state-integrity fresh-RHS contract, not byte equality.
assert(freshF.relativeInfinity<2e-10,'ipm_cr:ArchFreshRhs', ...
    'Fresh F relative difference %.17g exceeds the existing 2e-10 contract (absolute %.17g).', ...
    freshF.relativeInfinity,freshF.absoluteInfinity);
isolation=struct('otherContextBefore',before,'otherContextAfter',after, ...
    'ownContextBefore',ownBefore,'ownContextAfter',ownAfter,'WBitwiseIdentical',isequaln(Wa,s.W), ...
    'qBitwiseIdentical',isequaln(fa.q,s.rhs.q),'physicalBitwiseIdentical',isequaln(fa.physical,s.rhs.physical), ...
    'FBitwiseIdentical',isequaln(fa.F,s.rhs.F),'freshF',freshF);
% Mathematical input identities are checked independently of time stepping.
changed=a.opts;changed.background.widthX=changed.background.widthX*1.01;
c=ipm_cr.context(a.mesh,a.wall,changed);
assert(~isequaln(s.ctxKey,c.key),'ipm_cr:ArchBackgroundIdentity');
expect(@()ipm_cr.refresh(s,c),'ipm_cr:StateSetup');
% A same-sized grid with different coordinates has a different input key.
changed=a.opts;changed.geometry.coreSpacing=.061;
changed.geometry.coreHalfWidth=1.22;changed.geometry.coreHeight=1.22;
g=ipm_cr.mesh.build(changed.geometry);
c=ipm_cr.context(g,a.wall,changed);
assert(isequal(size(g.core.X),size(a.mesh.core.X))&& ...
    ~isequal(g.core.x,a.mesh.core.x)&&~isequaln(s.ctxKey,c.key),'ipm_cr:ArchGeometryIdentity');
changed=a.opts;changed.transport.smoothSpeedCoefficient=.11;
c=ipm_cr.context(a.mesh,a.wall,changed);
assert(~isequaln(s.ctxKey,c.key),'ipm_cr:ArchTransportIdentity');
% Time policy may vary without changing the mathematical RHS identity.
expectedDt=min([.0005,.95*.2/s.rhs.rate,.001-s.scale.canonicalTime]);
[s,info]=ipm_cr.step(s,a,.001,struct('maxDt',.0005,'cfl',.2));
assert(s.scale.canonicalTime==expectedDt&&info.dt==expectedDt&& ...
    info.cflStepFormulaApplied&&isequal(info.stageCfl,info.dt*info.stageRates)&& ...
    info.realizedMaximumStageCfl==max(info.stageCfl)&& ...
    strcmp(s.ctxId,a.id)&&isequaln(s.ctxKey,a.key), ...
    'ipm_cr:ArchTimeControls');
report=struct('kind','ipm_cr_public_architecture_v1','passed',true, ...
    'meshIsStruct',true,'publicStateUnaliased',true,'inputAndOutputIdentityMatched',true, ...
    'changedRhsKeysDistinguished',true,'changedSetupRefreshRejected',true, ...
    'contextsIsolated',true,'genericQueries',true,'initialCflFormulaExact',true, ...
    'rhsReevaluation',isolation, ...
    'timePolicyIndependent',true,'referenceWithoutTransportGrid',true, ...
    'authorityShaAndDatumRejected',true,'lateReferenceRejected',true,'seconds',toc(timer));
ipm_cr.io.json(fullfile(outputDir,'architecture.json'),report);
fprintf('CR ARCH: state/setup identity, isolated contexts and generic sampling passed (%.3fs)\n',report.seconds);
end

function d=field_difference(a,b)
d=struct('relativeInfinity',0,'absoluteInfinity',0,'patches',struct());
for patch={'core','polar'}
    name=patch{1};absolute=max(abs(a.(name)-b.(name)),[],'all');
    scale=max(1,max(abs(b.(name)),[],'all'));
    d.patches.(name)=struct('absoluteInfinity',absolute,'referenceScale',scale,'relativeInfinity',absolute/scale);
    d.relativeInfinity=max(d.relativeInfinity,absolute/scale);d.absoluteInfinity=max(d.absoluteInfinity,absolute);
end
end

function expect(callback,identifier)
try
    callback();error('ipm_cr:ArchExpectedFailure','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:ArchWrongFailure', ...
        'Expected %s, got %s: %s',identifier,problem.identifier,problem.message);
end
end
