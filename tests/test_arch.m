function report=test_arch(outputDir)
%TEST_ARCH Public-state identity, context isolation and generic field readout.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','architecture');end
if ~isfolder(outputDir),mkdir(outputDir);end
method=struct('mesh',struct('coreSpacing',.06, ...
    'grid',struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
    'poisson',struct('coreDegree',[48,48],'angularDegree',64,'polarDegree',16,'rOuter',1e10));
timer=tic;[s,a]=ipm_cr.init(method);
assert(isstruct(a.mesh)&&isa(a,'ipm_cr.context'),'ipm_cr:ArchMesh');
assert(all(isfield(s,{'W','eta','rhs','scale','step','ctxId','ctxKey'}))&& ...
    ~any(isfield(s,{'rho','residual','flow'})),'ipm_cr:ArchState');
canonical=ipm_cr.config(struct('geometry',method.mesh,'solverOptions',method.poisson));
assert(isequaln(a.opts,canonical),'ipm_cr:ArchConfig');
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
expect(@()ipm_cr.step(s,b,.001),'ipm_cr:StaleState');
expect(@()ipm_cr.sample(s,b,x,y),'ipm_cr:StaleState');
recovered=ipm_cr.refresh(s,b);assert(isequaln(s.W,recovered.W),'ipm_cr:ArchRefresh');
assert(all(isfield(recovered,{'W','eta','rhs'}))&& ...
    ~any(isfield(recovered,{'rho','residual','flow'})),'ipm_cr:ArchState');
before=b.stats();[fa,Wa]=ipm_cr.rhs(s,a);after=b.stats();
assert(isequaln(before,after)&&isequaln(fa.F,s.rhs.F)&&isequaln(Wa,s.W), ...
    'ipm_cr:ArchIsolation');
% A valid changed background cannot reuse the stored accepted initial RHS.
changed=a.opts;changed.background.widthX=changed.background.widthX*1.01;
c=ipm_cr.context(a.mesh,a.wall,changed);
expect(@()ipm_cr.refresh(s,c),'ipm_cr:StateSetup');
expect(@()ipm_cr.time.step(s,a.mesh,a.wall,changed,.001,a),'ipm_cr:StateContext');
% Geometry identity must also reject a same-sized but differently spaced grid.
changed=a.opts;changed.geometry.coreSpacing=.061;
changed.geometry.coreHalfWidth=1.22;changed.geometry.coreHeight=1.22;
g=ipm_cr.mesh.build(changed.geometry);
expect(@()ipm_cr.time.step(s,g,a.wall,changed,.001,a),'ipm_cr:StateContext');
changed=a.opts;changed.transport.smoothSpeedCoefficient=.11;
expect(@()ipm_cr.time.step(s,a.mesh,a.wall,changed,.001,a),'ipm_cr:StateContext');
% Time policy may vary without changing the mathematical RHS identity.
[s,info]=ipm_cr.step(s,a,.001,struct('maxDt',.0005,'cfl',.2));
assert(s.scale.canonicalTime==.0005&&info.realizedMaximumStageCfl<=.2*(1+128*eps), ...
    'ipm_cr:ArchTimeControls');
report=struct('kind','ipm_cr_public_architecture_v1','passed',true, ...
    'meshIsStruct',true,'publicStateUnaliased',true,'staleContextRejected',true, ...
    'changedRhsRejected',true,'contextsIsolated',true,'genericQueries',true, ...
    'timePolicyIndependent',true,'referenceWithoutTransportGrid',true, ...
    'authorityShaAndDatumRejected',true,'lateReferenceRejected',true,'seconds',toc(timer));
ipm_cr.io.json(fullfile(outputDir,'architecture.json'),report);
fprintf('CR ARCH: state/setup identity, isolated contexts and generic sampling passed (%.3fs)\n',report.seconds);
end

function expect(callback,identifier)
try
    callback();error('ipm_cr:ArchExpectedFailure','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:ArchWrongFailure', ...
        'Expected %s, got %s: %s',identifier,problem.identifier,problem.message);
end
end
