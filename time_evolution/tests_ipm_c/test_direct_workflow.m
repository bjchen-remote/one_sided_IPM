function report=test_direct_workflow(folder)
%TEST_DIRECT_WORKFLOW Actual CFL-limited explicit steps, checkpoint and replay.
assert(startsWith(folder,filesep)&&~isfolder(folder));mkdir(folder);
% Use a resolved production geometry so this lifecycle test also obeys the
% independent spatial-defect gate, rather than a 65-point stretched tail.
cfg=ipm_c.config();cfg.grid.nr=385;cfg.grid.ntheta=128;
cfg.remesh.minimumFirstSpacing=cfg.grid.firstSpacing;
cfg.targetTau=.005;cfg.time.maxDt=.2;cfg.maxSteps=2;cfg.wallSeconds=90;
cfg.runtime.threads=1;cfg.output.verbose=false;
direct=fullfile(folder,'direct');r=ipm_c.solve(cfg,direct);
assert(strcmp(r.status,'target_reached')&&r.steps==2&&~r.exponentialAuxiliarySystemUsed&& ...
    ~r.nonlinearStageSolverUsed&&~r.stepRetryUsed, ...
    'ipm_c:SharedNativeSteps');
[reference,~]=ipm_c.output.loadCheckpoint(direct);
assert(isfield(reference.state,'rho')&&~isfield(reference.state,'v')&& ...
    ~isfield(reference.state,'ops')&&~isfield(reference.state,'rhsCache'),'ipm_c:SharedPrimaryState');
split=fullfile(folder,'split');short=cfg;short.maxSteps=1;
first=ipm_c.solve(short,split);assert(first.steps==1&&strcmp(first.status,'segment_budget'));
assert(strcmp(first.lastTimeStep.activeLimiter,'cfl')&& ...
    abs(first.lastTimeStep.dt*first.lastTimeStep.rate-cfg.time.cfl)<1e-12, ...
    'ipm_c:DirectActualCflSelection');
save(fullfile(folder,'resume_request.mat'),'cfg','split','-v7.3');
% A same-sequence publication must always leave the pointed bytes intact.
[state,context]=ipm_c.initialize(cfg);authority=struct('cart',context.cart,'shared',context.shared, ...
    'authoritySha256',ipm_c.output.sha256(cfg.authorityFile));
io=fullfile(folder,'io_fixture');mkdir(io);history=zeros(1,13);history(1,4)=NaN;
p0=ipm_c.output.saveCheckpoint(io,state,authority,cfg,history,NaN,0);
p1=ipm_c.output.saveCheckpoint(io,state,authority,cfg,history,NaN,0);
assert(~strcmp(p0.path,p1.path)&&strcmp(ipm_c.output.sha256(p0.path),p0.sha256), ...
    'ipm_c:SharedRepeatedPublication');
p2=ipm_c.output.saveCheckpoint(io,state,authority,cfg,history,NaN,0);
assert(~strcmp(p1.path,p2.path)&&strcmp(ipm_c.output.sha256(p1.path),p1.sha256)&& ...
    numel(dir(fullfile(io,'rolling_*.mat')))==2,'ipm_c:SharedTwoSlots');
ipm_c.output.loadCheckpoint(io);
original=ipm_c.output.frame(io,state,history);
reconstructed=state;reconstructed.flow.fullQ=reconstructed.flow.fullQ+1e-12;
retained=ipm_c.output.frame(io,reconstructed,history);
index=ipm_c.data(fullfile(io,'data.mat'));
assert(isequaln(retained.Q,original.Q)&&isequaln(index.rho_xt,original.Q), ...
    'ipm_c:SharedRetainedFrameIndex','The latest index must use the actual retained immutable frame.');
remesh=test_remesh_transaction(state);
data=ipm_c.data(fullfile(direct,'data.mat'));
frame=load(fullfile(direct,data.latestFrame),'frame');
ct=cos(data.theta);ct([1,end])=[1,-1];
if mod(numel(ct),2)==1,ct((numel(ct)+1)/2)=0;end
assert(isequaln(data.rho_xt,frame.frame.Q)&& ...
    isequaln(data.X,data.r*ct)&& ...
    all(data.Y(:,[1,end])==0,'all')&&data.grid.originIsSinglePhysicalPoint, ...
    'ipm_c:SharedPolarOutput');
report=struct('kind','direct_cfl_polar_native_two_steps_and_io_v3','passed',true, ...
    'native',r,'splitFirst',first,'resumePending',true,'sameSequenceSlotSafety',true, ...
    'polarOutputCoordinatesVerified',true,'retainedFrameIndexVerified',true,'automaticRemeshPreserved',true,'remesh',remesh,'clockKind','canonical_tau', ...
    'spatialTrajectoryQualified',false);
ipm_c.output.atomicJson(fullfile(folder,'native_report.json'),report);
end
