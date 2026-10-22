function report=test_dataset_io(folder,state)
%TEST_DATASET_IO Portable trajectory, native slices and mesh-aware selection.
% A real audited remesh supplies the second geometry. Step numbers and
% clocks below are explicitly manufactured I/O labels; no time evolution or
% long-time accuracy is claimed by this storage-contract test.
assert(startsWith(folder,filesep)&&~isfolder(folder));mkdir(folder);
if nargin<2
    cfg=ipm_c.config();cfg.grid.nr=385;cfg.grid.ntheta=128;
    cfg.runtime.threads=1;cfg.output.verbose=false;state=ipm_c.initialize(cfg);
end
features=ipm_c.diagnostics.sourceFeatures(state.flow.source,state.ops.polarGrid,[.1,.2,.5]);
assert(features.usableForCdf,'ipm_c:DatasetTestSource');
original=state;original.config.remesh.adaptiveRemesh=true;original.config.remesh.maxRemeshes=Inf;
original.ops.rescaling.referenceLevelSetAreas=features.areas/1.3;
original.ops.rescaling.referenceLevelSetRadialBounds=features.radialProjectionBounds;
original.ops.rescaling.areaRemeshObserver=struct('checks',0,'failures',0,'retryAt',0,'features',[]);
r=original.ops.polarGrid.r;z=linspace(0,1,numel(r))';proposal=r.*(1+1e-6*sin(pi*z).^2);
proposal([1,end])=r([1,end]);
fixture=struct('manufactured',true,'proposalAxis',proposal,'allowIdentity',false,'maximumRetries',0);
[migrated,transfer]=ipm_c.remesh.apply(original,fixture);
assert(transfer.applied&&transfer.geometryChanged&&isequaln(migrated.ops.polarGrid.r,proposal), ...
    'ipm_c:DatasetTestRealMeshTransfer');

first=stamp(original,0,10,1,.2,-.1,.4);
second=stamp(migrated,1,20,2,.4,-.2,.6);
history=[row(first,NaN);row(second,10)];
native=fullfile(folder,'native');mkdir(native);
f0=ipm_c.output.frame(native,first,history(1,:));
f1=ipm_c.output.frame(native,second,history);
raw=load(fullfile(native,'data.mat'),'data');catalog=raw.data;
assert(catalog.schemaVersion==3&&size(catalog.history,2)==20&& ...
    numel(catalog.frames)==2&&numel(dir(fullfile(native,'grids','mesh_*.mat')))==2&& ...
    ~isfield(catalog,'rho')&&~isfield(catalog,'X'),'ipm_c:DatasetTestSmallCatalog');
assert(~strcmp(catalog.frames(1).meshId,catalog.frames(2).meshId)&& ...
    all(~startsWith({catalog.frames.file},filesep))&& ...
    all(~startsWith({catalog.frames.meshFile},filesep)),'ipm_c:DatasetTestRelativeCatalog');
paths=cellfun(@(p)fullfile(native,p),{catalog.frames.file},'UniformOutput',false);
hashes=cellfun(@ipm_c.output.sha256,paths,'UniformOutput',false);

d=ipm_c.data(native);assert(d.selected.tau==20&&d.selected.frameIndex==2,'ipm_c:DatasetTestLatest');
names={'step','tau','t','H','Qinf','Finf','cl','cw','cr','rhoInf','poisson','gauge', ...
    'remesh','logC_l','logC_omega','X_shift','rhoXinf','rhoYinf','uInf','mass'};
for k=1:numel(names)
    assert(iscolumn(d.(names{k}))&&isequaln(d.(names{k}),history(:,k)), ...
        'ipm_c:DatasetTestNamedTrajectory','Incorrect trajectory column: %s',names{k});
end
assert(isequaln(d.C_l,exp(history(:,14)))&&isequaln(d.C_omega,exp(history(:,15)))&& ...
    isequaln(d.c,history(:,7:9))&&isequaln(d.residual_inf,d.Qinf),'ipm_c:DatasetTestDerivedTrajectory');
check_fields(d,second,f1);
old=ipm_c.data(d,1,'index');assert(old.selected.tau==10,'ipm_c:DatasetTestIndex');
check_fields(old,first,f0);
exact=ipm_c.data(old,20,'exact');assert(exact.selected.frameIndex==2,'ipm_c:DatasetTestExact');
nearest=ipm_c.data(exact,11,'nearest');
assert(nearest.selected.frameIndex==1&&nearest.selected.requested==11&& ...
    ~nearest.selected.timeInterpolationUsed&&isequaln(nearest.tau,history(:,2)), ...
    'ipm_c:DatasetTestNearest');
expect_error(@()ipm_c.data(d,15,'exact'),'ipm_c:DatasetTimeNotSaved');

% A moved run uses only relative frame and mesh references. Check the copy
% before changing either source or copy, then detect corrupted mesh bytes.
portable=fullfile(folder,'portable');[ok,message]=copyfile(native,portable);
assert(ok,'ipm_c:DatasetTestCopy','%s',message);
copy=ipm_c.data(portable);check_fields(copy,second,f1);
assert(strcmp(copy.meta.baseFolder,portable),'ipm_c:DatasetTestPortableBase');
meshPath=fullfile(portable,catalog.frames(2).meshFile);saved=load(meshPath,'grid');grid=saved.grid;
grid.r(3)=grid.r(3)+1e-8*(grid.r(4)-grid.r(3));save(meshPath,'grid','-v7.3');
expect_error(@()ipm_c.data(portable),'ipm_c:DatasetMeshMismatch');
copyfile(fullfile(native,catalog.frames(2).meshFile),meshPath,'f');
copy=ipm_c.data(portable);check_fields(copy,second,f1);

% A checkpoint-only publication advances every scalar vector while spatial
% selection remains at its most recent actual saved slice.
third=stamp(second,2,25,3,.5,-.3,.8);longHistory=[history;row(third,5)];
beforeFrames=numel(dir(fullfile(native,'data_frames','frame_*.mat')));
ipm_c.output.frame(native,third,longHistory,false);after=ipm_c.data(native);
assert(numel(after.frames)==2&&numel(after.tau)==3&&after.tau(end)==25&& ...
    after.selected.tau==20&&after.selected.step==1&& ...
    numel(dir(fullfile(native,'data_frames','frame_*.mat')))==beforeFrames&& ...
    numel(dir(fullfile(native,'grids','mesh_*.mat')))==2, ...
    'ipm_c:DatasetTestScalarOnlyPublication');
check_fields(after,second,f1);

% Emulate publication immediately after restoring an older durable prefix.
% Immutable future files survive; only the catalog drops their references.
ipm_c.output.frame(native,first,history(1,:),true);restored=ipm_c.data(native);
assert(numel(restored.frames)==1&&restored.selected.step==0&& ...
    isequaln(restored.tau,history(1,2))&& ...
    numel(dir(fullfile(native,'data_frames','frame_*.mat')))==beforeFrames, ...
    'ipm_c:DatasetTestRestorePrefix');
for k=1:numel(paths)
    assert(strcmp(ipm_c.output.sha256(paths{k}),hashes{k}), ...
        'ipm_c:DatasetTestImmutableBytes','Restore changed a previously saved native slice.');
end
check_fields(restored,first,f0);
report=struct('kind','polar_dataset_native_mesh_portability_and_restore_prefix_io_v3', ...
    'passed',true,'actualAuditedMeshTransfer',transfer,'spatialFramesTested',2, ...
    'scalarColumnsTested',20,'meshRoundTrip',true,'portableRelativeReferences',true, ...
    'latestIndexNearestExactSelection',true,'physicalScalingVerified',true, ...
    'checkpointOnlyScalarsVerified',true,'meshCorruptionDetected',true, ...
    'restorePrefixCatalogVerified',true,'immutableFutureBytesPreserved',true, ...
    'clocksAreManufacturedIoLabels',true,'timeEvolutionTested',false, ...
    'spatialTrajectoryQualified',false);
end

function state=stamp(state,step,tau,t,ll,lo,shift)
state.step=step;state.scale.canonicalTime=tau;state.scale.physicalTime=t;
state.scale.logC_l=ll;state.scale.logC_omega=lo;state.scale.X_shift=shift;
state.runMetadata.datasetIoManufactured=true;
end

function values=row(s,h)
f=s.flow;d=ipm_c.field.derivatives(s.rho,s.ops.polarGrid,s.config.rhsOptions);
values=[s.step,s.scale.canonicalTime,s.scale.physicalTime,h, ...
    max(abs(f.fullQ),[],'all'),max(abs(s.rhsCache.rhoRate),[],'all'), ...
    f.c_l,f.c_omega,f.c_r,max(abs(s.rho),[],'all'),f.poissonResidual, ...
    max(abs(f.gaugeResidual)),s.ops.remeshCount,s.scale.logC_l,s.scale.logC_omega,s.scale.X_shift, ...
    max(abs(f.source),[],'all'),max(abs(d.ry),[],'all'), ...
    max(hypot(f.u1,f.u2),[],'all'),sum(s.rho.*s.ops.integrationWeights,'all')];
end

function check_fields(d,s,f)
g=s.ops.polarGrid;
assert(isequaln(d.r,g.r(:))&&isequaln(d.theta,g.theta(:).')&& ...
    isequaln(d.X,g.X)&&isequaln(d.Y,g.Y)&&d.grid.originIsSinglePhysicalPoint, ...
    'ipm_c:DatasetTestNativeMesh');
names={'rho','rho_X','rho_Y','rho_tau','Q','psi','u1','u2','vx','vy'};
for k=1:numel(names)
    assert(isequaln(d.(names{k}),f.(names{k})),'ipm_c:DatasetTestNativeField', ...
        'The selected native field changed: %s',names{k});
end
assert(isequaln(d.rho_x,d.rho_X)&&isequaln(d.rho_y,d.rho_Y)&& ...
    isequaln(d.F,d.rho_tau)&&isequaln(d.rho_xt,d.Q)&& ...
    isequaln(d.vmag,hypot(d.u1,d.u2))&&isequaln(d.transport_speed,hypot(d.vx,d.vy)), ...
    'ipm_c:DatasetTestFieldAliases');
Cl=exp(s.scale.logC_l);Co=exp(s.scale.logC_omega);
assert(isequaln(d.physical_X,(g.X-s.scale.X_shift)/Cl)&& ...
    isequaln(d.physical_Y,g.Y/Cl)&&isequaln(d.rho_physical,s.rho/Co)&& ...
    isequaln(d.rho_x_physical,(Cl/Co)*s.flow.source)&& ...
    isequaln(d.u1_physical,s.flow.u1/Co)&&isequaln(d.u2_physical,s.flow.u2/Co)&& ...
    isequaln(d.psi_physical,s.flow.psi/(Cl*Co)), ...
    'ipm_c:DatasetTestPhysicalUnits');
end

function expect_error(action,identifier)
try
    action();error('ipm_c:DatasetTestExpectedFailure','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_c:DatasetTestWrongFailure', ...
        'Expected %s, got %s.',identifier,problem.identifier);
end
end
