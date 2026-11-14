function varargout=checkpoint(action,folder,varargin)
%CHECKPOINT Accepted-only v2 checkpoint and mesh transactions.
% [cp,pointer]=checkpoint('load',folder) reads portable v2 records only.
% [pointer,cp]=checkpoint('save',folder,cp,previousPointer) audits and commits.
% checkpoint('recordMesh',folder,before,after,info) stages a mesh event.
% checkpoint('verifyMesh',folder,history) verifies referenced event identities.
% The writer must hold the run's owner lock. Stale pointer chains are rejected.
switch action
    case 'seed',[varargout{1:nargout}]=import_seed(folder,varargin{:});
    case 'load',[varargout{1:nargout}]=read_v2(folder);
    case 'save',[varargout{1:nargout}]=write_v2(folder,varargin{:});
    case 'recordMesh',varargout{1}=record_mesh(folder,varargin{:});
    case 'verifyMesh',verify_mesh(folder,varargin{:});
    otherwise,error('ipm_cr:CheckpointAction','Unknown checkpoint action %s.',action);
end
end

function [pointer,cp]=write_v2(folder,cp,previous)
timer=tic;cursor=fullfile(folder,'latest.json');
if isfile(cursor)
    durable=jsondecode(fileread(cursor));
    assert(~isempty(previous)&&isequaln(durable,previous),'ipm_cr:CheckpointStale', ...
        'The durable accepted pointer changed; reload it before another commit.');
    assert(cp.state.step>=previous.sequence&&cp.state.scale.canonicalTime>=previous.tau, ...
        'ipm_cr:CheckpointSequence','An accepted commit cannot move the trajectory backwards.');
else
    assert(isempty(previous),'ipm_cr:CheckpointStale','The previous accepted pointer is missing.');
end
assert(isequaln(cp.sourceInventory,ipm_cr.io.identity('source')), ...
    'ipm_cr:CheckpointSourceChanged','Source changed before an accepted commit.');
initial='data/initial.json';
if isfield(cp,'explicitOptions')&&isfield(cp.explicitOptions,'initial'),initial=cp.explicitOptions.initial;end
assert(isequaln(cp.authority,ipm_cr.io.identity('authority',[],initial)), ...
    'ipm_cr:CheckpointAuthorityChanged','Initial authority changed before an accepted commit.');
cp.kind='ipm_cr_accepted_checkpoint_v2';cp.schemaVersion=2;
cp.stateRepresentation='W_centered_density_eta_analytic_residual';
cp.acceptedOnly=true;cp.solverFactorCacheSaved=false;cp.sequence=cp.state.step;
cp.retainedStartCheckpoint=true;
cp.contextKey=cp.state.ctxKey;
cp.spatialAudit=ipm_cr.diagnostics.overlap(cp.state,cp.geometry);
entry=struct('step',cp.state.step,'tau',cp.state.scale.canonicalTime,'audit',cp.spatialAudit);
if isempty(cp.auditHistory),cp.auditHistory=entry;else,cp.auditHistory(end+1)=entry;end
slot=0;generation=0;
if ~isempty(previous)
    assert(strcmp(previous.kind,'ipm_cr_accepted_pointer_v2'),'ipm_cr:CheckpointVersion','Only v2 pointers may be advanced.');
    slot=double(strcmp(previous.path,'accepted_0.mat'));generation=previous.generation+1;
end
name=sprintf('accepted_%d.mat',slot);path=fullfile(folder,name);temporary=[path,'.writing.mat'];
checkpoint=cp; %#ok<NASGU>
save(temporary,'checkpoint','-v7.3');check=load(temporary,'checkpoint');
assert(isequaln(check.checkpoint,cp),'ipm_cr:CheckpointReadback', ...
    'Accepted checkpoint readback differs from its complete payload.');
[ok,message]=movefile(temporary,path,'f');assert(ok,'ipm_cr:CheckpointRename','%s',message);
pointer=struct('kind','ipm_cr_accepted_pointer_v2','schemaVersion',2,'path',name, ...
    'generation',generation,'sequence',cp.sequence,'tau',cp.state.scale.canonicalTime, ...
    'sha256',ipm_cr.io.identity('sha',path),'accumulatedWallSeconds',cp.accumulatedWallSeconds+toc(timer));
if isempty(previous)
    % Retain the actual start state when the two rolling slots advance.
    % Publish it before the first accepted cursor, so a durable run always
    % has both a hash-bound start and its latest accepted endpoint.
    initial=fullfile(folder,'initial.mat');pending=[initial,'.writing'];
    [ok,message]=copyfile(path,pending,'f');assert(ok,'ipm_cr:CheckpointInitial','%s',message);
    assert(strcmp(ipm_cr.io.identity('sha',pending),pointer.sha256),'ipm_cr:CheckpointInitial', ...
        'Initial checkpoint copy differs from the accepted start state.');
    [ok,message]=movefile(pending,initial,'f');assert(ok,'ipm_cr:CheckpointInitial','%s',message);
    start=struct('kind','ipm_cr_initial_pointer_v1','path','initial.mat', ...
        'tau',pointer.tau,'sequence',pointer.sequence,'sha256',pointer.sha256);
    ipm_cr.io.identity('json',fullfile(folder,'initial.json'),start);
end
ipm_cr.io.identity('json',fullfile(folder,'latest.json'),pointer);
end

function [cp,pointer]=read_v2(folder)
pointer=jsondecode(fileread(fullfile(folder,'latest.json')));
assert(isstruct(pointer)&&isscalar(pointer)&&all(isfield(pointer, ...
    {'kind','schemaVersion','path','generation','sequence','tau','sha256'}))&& ...
    strcmp(pointer.kind,'ipm_cr_accepted_pointer_v2')&&pointer.schemaVersion==2&& ...
    any(strcmp(pointer.path,{'accepted_0.mat','accepted_1.mat'})), ...
    'ipm_cr:CheckpointPointer','Invalid v2 accepted-state pointer.');
path=fullfile(folder,pointer.path);
assert(isfile(path)&&strcmp(ipm_cr.io.identity('sha',path),pointer.sha256), ...
    'ipm_cr:CheckpointIdentity','Accepted checkpoint SHA does not match the durable pointer.');
loaded=load(path,'checkpoint');assert(isfield(loaded,'checkpoint'),'ipm_cr:CheckpointPayload','Missing checkpoint payload.');
cp=loaded.checkpoint;
required={'kind','schemaVersion','stateRepresentation','acceptedOnly','state','geometry','wallOps', ...
    'explicitOptions','sourceInventory','authority','history','historyColumns','auditHistory','sequence','contextKey'};
assert(isstruct(cp)&&isscalar(cp)&&all(isfield(cp,required))&& ...
    strcmp(cp.kind,'ipm_cr_accepted_checkpoint_v2')&&cp.schemaVersion==2&&cp.acceptedOnly&& ...
    strcmp(cp.stateRepresentation,'W_centered_density_eta_analytic_residual')&& ...
    all(isfield(cp.state,{'W','rhs','scale','step','ctxKey'}))&& ...
    cp.sequence==pointer.sequence&&cp.state.step==pointer.sequence&& ...
    cp.state.scale.canonicalTime==pointer.tau&&isequaln(cp.state.ctxKey,cp.contextKey), ...
    'ipm_cr:CheckpointState','Accepted checkpoint schema or sequence is inconsistent.');
% Earlier v2 checkpoints remain readable without inventing a saved start.
% Every checkpoint written by this implementation declares and requires it.
if isfield(cp,'retainedStartCheckpoint')&&cp.retainedStartCheckpoint
    initialFile=fullfile(folder,'initial.json');
    assert(isfile(initialFile),'ipm_cr:CheckpointInitial','The retained starting checkpoint manifest is missing.');
    start=jsondecode(fileread(initialFile));
    assert(isstruct(start)&&isscalar(start)&&all(isfield(start,{'kind','path','tau','sequence','sha256'}))&& ...
        strcmp(start.kind,'ipm_cr_initial_pointer_v1')&&strcmp(start.path,'initial.mat')&& ...
        pointer.tau>=start.tau&&pointer.sequence>=start.sequence&& ...
        isfile(fullfile(folder,start.path))&&strcmp(ipm_cr.io.identity('sha',fullfile(folder,start.path)),start.sha256), ...
        'ipm_cr:CheckpointInitial','Retained starting checkpoint identity is missing or changed.');
end
end

function entry=record_mesh(folder,before,after,info)
% The accepted cursor owns this lineage. Orphan events after a crash cannot
% overwrite an earlier event, and are not mistaken for committed geometry.
directory=fullfile(folder,'mesh');if ~isfolder(directory),mkdir(directory);end
temporary=[tempname(directory),'.mat'];
event=struct('format','ipm_cr_remesh_event_v1','info',info, ...
    'before',before.geometryKey,'after',after.geometryKey);
save(temporary,'event','-v7.3');readback=load(temporary,'event');
assert(isequaln(event,readback.event),'ipm_cr:RemeshStorage','Mesh event readback differs.');
entry=struct('step',info.step,'tau',info.tau,'status',info.status, ...
    'reason',info.reason,'file','','sha256',ipm_cr.io.identity('sha',temporary));
[~,name,extension]=fileparts(temporary);entry.file=['mesh/',name,extension];
end

function verify_mesh(folder,history)
% Verify the committed mesh event file identities.
for k=1:numel(history)
    entry=history(k);
    assert(~isempty(regexp(entry.file,'^mesh/[A-Za-z0-9_]+\.mat$','once')), ...
        'ipm_cr:RemeshStorage','Mesh lineage must reference a local mesh event basename.');
    file=fullfile(folder,strrep(entry.file,'/',filesep));
    assert(isfile(file)&&strcmp(ipm_cr.io.identity('sha',file),entry.sha256), ...
        'ipm_cr:RemeshStorage','Missing or changed mesh lineage file: %s.',entry.file);
end
end

function [state,ctx]=import_seed(method,seed)
% Initialize a new run from an accepted state at its recorded time.
assert(isfile(seed.file)&&strcmp(ipm_cr.io.identity('sha',seed.file),seed.sha256), ...
    'ipm_cr:SeedIdentity','Seed file does not match its declared SHA256.');
loaded=load(seed.file,'checkpoint');
assert(strcmp(ipm_cr.io.identity('sha',seed.file),seed.sha256),'ipm_cr:SeedIdentity', ...
    'Seed changed during read.');
assert(isfield(loaded,'checkpoint'),'ipm_cr:Seed','A complete accepted checkpoint is required.');
parent=loaded.checkpoint;
assert(all(isfield(parent,{'state','geometry','wallOps','acceptedOnly'}))&&parent.acceptedOnly&& ...
    isfield(parent.state,'W')&&parent.state.scale.canonicalTime==seed.tau, ...
    'ipm_cr:Seed','Seed is not an accepted state at the declared actual tau.');
[~,cfg]=ipm_cr.config(method);[wall,~,~]=ipm_cr.model.reference(cfg);
assert(isequaln(wall,parent.wallOps),'ipm_cr:SeedAuthority','Seed original wall differs.');
assert(isfield(parent,'explicitOptions')&& ...
    isequaln(parent.explicitOptions.background,method.background), ...
    'ipm_cr:SeedBackground','Import cannot reinterpret the residual under a different background.');
assert(isfield(parent,'authority')&&isfield(parent.authority,'sha256')&& ...
    strcmp(parent.authority.sha256,cfg.authoritySha), ...
    'ipm_cr:SeedAuthority','Seed initial authority differs.');
g=ipm_cr.mesh.build(method.geometry);ctx=ipm_cr.context(g,wall,method);
[state,transfer]=ipm_cr.mesh.transfer(parent.state,parent.geometry,g,true);
if isfield(state,'remesh'),state=rmfield(state,'remesh');end
state.ctxId=ctx.id;state.ctxKey=ctx.key;
[state.rhs,state.W]=ipm_cr.rhs(state,ctx);
comparison=seed_comparison(parent.state,parent.geometry,state,g,method);
state.import=struct('kind','explicit_cross_method_or_grid_seed_v1', ...
    'file',seed.file,'sha256',seed.sha256,'tau',seed.tau,'parentStep',state.step, ...
    'transfer',transfer,'comparison',comparison, ...
    'sameMethodContinuationClaimed',false,'timeRelabeled',false, ...
    'originalGeometryKind',parent.geometry.kind,'initialGeometryKind',g.kind, ...
    'gaugeResidual',state.rhs.gaugeResidual,'poisson',state.rhs.poisson);
end

function result=seed_comparison(a,ga,b,gb,method)
[x,y]=meshgrid(linspace(gb.core.activeXBounds(1),gb.core.activeXBounds(2),17), ...
    linspace(0,gb.core.activeYBounds(2),17));
r=exp(linspace(log(gb.polar.innerActiveRadius),log(gb.polar.r(end)),17))';r(end)=gb.polar.r(end);
theta=linspace(0,pi,17);px=r*cos(theta);py=r*sin(theta);py(:,[1,end])=0;
X=[x(:);px(:)];Y=[y(:);py(:)];result=struct('queryCount',numel(X));
% Honor the same analytic-background representation as the public sampler.
% Interpolating complete W would instead measure coarse interpolation of H,
% even when the actual transferred residual is exactly unchanged.
av=ipm_cr.diagnostics.fields(a,ga,method,X,Y);bv=ipm_cr.diagnostics.fields(b,gb,method,X,Y);
for name={'W','q','F','Q'}
    key=name{1};aa=av.(key);bb=bv.(key);
    result.(key)=max(abs(aa-bb)./max(1,max(abs(aa),abs(bb))));
end
result.scope='saved parent and fresh import using analytic background plus sampled residual at common physical points; not continuum truth';
end
