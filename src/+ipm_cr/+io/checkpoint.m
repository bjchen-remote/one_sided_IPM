function varargout=checkpoint(action,folder,varargin)
%CHECKPOINT Accepted-only v2 transactions and explicit read-only v1 import.
% [cp,pointer]=checkpoint('load',folder) reads portable v2 records only.
% [pointer,cp]=checkpoint('save',folder,cp,previousPointer) audits and commits.
% The writer must hold the run's owner lock. Stale pointer chains are rejected.
% checkpoint('legacy',folder) returns an untouched historical v1 payload;
% it does not migrate it, rebind a context, or authorize v2 continuation.
% checkpoint('import',folder) reads the same v1 data with W/eta/rhs state
% names for diagnostics only; it never writes or authorizes continuation.
switch action
    case 'load',[varargout{1:nargout}]=read_v2(folder);
    case 'save',[varargout{1:nargout}]=write_v2(folder,varargin{:});
    case 'legacy',[varargout{1:nargout}]=read_legacy(folder);
    case 'import'
        [cp,pointer]=read_legacy(folder);
        s=cp.state;old={'rho','residual','flow'};names={'W','eta','rhs'};
        for k=1:numel(old)
            if isfield(s,old{k})
                assert(~isfield(s,names{k}),'ipm_cr:StateAlias','Duplicate state representations.');
                s.(names{k})=s.(old{k});s=rmfield(s,old{k});
            end
        end
        cp.state=s;varargout={cp,pointer};
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
assert(isequaln(cp.authority,ipm_cr.io.identity('authority')), ...
    'ipm_cr:CheckpointAuthorityChanged','Initial authority changed before an accepted commit.');
cp.kind='ipm_cr_accepted_checkpoint_v2';cp.schemaVersion=2;
cp.stateRepresentation='W_centered_density_eta_analytic_residual';
cp.acceptedOnly=true;cp.solverFactorCacheSaved=false;cp.sequence=cp.state.step;
cp.contextKey=cp.state.ctxKey;
cp.spatialAudit=ipm_cr.diag.overlap(cp.state,cp.geometry);
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
    'sha256',ipm_cr.io.sha(path),'accumulatedWallSeconds',cp.accumulatedWallSeconds+toc(timer));
ipm_cr.io.json(fullfile(folder,'latest.json'),pointer);
end

function [cp,pointer]=read_v2(folder)
pointer=jsondecode(fileread(fullfile(folder,'latest.json')));
if isfield(pointer,'kind')&&strcmp(pointer.kind,'ipm_cr_explicit_accepted_pointer_v1')
    error('ipm_cr:CheckpointLegacy','Use checkpoint(''legacy'',folder) for an explicit read-only v1 import.');
end
assert(isstruct(pointer)&&isscalar(pointer)&&all(isfield(pointer, ...
    {'kind','schemaVersion','path','generation','sequence','tau','sha256'}))&& ...
    strcmp(pointer.kind,'ipm_cr_accepted_pointer_v2')&&pointer.schemaVersion==2&& ...
    any(strcmp(pointer.path,{'accepted_0.mat','accepted_1.mat'})), ...
    'ipm_cr:CheckpointPointer','Invalid v2 accepted-state pointer.');
path=fullfile(folder,pointer.path);
assert(isfile(path)&&strcmp(ipm_cr.io.sha(path),pointer.sha256), ...
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
end

function [cp,pointer]=read_legacy(folder)
pointer=jsondecode(fileread(fullfile(folder,'latest.json')));
assert(strcmp(pointer.kind,'ipm_cr_explicit_accepted_pointer_v1')&& ...
    any(strcmp(pointer.path,{fullfile(folder,'accepted_0.mat'),fullfile(folder,'accepted_1.mat')}))&& ...
    strcmp(ipm_cr.io.sha(pointer.path),pointer.sha256),'ipm_cr:LegacyCheckpointIdentity');
loaded=load(pointer.path,'checkpoint');cp=loaded.checkpoint;
assert(strcmp(cp.kind,'ipm_cr_explicit_accepted_checkpoint_v1')&&cp.acceptedOnly&& ...
    cp.sequence==pointer.sequence&&cp.state.step==pointer.sequence&& ...
    cp.state.scale.canonicalTime==pointer.tau&&isequaln(cp.rho,cp.state.rho), ...
    'ipm_cr:LegacyCheckpointState','The original v1 state is inconsistent.');
end
