function pointer=saveCheckpoint(folder,state,authority,cfg,history,holdStart,sequence)
%SAVECHECKPOINT Accepted W=R-a and a; derived operators are rebuilt.
names={'rho','scale','step','normalizedTime','runMetadata','mass0','rhoRange0'};
payload=struct();for k=1:numel(names),payload.(names{k})=state.(names{k});end
checkpoint=struct('kind','ipm_c_centered_density_checkpoint_v4','sequence',sequence, ...
    'stateRepresentation','density_minus_origin','state',payload,'grid',state.ops.polarGrid,'rescaling',state.ops.rescaling, ...
    'authority',authority,'cfg',cfg,'history',history,'holdStart',holdStart);
% A resume/final/error publication can repeat an accepted sequence. Always
% write the UNPOINTED slot so a crash cannot invalidate the published hash.
cursor=fullfile(folder,'latest.json');slotIndex=0;generation=0;
if isfile(cursor)
    [~,previous]=ipm_c.output.loadCheckpoint(folder);
    slotIndex=double(strcmp(previous.path,fullfile(folder,'rolling_0.mat')));
    if isfield(previous,'generation'),generation=previous.generation+1;else,generation=1;end
end
slot=fullfile(folder,sprintf('rolling_%d.mat',slotIndex));
temporary=[slot,'.writing.mat'];save(temporary,'checkpoint','-v7.3');
loaded=load(temporary,'checkpoint');assert(isequaln(checkpoint,loaded.checkpoint),'ipm_c:SharedCheckpointReadback');
[ok,message]=movefile(temporary,slot,'f');assert(ok,'ipm_c:SharedCheckpointRename','%s',message);
pointer=struct('kind','ipm_c_centered_density_pointer_v4','sequence',sequence,'path',slot, ...
    'generation',generation,'sha256',ipm_c.output.sha256(slot), ...
    'tau',state.scale.canonicalTime,'physicalTime',state.scale.physicalTime, ...
    'physicalTimeLow',state.scale.physicalTimeLow,'densityOrigin',state.scale.densityOrigin);
ipm_c.output.atomicJson(fullfile(folder,'latest.json'),pointer);
end
