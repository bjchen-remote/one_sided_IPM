function memory=controllerTelemetry(memory,state,policy)
%IPM.REMESH.CONTROLLERTELEMETRY Retain a restartable, finite-time mesh trend.
% Every datum comes from the accepted state; no rate or clock is modified.
% The time horizon does not yet impose a fixed bound on the sample count.
if isempty(memory)
    assert(state.step==0 && state.scale.canonicalTime==0 && state.scale.physicalTime==0, ...
        'ipm:AutonomousMeshMemory','New controller memory must start at the original zero-time state.');
    memory=struct('version',1,'policy',policy,'transactions',struct([]), ...
        'cumulativeAbsolutePeakJump',0,'window',struct('time',[],'step',[], ...
        'remeshCount',[],'coreCells',zeros(0,2),'safety',[]), ...
        'initialization',struct(),'lastDecision',struct(),'lastFailure',struct());
    if any(policy.version == [2,3,4,5])
        memory.version=policy.version;
        memory.referenceFamily=struct();memory.currentLevelId=1;
        memory.window.nodeCount=zeros(0,2);memory.window.levelId=[];
    end
end
assert(memory.version==policy.version && isequaln(memory.policy,policy) && ...
    numel(memory.transactions)==state.ops.remeshCount, ...
    'ipm:AutonomousMeshMemory','Policy and transaction history must match the current native state.');
if isempty(memory.transactions),total=0;else,total=sum([memory.transactions.relativePeakJump]);end
assert(total==memory.cumulativeAbsolutePeakJump && total<=policy.maximumCumulativeAbsolutePeakJump, ...
    'ipm:AutonomousMeshMemory','The actual transaction ledger must reconstruct the cumulative budget.');
w=memory.window;t=state.scale.canonicalTime;n=state.step;epoch=state.ops.remeshCount;
if ~isempty(w.time)
    assert(t>=w.time(end) && n>=w.step(end) && epoch>=w.remeshCount(end), ...
        'ipm:AutonomousMeshMemory','Accepted mesh telemetry cannot move backwards.');
    if t==w.time(end) || n==w.step(end)
        assert(t==w.time(end) && n==w.step(end), ...
            'ipm:AutonomousMeshMemory','A same-time transaction must preserve the accepted step.');
        names=fieldnames(w);for k=1:numel(names),w.(names{k})(end,:)=[];end
    end
end
w.time(end+1,1)=t;w.step(end+1,1)=n;w.remeshCount(end+1,1)=epoch;
w.coreCells(end+1,:)=[state.flow.coreGridPoints,state.flow.verticalCoreGridPoints];
w.safety(end+1,1)=state.flow.safetyFactor;
if any(policy.version == [2,3,4,5])
    actual=[numel(state.ops.x),numel(state.ops.y)];
    assert(isequal(actual,[state.ops.nx,state.ops.ny]), ...
        'ipm:AutonomousMeshMemory','Current node counts must match the actual operators.');
    if isempty(fieldnames(memory.referenceFamily))
        assert(state.step==0 && epoch==0 && memory.currentLevelId==1 && ...
            isequal(actual,[state.config.grid.nx,state.config.grid.ny]), ...
            'ipm:AutonomousMeshMemory','Only initial planning may precede family registration.');
    else
        member=memory.referenceFamily.members(memory.currentLevelId);
        assert(member.resourceAdmitted && member.qualityPassed && ...
            isequal(actual,member.nodeCount) && ...
            isequal(state.ops.baseX,member.baseX) && isequal(state.ops.baseY,member.baseY), ...
            'ipm:AutonomousMeshMemory','The active reference member must match native operators.');
    end
    w.nodeCount(end+1,:)=actual;w.levelId(end+1,1)=memory.currentLevelId;
end
assert(all(isfinite([w.time;w.step;w.remeshCount;w.coreCells(:);w.safety])) && ...
    all(w.coreCells(:)>0),'ipm:AutonomousMeshMemory','Controller observations must be finite.');
keep=w.remeshCount==epoch & w.time>=t-policy.trendWindow;
names=fieldnames(w);for k=1:numel(names),w.(names{k})=w.(names{k})(keep,:);end
memory.window=w;
end
