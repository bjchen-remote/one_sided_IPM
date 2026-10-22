function report = validateController(metadata,config,current)
%IPM.REMESH.VALIDATECONTROLLER Validate persisted autonomous mesh evidence.
% CURRENT projects state/checkpoint/result into step/canonicalTime/
% physicalTime/normalizedTime/remeshCount/coreCells/safety/x/y/history.
% baseX/baseY may additionally be supplied by state/checkpoint callers.
% Version two also requires levelId/nodeCount and snapshots(rho/x/y cells).
% This is consistency validation, not replay of historical field transfers.
id = 'ipm:AutonomousMeshController';
require(isstruct(metadata) && isscalar(metadata),id,'Run metadata must be scalar.');
hasPolicy = isfield(config,'remesh') && isfield(config.remesh,'autonomousMesh');
hasMemory = isfield(metadata,'autonomousMesh');
report = struct('checked',false,'enabled',false,'transactions',0,'windowRecords',0, ...
    'historicalFieldTransfersReplayed',false,'inputModified',false);
if ~hasPolicy
    require(~hasMemory,id,'A controller ledger requires an enabled matching policy.');
    return;
end
policy = config.remesh.autonomousMesh;
schema = ipm.config.schema();choices = struct();
for k = 1:numel(schema.domainNames)
    group = config.(schema.domainNames{k});names = fieldnames(group);
    for j = 1:numel(names),choices.(names{j}) = group.(names{j});end
end
normalizedPolicy = ipm.config.autonomousMeshPolicy(policy,choices);
require(config.schemaVersion == 4 && same_value(policy,normalizedPolicy),id, ...
    'The stored policy must be a complete canonical registered policy.');
report.enabled = policy.enabled;
if ~policy.enabled
    require(~hasMemory,id,'A controller ledger requires an enabled matching policy.');
    return;
end
require(hasMemory,id,'An enabled policy requires its original controller ledger.');
memory = metadata.autonomousMesh;
fields(memory,{'version','policy','transactions','cumulativeAbsolutePeakJump', ...
    'window','initialization','lastDecision','lastFailure'},id,'controller');
require(same_value(memory.version,policy.version) && same_value(memory.policy,policy),id, ...
    'Controller version and policy must exactly match the configured policy.');
require(isstruct(memory.lastDecision) && isscalar(memory.lastDecision) && ...
    isstruct(memory.lastFailure) && isscalar(memory.lastFailure),id, ...
    'Decision/failure diagnostics must be scalar structures.');
fields(current,{'step','canonicalTime','physicalTime','normalizedTime','remeshCount', ...
    'coreCells','safety','x','y','history'},id,'current projection');
require(integer(current.step) && integer(current.remeshCount) && ...
    number(current.canonicalTime) && current.canonicalTime >= 0 && ...
    number(current.physicalTime) && current.physicalTime >= 0 && ...
    number(current.normalizedTime) && current.normalizedTime == current.canonicalTime && ...
    pair(current.coreCells) && number(current.safety) && current.safety >= 0,id, ...
    'Current projected counts, clocks, core or safety are invalid.');
versionTwo = any(policy.version == [2,3,4,5]);
if versionTwo
    fields(current,{'levelId','nodeCount','snapshots'},id,'version-two current projection');
    require(same_value(current.nodeCount,[numel(current.x),numel(current.y)]),id, ...
        'Projected node counts must exactly match the actual axes.');
    axis_pair(current.x,current.y,current.nodeCount(1),current.nodeCount(2),id);
else
    axis_pair(current.x,current.y,config.grid.nx,config.grid.ny,id);
end
init = memory.initialization;
fields(init,{'sourceBaseX','sourceBaseY','selectedBaseX','selectedBaseY', ...
    'selectedCandidateIndex','audit','attempts','analyticDatumResampled', ...
    'originalStep','originalCanonicalTime','originalPhysicalTime','originalRemeshCount', ...
    'newPhysicalEpochCreated','operatorGridRepresentation'},id,'initialization');
require(integer(init.originalStep) && init.originalStep == 0 && ...
    number(init.originalCanonicalTime) && init.originalCanonicalTime == 0 && ...
    number(init.originalPhysicalTime) && init.originalPhysicalTime == 0 && ...
    integer(init.originalRemeshCount) && init.originalRemeshCount == 0 && ...
    yes(init.analyticDatumResampled) && no(init.newPhysicalEpochCreated),id, ...
    'Initialization must resample the configured datum at original step/time/remesh zero.');
require(ischar(init.operatorGridRepresentation) && isrow(init.operatorGridRepresentation) && ...
    strcmp(init.operatorGridRepresentation,'custom_axes'),id, ...
    'Enabled controllers must retain the custom-axis operator construction used from initialization.');
if isfield(metadata,'caseMetadata') && isstruct(metadata.caseMetadata)
    require(~isfield(metadata.caseMetadata,'latePhysicalBoxBranch'),id, ...
        'A late fresh physical epoch cannot be relabelled original zero-time initialization.');
end
axis_pair(init.sourceBaseX,init.sourceBaseY,config.grid.nx,config.grid.ny,id);
axis_pair(init.selectedBaseX,init.selectedBaseY,config.grid.nx,config.grid.ny,id);
require(isequal(init.sourceBaseX([1,end]),current.x([1,end])) && ...
    isequal(init.sourceBaseY([1,end]),current.y([1,end])) && ...
    isequal(init.selectedBaseX([1,end]),current.x([1,end])) && ...
    isequal(init.selectedBaseY([1,end]),current.y([1,end])) && ...
    any(init.selectedBaseX == config.scaling.transportAnchorX) && ...
    any(init.selectedBaseX == -config.scaling.transportAnchorX),id, ...
    'Initial source/selected axes must preserve the domain, node budget and exact anchors.');
family = struct();
if versionTwo
    family = reference_family(memory,init,policy,current,config.scaling.transportAnchorX,id);
end
if isfield(current,'baseX') || isfield(current,'baseY')
    baseX=init.selectedBaseX;baseY=init.selectedBaseY;
    if versionTwo
        baseX=family.members(current.levelId).baseX;
        baseY=family.members(current.levelId).baseY;
    end
    sameBase=isfield(current,'baseX') && isfield(current,'baseY');
    if sameBase && versionTwo
        sameBase=same_value(current.baseX,baseX) && same_value(current.baseY,baseY);
    elseif sameBase
        sameBase=isequal(current.baseX,baseX) && isequal(current.baseY,baseY);
    end
    require(sameBase,id, ...
        'The native base axes must equal their registered immutable reference exactly.');
end
if isfield(config.grid,'customX') && ~isempty(config.grid.customX)
    require(isequal(config.grid.customX,init.sourceBaseX) && ...
        isequal(config.grid.customY,init.sourceBaseY),id, ...
        'Initial source axes must match the frozen explicit initial grid.');
end
initial_audit(init.audit,policy,id);
require(integer(init.selectedCandidateIndex) && init.selectedCandidateIndex >= 1 && ...
    init.selectedCandidateIndex <= policy.search.maximumPairCandidates && ...
    isstruct(init.attempts) && numel(init.attempts) == init.selectedCandidateIndex,id, ...
    'Initial selected candidate must identify its complete bounded attempt sequence.');
for k = 1:numel(init.attempts)
    a = init.attempts(k);fields(a,{'candidateIndex','passed','audit','errorIdentifier','errorMessage'},id,'initial attempt');
    require(integer(a.candidateIndex) && a.candidateIndex == k && ...
        islogical(a.passed) && isscalar(a.passed) && a.passed == (k == numel(init.attempts)),id, ...
        'Only the final selected initial attempt may be accepted.');
end
require(same_value(init.attempts(end).audit,init.audit) && ...
    isempty(init.attempts(end).errorIdentifier) && isempty(init.attempts(end).errorMessage),id, ...
    'The initial success audit must equal its selected successful attempt.');

initial_observation(init,config,id);

transactions = memory.transactions;n = numel(transactions);
require(isstruct(transactions) && n == current.remeshCount && n <= config.remesh.maxRemeshes,id, ...
    'Every native remesh must have exactly one committed transaction audit.');
jumps = zeros(1,n);steps = zeros(n,1);times = steps;physical = steps;
levelIds = ones(n+1,1);growthCount=0;
for k = 1:n
    a = transactions(k);
    fields(a,{'version','kind','passed','reasons','sourceStep','sourceCanonicalTime','sourcePhysicalTime', ...
        'sourceRemeshCount','coreCells','leftFrontCells','relativePeakJump', ...
        'priorCumulativeAbsolutePeakJump','cumulativeAbsolutePeakJump','massRelativeDefect', ...
        'relativeRangeViolation','xQuality','yQuality','clocksAndReferencesPreserved','sameBoxAndNodeCount'},id,'transaction');
    require(same_value(a.version,policy.version) && strcmp(a.kind,'autonomous_actual_candidate_audit') && ...
        yes(a.passed) && iscell(a.reasons) && isempty(a.reasons) && ...
        yes(a.clocksAndReferencesPreserved),id, ...
        'Only passed actual candidate audits with preserved native invariants may be committed.');
    if versionTwo
        [levelIds(k+1),grew] = level_transaction(a,family,levelIds(k),policy,id);
        growthCount=growthCount+grew;
        require(growthCount<=policy.nodeFamily.maximumAcceptedGrowthTransitions,id, ...
            'The committed growth transitions exceed their registered budget.');
        transaction_decision(a,family,policy,current.history,id);
    else
        require(yes(a.sameBoxAndNodeCount),id, ...
            'Version-one transactions must preserve the box and node count.');
    end
    require(integer(a.sourceStep) && a.sourceStep >= 1 && ...
        integer(a.sourceRemeshCount) && a.sourceRemeshCount == k-1 && ...
        number(a.sourceCanonicalTime) && a.sourceCanonicalTime > 0 && ...
        number(a.sourcePhysicalTime) && a.sourcePhysicalTime >= 0,id, ...
        'Transaction source counts and clocks must identify their native accepted step.');
    require(pair(a.coreCells) && all(a.coreCells >= policy.transactionMinimumCoreCells) && ...
        number(a.leftFrontCells) && a.leftFrontCells >= policy.minimumFrontCells && ...
        number(a.relativePeakJump) && a.relativePeakJump >= 0 && a.relativePeakJump <= policy.maximumSinglePeakJump && ...
        number(a.massRelativeDefect) && a.massRelativeDefect >= 0 && a.massRelativeDefect <= policy.maximumMassRelativeDefect && ...
        number(a.relativeRangeViolation) && a.relativeRangeViolation >= 0 && a.relativeRangeViolation <= policy.maximumRelativeRangeViolation,id, ...
        'A committed transaction violates the original core/front/peak/mass/range gates.');
    quality(a.xQuality,policy.qualityLimits,id);quality(a.yQuality,policy.qualityLimits,id);
    jumps(k) = a.relativePeakJump;
    prior = sum(jumps(1:k-1));
    require(number(a.priorCumulativeAbsolutePeakJump) && a.priorCumulativeAbsolutePeakJump == prior && ...
        number(a.cumulativeAbsolutePeakJump) && a.cumulativeAbsolutePeakJump == prior+a.relativePeakJump && ...
        a.cumulativeAbsolutePeakJump <= policy.maximumCumulativeAbsolutePeakJump && ...
        sum(jumps(1:k)) <= policy.maximumCumulativeAbsolutePeakJump,id, ...
        'Each transaction must retain its exact prior sum and actual added peak jump.');
    steps(k) = a.sourceStep;times(k) = a.sourceCanonicalTime;physical(k) = a.sourcePhysicalTime;
end
require(all(diff(steps)>0) && all(diff(times)>0) && all(diff(physical)>=0) && ...
    all(steps<=current.step) && all(times<=current.canonicalTime) && all(physical<=current.physicalTime),id, ...
    'Transaction chronology cannot repeat, move backwards or exceed the current accepted state.');
require(number(memory.cumulativeAbsolutePeakJump) && memory.cumulativeAbsolutePeakJump == sum(jumps) && ...
    memory.cumulativeAbsolutePeakJump <= policy.maximumCumulativeAbsolutePeakJump,id, ...
    'The controller cumulative budget must be the exact sum of actual committed audit jumps.');
if versionTwo
    require(same_value(levelIds(end),current.levelId),id, ...
        'The terminal level must equal the committed transaction prefix.');
end

w = memory.window;fields(w,{'time','step','remeshCount','coreCells','safety'},id,'window');
m = numel(w.time);
require(m >= 1 && column(w.time,m) && column(w.step,m) && column(w.remeshCount,m) && ...
    column(w.safety,m) && matrix(w.coreCells,m,2) && ...
    all(diff(w.time)>0) && all(diff(w.step)>0) && ...
    all(w.step>=0 & mod(w.step,1)==0) && all(w.remeshCount==current.remeshCount) && ...
    all(w.time>=0 & w.time>=current.canonicalTime-policy.trendWindow) && ...
    all(w.coreCells(:)>=policy.historyCoreFloor) && all(w.safety>=0 & w.safety<policy.maximumSafety),id, ...
    'The retained telemetry must be monotone and within the current mesh epoch/time window and hard floors.');
require(same_value(w.time(end),current.canonicalTime) && same_value(w.step(end),current.step) && ...
    same_value(w.coreCells(end,:),current.coreCells) && same_value(w.safety(end),current.safety),id, ...
    'The terminal window observation must exactly match the current native projection.');
if n>0
    require(all(w.time>=times(end)) && all(w.step>=steps(end)),id, ...
        'Current-epoch observations cannot precede its accepted transaction.');
end
history(current.history,current,w,steps,times,physical,id);
if versionTwo
    level_observations(current,w,family,levelIds,config.output.storeSnapshots,id);
    last_decision(memory.lastDecision,current,family,levelIds,transactions,policy,id);
end
if any(policy.version==[4,5])
    failure_search(memory.lastFailure,policy,family,current,levelIds,transactions,config.time,id);
end
report.checked = true;report.transactions = n;report.windowRecords = m;
end

function family = reference_family(memory,init,policy,current,anchor,id)
fields(memory,{'referenceFamily','currentLevelId'},id,'version-two controller');
try
    expected=ipm.remesh.referenceAxisFamily(init.selectedBaseX,init.selectedBaseY,anchor,policy);
catch exception
    error(id,'The registered reference family cannot be reconstructed: %s',exception.message);
end
family=memory.referenceFamily;
require(same_value(family,expected),id, ...
    'The persisted family must exactly match its original selected-base generator.');
require(level_id(memory.currentLevelId,numel(family.members)) && ...
    same_value(current.levelId,memory.currentLevelId),id,'Current level identifiers must match exactly.');
member=family.members(current.levelId);
require(yes(member.resourceAdmitted) && yes(member.qualityPassed) && ...
    same_value(current.nodeCount,member.nodeCount) && ...
    same_value(member.index,current.levelId),id,'The current level is not an admitted registered node budget.');
require(yes(family.members(1).resourceAdmitted) && yes(family.members(1).qualityPassed) && ...
    same_value(family.members(1).baseX,init.selectedBaseX) && ...
    same_value(family.members(1).baseY,init.selectedBaseY),id, ...
    'Level one must be the unchanged admitted initial selected reference.');
end

function [target,grew] = level_transaction(a,family,source,policy,id)
fields(a,{'sourceLevelId','targetLevelId','sourceNodeCount','targetNodeCount', ...
    'sourceCellFactors','targetCellFactors','sameBox','registeredNodeCounts', ...
    'matchingReferenceMembers','admittedLevelTransition'},id,'version-two transaction');
require(same_value(a.sourceLevelId,source) && level_id(a.targetLevelId,numel(family.members)),id, ...
    'Transaction levels must follow the actual committed prefix.');
old=family.members(source);target=a.targetLevelId;new=family.members(target);
require(same_value(a.sourceNodeCount,old.nodeCount) && same_value(a.targetNodeCount,new.nodeCount) && ...
    same_value(a.sourceCellFactors,old.cellFactors) && same_value(a.targetCellFactors,new.cellFactors),id, ...
    'Transaction counts and factors must exactly identify their reference members.');
require(yes(a.sameBox) && yes(a.registeredNodeCounts) && yes(a.matchingReferenceMembers) && ...
    yes(a.admittedLevelTransition) && yes(new.resourceAdmitted) && yes(new.qualityPassed) && ...
    all(new.cellFactors>=old.cellFactors) && prod(new.nodeCount)<=policy.nodeFamily.maximumTotalNodes && ...
    same_value(a.sameBoxAndNodeCount,isequal(old.nodeCount,new.nodeCount)),id, ...
    'A level transition must preserve the box and obey its registered monotone resource boundary.');
grew=target~=source;
end

function transaction_decision(a,family,policy,h,id)
fields(a,{'controllerDecision','candidateIndex','attemptSummary'},id,'transaction decision evidence');
d=a.controllerDecision;
decision_observation(d,policy,h,id);
names={'sourceStep','sourceCanonicalTime','sourcePhysicalTime','sourceRemeshCount','sourceLevelId','sourceNodeCount'};
for k=1:numel(names)
    require(same_value(d.(names{k}),a.(names{k})),id,'Decision source must equal the actual pre-transaction source.');
    require(isa(d.(names{k}),'double'),id,'Decision source fields must retain their native double type.');
end
require(yes(d.requested) && no(d.initial),id,'Committed transaction decisions must request an evolved remesh.');
require(isempty(d.stopReason) && same_value(d.predictedCoreCells, ...
    d.coreCells.*exp(-policy.maximumReviewInterval*d.decayEstimate)),id, ...
    'An accepted decision must retain the exact registered forecast and no stop reason.');
if any(d.coreCells<policy.regridCoreTrigger)
    require(strcmp(d.reason,'current_core_trigger'),id,'The actual core trigger has priority.');
else
    require(d.trendAvailable && any(d.predictedCoreCells<policy.predictedCoreBuffer & ...
        d.coreCells<policy.targetCoreCells) && strcmp(d.reason,'forecast_core_trigger'),id, ...
        'An accepted decision without a current-core trigger requires the actual forecast trigger.');
end
source=family.members(a.sourceLevelId);eligible=0;
for k=1:numel(family.members)
    member=family.members(k);
    eligible=eligible+(member.resourceAdmitted && member.qualityPassed && all(member.cellFactors>=source.cellFactors));
end
require(isa(a.candidateIndex,'double') && integer(a.candidateIndex) && a.candidateIndex>=1 && ...
    a.candidateIndex<=eligible*policy.search.maximumPairCandidates && ...
    isstruct(a.attemptSummary) && isrow(a.attemptSummary) && numel(a.attemptSummary)==a.candidateIndex,id, ...
    'The accepted candidate must identify its bounded complete attempt sequence.');
for k=1:numel(a.attemptSummary)
    row=a.attemptSummary(k);fields(row,{'candidateIndex','passed','errorIdentifier','auditReasons'},id,'attempt summary');
    require(same_value(row.candidateIndex,k) && same_value(row.passed,k==a.candidateIndex) && ...
        ischar(row.errorIdentifier) && (isrow(row.errorIdentifier)||isempty(row.errorIdentifier)) && ...
        iscell(row.auditReasons) && (isrow(row.auditReasons)||isempty(row.auditReasons)) && ...
        all(cellfun(@(s)ischar(s) && isrow(s),row.auditReasons)),id,'Invalid ordered candidate attempt summary.');
end
require(isempty(row.errorIdentifier) && isempty(row.auditReasons),id, ...
    'An accepted candidate must have no error or audit rejection reasons.');
if any(policy.version == [3,4,5])
    % V3 records each attempted level and its original local pair position.
    % Only the source level can omit local pair one (a filtered valid keep).
    factors=vertcat(family.members.cellFactors);counts=vertcat(family.members.nodeCount);
    ids=find(all(factors>=source.cellFactors,2) & [family.members.resourceAdmitted]' & ...
        [family.members.qualityPassed]');
    [~,order]=sortrows([prod(counts(ids,:),2),ids],[1,2]);ids=ids(order);
    lastRank=0;lastLocal=0;perLevel=zeros(numel(ids),1);
    for k=1:numel(a.attemptSummary)
        attempt=a.attemptSummary(k);
        fields(attempt,{'targetLevelId','localPairIndex'},id,'v3 per-level attempt');
        require(level_id(attempt.targetLevelId,numel(family.members)) && ...
            isa(attempt.localPairIndex,'double') && integer(attempt.localPairIndex) && ...
            attempt.localPairIndex>=1 && attempt.localPairIndex<=policy.search.maximumPairCandidates,id, ...
            'Each v3 attempt must retain a valid double target and local pair index.');
        rank=find(ids==attempt.targetLevelId);
        require(isscalar(rank) && rank>=lastRank,id,'V3 attempts must follow the admitted node-product member order.');
        if rank==lastRank
            require(attempt.localPairIndex==lastLocal+1,id,'A member cannot repeat or skip an attempted local pair.');
        else
            allowed=attempt.localPairIndex==1 || ...
                (attempt.targetLevelId==a.sourceLevelId && attempt.localPairIndex==2);
            require(allowed,id,'A new member starts at pair one, except the filtered source keep.');
        end
        perLevel(rank)=perLevel(rank)+1;
        require(perLevel(rank)<=policy.search.maximumPairCandidates,id,'A member exceeded its registered three-pair budget.');
        lastRank=rank;lastLocal=attempt.localPairIndex;
    end
    require(same_value(a.attemptSummary(end).targetLevelId,a.targetLevelId),id, ...
        'The final successful attempt must identify the actually committed target member.');
end
if any(policy.version==[4,5])
    fields(d,{'axisSearchEvidence'},id,'v4 committed search evidence');
    ipm.remesh.searchEvidence('validate',d.axisSearchEvidence,policy,family,a.sourceLevelId,false,true);
    ipm.remesh.searchEvidence('attempts',a.attemptSummary,d.axisSearchEvidence,a.targetLevelId);
end

end

function decision_observation(d,policy,h,id,unacceptedSource)
if nargin<5,unacceptedSource=false;end
fields(d,{'sourceStep','sourceCanonicalTime','sourcePhysicalTime','sourceNormalizedTime', ...
    'sourceRemeshCount','sourceLevelId','sourceNodeCount','requested','initial','coreCells', ...
    'predictedCoreCells','decayEstimate','reason','stopReason','trendAvailable','decayEvidence'},id,'controller decision');
names={'sourceStep','sourceCanonicalTime','sourcePhysicalTime','sourceNormalizedTime','sourceRemeshCount','sourceLevelId'};
for j=1:numel(names),require(isa(d.(names{j}),'double') && number(d.(names{j})),id,'Decision source clocks and indices must be finite doubles.');end
require(integer(d.sourceStep) && integer(d.sourceRemeshCount) && d.sourceCanonicalTime>=0 && d.sourcePhysicalTime>=0,id,'Invalid decision source counts or clocks.');
require(same_value(d.sourceNormalizedTime,d.sourceCanonicalTime) && ...
    islogical(d.requested) && isscalar(d.requested) && islogical(d.initial) && isscalar(d.initial) && ...
    isa(d.coreCells,'double') && pair(d.coreCells) && ...
    isa(d.predictedCoreCells,'double') && matrix(d.predictedCoreCells,1,2) && all(d.predictedCoreCells>=0) && ...
    isa(d.decayEstimate,'double') && matrix(d.decayEstimate,1,2) && all(d.decayEstimate>=0) && ...
    ischar(d.reason) && isrow(d.reason) && ~isempty(d.reason) && ...
    ischar(d.stopReason) && (isrow(d.stopReason)||isempty(d.stopReason)) && ...
    islogical(d.trendAvailable) && isscalar(d.trendAvailable),id,'Invalid transaction decision values or types.');
e=d.decayEvidence;
require(isstruct(e) && isequal(size(e),[1,2]),id,'Two ordered directional decay evidence records are required.');
for k=1:2
    fields(e(k),{'axis','fitDecay','maximumSecant','secantFromStep','secantToStep','secantFromTime','secantToTime'},id,'decay evidence');
    names={'fitDecay','maximumSecant','secantFromStep','secantToStep','secantFromTime','secantToTime'};
    for j=1:numel(names),require(isa(e(k).(names{j}),'double') && number(e(k).(names{j})),id,'Decay evidence scalars must be finite doubles.');end
    require(same_value(e(k).axis,k) && integer(e(k).secantFromStep) && integer(e(k).secantToStep) && ...
        e(k).secantFromStep<=e(k).secantToStep && e(k).secantToStep<=d.sourceStep && ...
        e(k).secantFromTime>=max(0,d.sourceCanonicalTime-policy.trendWindow) && ...
        e(k).secantFromTime<=e(k).secantToTime && e(k).secantToTime<=d.sourceCanonicalTime && ...
        same_value(d.decayEstimate(k),max([0,e(k).fitDecay,e(k).maximumSecant])),id, ...
        'Directional decay evidence is outside its source window or disagrees with the stored maximum.');
    if d.trendAvailable
        require(e(k).secantFromStep<e(k).secantToStep && e(k).secantFromTime<e(k).secantToTime,id, ...
            'An observed secant must span two distinct accepted observations.');
    else
        require(same_value([e(k).fitDecay,e(k).maximumSecant],[0,0]) && ...
            same_value([e(k).secantFromStep,e(k).secantToStep],[d.sourceStep,d.sourceStep]) && ...
            same_value([e(k).secantFromTime,e(k).secantToTime],[d.sourceCanonicalTime,d.sourceCanonicalTime]),id, ...
            'Unavailable trends must retain the explicit zero/current-source evidence.');
    end
    for side={'From','To'}
        step=e(k).(['secant',side{1},'Step']);time=e(k).(['secant',side{1},'Time']);
        matched=find(h.common.acceptedStep==step);
        if unacceptedSource && step==d.sourceStep,matched=[];end
        require(isempty(matched) || all(h.common.canonicalTau(matched)==time),id, ...
            'Decay evidence clocks must agree with any matching persisted history observation.');
    end
end
end

function last_decision(d,c,family,levels,transactions,policy,id)
decision_observation(d,policy,c.history,id);
if any(policy.version==[4,5])
    fields(d,{'axisSearchEvidence'},id,'v4 current search evidence');
    ipm.remesh.searchEvidence('validate',d.axisSearchEvidence,policy,family,d.sourceLevelId,d.initial,d.requested);
end
require(same_value(d.sourceStep,c.step) && same_value(d.sourceCanonicalTime,c.canonicalTime) && ...
    same_value(d.sourcePhysicalTime,c.physicalTime) && d.sourceRemeshCount<=c.remeshCount && ...
    same_value(d.sourceLevelId,levels(d.sourceRemeshCount+1)),id, ...
    'The last accepted decision must retain the current clocks and its correct pre-transaction level.');
member=family.members(d.sourceLevelId);
require(same_value(d.sourceNodeCount,member.nodeCount),id,'Last decision counts must match its source reference member.');
if ~isempty(transactions) && transactions(end).sourceStep==c.step
    require(same_value(d,transactions(end).controllerDecision),id, ...
        'Immediately after commit, the last decision must exactly retain the actual transaction source decision.');
else
    require(same_value(d.sourceRemeshCount,c.remeshCount) && same_value(d.sourceLevelId,c.levelId) && ...
        same_value(d.sourceNodeCount,c.nodeCount),id, ...
        'A later accepted step must record its current mesh epoch and reference member.');
end
if d.initial
    require(d.sourceStep==0 && d.sourceCanonicalTime==0 && d.sourcePhysicalTime==0 && d.sourceRemeshCount==0 && ...
        same_value(d.sourceLevelId,1) && yes(d.requested),id,'Initial decisions must describe the original zero-time source.');
end
matched=find(c.history.common.acceptedStep==d.sourceStep);
require(isempty(matched) || (all(c.history.common.canonicalTau(matched)==d.sourceCanonicalTime) && ...
    all(c.history.common.physicalTime(matched)==d.sourcePhysicalTime)),id,'Last decision must match any recorded source clocks.');
end

function failure_search(f,p,family,c,levels,transactions,time,id)
h=c.history;
% Failed candidate prefixes belong to their attempted source, which can be a
% rejected future step. Do not rewrite them to the rolled-back accepted clock.
if isempty(fieldnames(f)),return;end
if isfield(f,'controllerDecision')
    fields(f,{'attempts','kind'},id,'v4 transaction failure');
    require(strcmp(f.kind,'candidate_transaction_failure'),id,'Unexpected search failure kind.');
    d=f.controllerDecision;decision_observation(d,p,h,id);
    failure_source(d,c,family,levels,transactions,time,false,id);
    if isfield(f,'step'),require(same_value(f.step,d.sourceStep),id,'Failed attempt step disagrees with its decision.');end
    if isfield(f,'canonicalTime'),require(same_value(f.canonicalTime,d.sourceCanonicalTime),id,'Failed attempt clock disagrees with its decision.');end
    fields(d,{'axisSearchEvidence'},id,'v4 failed search evidence');
    require(no(d.initial)&&yes(d.requested),id,'A failed evolved transaction retains its requested plan.');
    ipm.remesh.searchEvidence('validate',d.axisSearchEvidence,p,family,d.sourceLevelId,false,true);
    ipm.remesh.searchEvidence('attempts',f.attempts,d.axisSearchEvidence,[]);
elseif isfield(f,'axisSearchEvidence')
    decision_observation(f,p,h,id,true);
    failure_source(f,c,family,levels,transactions,time,true,id);
    ipm.remesh.searchEvidence('validate',f.axisSearchEvidence,p,family,f.sourceLevelId,f.initial,f.requested);
else
    % Planning exceptions have no completed search to certify. They cannot
    % masquerade as a candidate attempt sequence or accepted plan.
    fields(f,{'kind','identifier','message'},id,'v4 planning exception');
    require(strcmp(f.kind,'planning_exception')&&~isfield(f,'attempts'),id, ...
        'An incomplete planning exception cannot claim a qualified search prefix.');
    fields(f,{'attemptedStep','sourceCanonicalTime','sourcePhysicalTime','sourceNormalizedTime'},id,'v4 planning exception clocks');
    failure_clock(f.attemptedStep,f.sourceCanonicalTime,f.sourcePhysicalTime,f.sourceNormalizedTime,c,time,true,id);
end
end

function failure_source(d,c,family,levels,transactions,time,mayBePending,id)
% Failed observations may survive a later resume. Bind them to the recorded
% historical mesh epoch; only the direct rollback path may be the next step.
require(no(d.initial)&&isa(d.sourceRemeshCount,'double')&&integer(d.sourceRemeshCount)&& ...
    d.sourceRemeshCount>=0&&d.sourceRemeshCount<=c.remeshCount,id,'Failure source epoch is outside the committed ledger.');
epoch=d.sourceRemeshCount;
require(same_value(d.sourceLevelId,levels(epoch+1))&& ...
    same_value(d.sourceNodeCount,family.members(d.sourceLevelId).nodeCount),id, ...
    'Failure source level and exact double node-count pair must match its immutable reference member.');
failure_clock(d.sourceStep,d.sourceCanonicalTime,d.sourcePhysicalTime,d.sourceNormalizedTime,c,time,mayBePending,id);
if d.sourceStep==c.step+1
    require(same_value(epoch,c.remeshCount)&&same_value(d.sourceLevelId,c.levelId)&& ...
        same_value(d.sourceNodeCount,c.nodeCount),id,'A pending rollback observation uses the current accepted grid epoch.');
end
if epoch>0
    previous=transactions(epoch);
    require(d.sourceStep>=previous.sourceStep&&d.sourceCanonicalTime>=previous.sourceCanonicalTime&& ...
        d.sourcePhysicalTime>=previous.sourcePhysicalTime,id,'Failure observation precedes the start of its mesh epoch.');
end
if epoch<numel(transactions)
    next=transactions(epoch+1);
    require(d.sourceStep<=next.sourceStep,id,'Failure source step exceeds the end of its mesh epoch.');
    if ~mayBePending
        require(d.sourceCanonicalTime<=next.sourceCanonicalTime&&d.sourcePhysicalTime<=next.sourcePhysicalTime,id, ...
            'An accepted failure source exceeds the end clock of its mesh epoch.');
    end
end
if mayBePending
    predecessor=find(c.history.common.acceptedStep==d.sourceStep-1);
    require(isempty(predecessor)||all(c.history.mesh.remeshCount(predecessor)==epoch),id, ...
        'A stored accepted predecessor must use the failed attempt source mesh epoch.');
end
end

function failure_clock(step,t,physical,normalized,c,time,mayBePending,id)
require(isa(step,'double')&&integer(step)&&step>=1&&step<=c.step+double(mayBePending)&& ...
    isa(t,'double')&&number(t)&&isa(physical,'double')&&number(physical)&& ...
    same_value(normalized,t)&&t>0&&physical>=0,id,'Failure attempt counts or clocks are invalid or beyond the next attempted step.');
if mayBePending
    % An old rejected attempt is not the later accepted observation bearing
    % the same step number. Terminal horizons can shorten on legal resume.
    require(t<=c.canonicalTime+time.maxDt,id,'Unaccepted attempt clock is beyond a conservative immutable step-size bound.');
    predecessor=find(c.history.common.acceptedStep==step-1);
    if ~isempty(predecessor)
        previousTime=c.history.common.canonicalTau(predecessor);previousPhysical=c.history.common.physicalTime(predecessor);
        require(all(t>previousTime & t<=previousTime+time.maxDt)&&all(physical>=previousPhysical),id, ...
            'Unaccepted attempt clock is outside the immutable step-size range of its stored predecessor.');
    end
    if step==c.step+1
        require(t>c.canonicalTime&&physical>=c.physicalTime,id,'A pending next attempt must follow the accepted endpoint.');
    end
    return
end
if step==c.step+1
    require(mayBePending&&t>c.canonicalTime&&t<=c.canonicalTime+time.maxDt&& ...
        physical>=c.physicalTime,id,'A rolled-back next step must retain a causal clock within the immutable canonical step-size limit.');
elseif step==c.step
    require(same_value(t,c.canonicalTime)&&same_value(physical,c.physicalTime),id,'A current accepted failure source must retain the exact current clocks.');
else
    require(t<c.canonicalTime&&physical<=c.physicalTime,id,'Historical failure clocks cannot exceed the accepted endpoint.');
end
h=c.history.common;matched=find(h.acceptedStep==step);
require(isempty(matched)||(all(h.canonicalTau(matched)==t)&&all(h.physicalTime(matched)==physical)),id, ...
    'Failure clocks must agree with every matching stored accepted history observation.');
end

function level_observations(c,w,family,levels,storeSnapshots,id)
fields(w,{'nodeCount','levelId'},id,'version-two window');m=numel(w.time);
require(same_value(w.nodeCount,repmat(c.nodeCount,m,1)) && ...
    same_value(w.levelId,repmat(c.levelId,m,1)),id,'The current-epoch window must retain actual registered dimensions and level.');
h=c.history.common;fields(h,{'nodeCountX','nodeCountY','meshLevelId'},id,'version-two history');
n=numel(h.acceptedStep);
require(isa(h.nodeCountX,'double') && column(h.nodeCountX,n) && ...
    isa(h.nodeCountY,'double') && column(h.nodeCountY,n) && ...
    isa(h.meshLevelId,'double') && column(h.meshLevelId,n),id,'History node counts and levels must be double columns.');
for k=1:n
    level=levels(c.history.mesh.remeshCount(k)+1);member=family.members(level);
    require(same_value(h.meshLevelId(k),level) && ...
        same_value([h.nodeCountX(k),h.nodeCountY(k)],member.nodeCount),id, ...
        'Every history record must match the node family at its actual transaction prefix.');
end
s=c.snapshots;fields(s,{'rho','x','y'},id,'version-two snapshots');
count=numel(s.rho);
require(iscell(s.rho) && iscell(s.x) && iscell(s.y) && numel(s.x)==count && numel(s.y)==count && ...
    ((storeSnapshots && count==n)||(~storeSnapshots && count==0)),id, ...
    'Snapshot counts must obey the frozen output policy and recorded history.');
for k=1:count
    member=family.members(h.meshLevelId(k));
    axis_pair(s.x{k},s.y{k},member.nodeCount(1),member.nodeCount(2),id);
    require(same_value(s.x{k}([1,end]),family.rootX([1,end])) && ...
        same_value(s.y{k}([1,end]),family.rootY([1,end])) && ...
        any(s.x{k}==family.anchor) && any(s.x{k}==-family.anchor) && ...
        matrix(s.rho{k},member.nodeCount(2),member.nodeCount(1)),id, ...
        'Each snapshot must retain its own actual registered dimensions, box and exact anchors.');
end
end

function valid = level_id(value,count)
valid=isa(value,'double') && integer(value) && value>=1 && value<=count;
end

function history(h,c,w,steps,times,physical,id)
fields(h,{'common','mesh'},id,'history');
fields(h.common,{'acceptedStep','canonicalTau','physicalTime','t'},id,'common history');
fields(h.mesh,{'remeshCount','coreGridPoints','verticalCoreGridPoints','safetyFactor'},id,'mesh history');
hs = h.common.acceptedStep(:);ht = h.common.canonicalTau(:);hp = h.common.physicalTime(:);hn = h.common.t(:);
he = h.mesh.remeshCount(:);hc = [h.mesh.coreGridPoints(:),h.mesh.verticalCoreGridPoints(:)];hz = h.mesh.safetyFactor(:);n = numel(hs);
require(n>=1 && column(hs,n) && column(ht,n) && column(hp,n) && column(hn,n) && ...
    column(he,n) && matrix(hc,n,2) && column(hz,n) && all(diff(hs)>0) && all(diff(ht)>0) && ...
    all(diff(hp)>=0) && all(mod(hs,1)==0) && all(mod(he,1)==0) && all(diff(he)>=0),id, ...
    'Native controller history must have aligned monotone accepted steps, clocks and remesh counts.');
require(hs(1)==0 && ht(1)==0 && hp(1)==0 && hn(1)==0 && he(1)==0 && ...
    isequal(ht,hn) && hs(end)==c.step && ht(end)==c.canonicalTime && hp(end)==c.physicalTime && ...
    he(end)==c.remeshCount && isequal(hc(end,:),c.coreCells) && hz(end)==c.safety,id, ...
    'Original zero-time history and terminal native projection must agree exactly.');
count = 0;
for k = 1:n
    while count<numel(steps) && steps(count+1)<=hs(k),count=count+1;end
    require(he(k)==count,id,'History remesh counts must equal the actual transaction prefix count.');
    if count>0
        require(times(count)<=ht(k) && physical(count)<=hp(k),id,'A committed transaction is ahead of its history record.');
        if steps(count)==hs(k)
            require(times(count)==ht(k) && physical(count)==hp(k),id,'Matching transaction/history steps must have identical clocks.');
        end
    end
    if count<numel(steps)
        require(times(count+1)>ht(k) && physical(count+1)>=hp(k),id,'A future transaction precedes an earlier history record.');
    end
end
[~,iw,ih] = intersect(w.step,hs);
require(isequal(w.time(iw),ht(ih)) && isequal(w.coreCells(iw,:),hc(ih,:)) && ...
    isequal(w.safety(iw),hz(ih)),id,'Shared window/history accepted steps must retain identical observations.');
end

function initial_audit(a,p,id)
fields(a,{'kind','passed','reasons','coreCells','leftFrontCells','qualityPassed','exactZeroTime', ...
    'xQuality','yQuality','massMigrationClaim','initialDensityExactlyResampled'},id,'initial audit');
require(strcmp(a.kind,'analytic_zero_time_axis_selection_audit') && yes(a.passed) && ...
    iscell(a.reasons) && isempty(a.reasons) && pair(a.coreCells) && ...
    all(a.coreCells>=p.transactionMinimumCoreCells) && number(a.leftFrontCells) && ...
    a.leftFrontCells>=p.minimumFrontCells && yes(a.qualityPassed) && yes(a.exactZeroTime) && ...
    no(a.massMigrationClaim) && yes(a.initialDensityExactlyResampled),id,'The initial accepted analytic selection audit is invalid.');
quality(a.xQuality,p.qualityLimits,id);quality(a.yQuality,p.qualityLimits,id);
end
function quality(q,p,id)
names = {'maximumAdjacentCellRatio','maximumLogSpacingCurvature','minimumStencilRcond', ...
    'minimumQuadratureWeightRatio','minimumQuadratureWeightToControlWidthRatio', ...
    'maximumQuadratureWeightToControlWidthRatio','quadratureWeightsStrictlyPositive'};
fields(q,names,id,'axis quality');
for k = 1:6,require(number(q.(names{k})) && q.(names{k})>=0,id,'Axis quality scalars must be finite and nonnegative.');end
require(yes(q.quadratureWeightsStrictlyPositive) && q.maximumAdjacentCellRatio<=p.maxAdjacentCellRatio && ...
    q.maximumLogSpacingCurvature<=p.maxLogSpacingCurvature && q.minimumStencilRcond>=p.minStencilRcond && ...
    q.minimumQuadratureWeightRatio>=p.minQuadratureWeightRatio && ...
    q.minimumQuadratureWeightToControlWidthRatio>=p.minWeightToControlWidth && ...
    q.maximumQuadratureWeightToControlWidthRatio<=p.maxWeightToControlWidth,id,'A persisted accepted audit violates an original mesh-quality gate.');
end
function axis_pair(x,y,nx,ny,id)
require(isnumeric(x) && isreal(x) && isrow(x) && numel(x)==nx && all(isfinite(x)) && ...
    isnumeric(y) && isreal(y) && iscolumn(y) && numel(y)==ny && all(isfinite(y)) && ...
    all(diff(x)>0) && all(diff(y)>0) && isequal(x,-fliplr(x)) && y(1)==0,id,'Invalid paired symmetric axes or node budget.');
end
function fields(value,names,id,label)
require(isstruct(value) && isscalar(value) && all(isfield(value,names)),id,['Incomplete ',label,'.']);
end
function yesValue = number(x)
yesValue = isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x);
end
function yesValue = integer(x)
yesValue = number(x) && x>=0 && mod(x,1)==0;
end
function yesValue = pair(x)
yesValue = matrix(x,1,2) && all(x>0);
end
function yesValue = column(x,n)
yesValue = matrix(x,n,1);
end
function yesValue = matrix(x,m,n)
yesValue = isnumeric(x) && isreal(x) && isequal(size(x),[m,n]) && all(isfinite(x),'all');
end
function value = yes(x)
value = islogical(x) && isscalar(x) && x;
end
function value = no(x)
value = islogical(x) && isscalar(x) && ~x;
end
function require(condition,id,message)
if ~condition,error(id,'%s',message);end
end
function equal = same_value(a,b)
equal = strcmp(class(a),class(b)) && isequal(size(a),size(b));
if ~equal,return;end
if isstruct(a)
    names = fieldnames(a);equal = isequal(sort(names),sort(fieldnames(b)));
    if ~equal,return;end
    for j = 1:numel(a)
        for k = 1:numel(names)
            if ~same_value(a(j).(names{k}),b(j).(names{k})),equal=false;return;end
        end
    end
else
    equal = isequaln(a,b);
end
end

function initial_observation(init,config,id)
hasPolicy=isfield(config.remesh,'initialMeshObservationFallback');
hasEvidence=isfield(init,'observationFallback');
if ~hasPolicy
    require(~hasEvidence,id,'Initial observation evidence requires its explicit policy.');return;
end
schema=ipm.config.schema();choices=struct();
for d=1:numel(schema.domainNames)
    group=config.(schema.domainNames{d});names=fieldnames(group);
    for j=1:numel(names),choices.(names{j})=group.(names{j});end
end
rule=ipm.config.initialMeshObservationPolicy(config.remesh.initialMeshObservationFallback,choices);
require(same_value(rule,config.remesh.initialMeshObservationFallback),id,'Initial observation policy must be exactly canonical.');
if ~rule.enabled,require(~hasEvidence,id,'Disabled initial observation cannot have a reserved ledger.');return;end
require(hasEvidence,id,'Enabled initial observation must record whether it was used.');
e=init.observationFallback;
fields(e,{'version','policy','used','attemptCount','selectedPhase'},id,'initial observation');
require(same_value(e.version,1)&&same_value(e.policy,rule)&&islogical(e.used)&&isscalar(e.used), ...
    id,'Initial observation rule/version/usage must be canonical.');
if ~e.used
    expected=struct('version',1,'policy',rule,'used',false,'attemptCount',0,'selectedPhase','original');
    require(same_value(e,expected),id,'An unused observer must retain exactly the zero-attempt original phase.');return;
end
fields(e,{'status','originalAxes','originalAttempt','sourceNodeCount','sourceCoreCells', ...
    'originalFailureClass','observation','candidateAudits','admittedOriginalPairIndices', ...
    'selectedLocalCandidateIndex','selectedOriginalPairIndex','noSourceFieldInterpolation', ...
    'noNativeQualificationClaim','axisReport'},id,'used initial observer');
p=config.remesh.autonomousMesh;
require(same_value(e.attemptCount,1)&&strcmp(e.selectedPhase,'fixed_observation')&& ...
    strcmp(e.status,'native_initialization_accepted')&&yes(e.noSourceFieldInterpolation)&& ...
    yes(e.noNativeQualificationClaim)&&same_value(e.sourceNodeCount,[config.grid.nx,config.grid.ny])&& ...
    same_value(e.originalAxes.x,init.sourceBaseX)&&same_value(e.originalAxes.y,init.sourceBaseY), ...
    id,'Used initial observer must preserve one attempt and the genuine original solver axes/budget.');
a=e.originalAttempt;fields(a,{'plan','attempts'},id,'original initial failure');d=a.plan;
fields(d,{'initial','requested','sourceStep','sourceCanonicalTime','sourcePhysicalTime','sourceNormalizedTime', ...
    'sourceRemeshCount','sourceLevelId','sourceNodeCount','coreCells','candidates','stopReason','axisReport'},id,'original initial plan');
require(yes(d.initial)&&yes(d.requested)&&same_value(d.sourceStep,0)&&same_value(d.sourceCanonicalTime,0)&& ...
    same_value(d.sourcePhysicalTime,0)&&same_value(d.sourceNormalizedTime,0)&& ...
    same_value(d.sourceRemeshCount,0)&&same_value(d.sourceLevelId,1)&& ...
    same_value(d.sourceNodeCount,e.sourceNodeCount)&&same_value(d.coreCells,e.sourceCoreCells)&& ...
    ~isfield(d,'observationFallback')&&numel(d.candidates)<=p.search.maximumPairCandidates&& ...
    numel(a.attempts)<=p.search.maximumPairCandidates,id,'Original failure must precede any observation or evolution.');
capacity=strcmp(d.stopReason,'autonomous_mesh_axis_capacity')&&isempty(d.candidates)&&isempty(a.attempts);
if capacity
    require(strcmp(e.originalFailureClass,'finite_family_exhaustion'),id,'Wrong initial failure classification.');
else
    require(strcmp(e.originalFailureClass,'actual_analytic_core_or_front')&&isempty(d.stopReason)&& ...
        ~isempty(d.candidates)&&numel(a.attempts)==numel(d.candidates),id,'A fallback needs a complete rejected original stage.');
    for k=1:numel(a.attempts)
        z=a.attempts(k);fields(z,{'candidateIndex','passed','audit','errorIdentifier','errorMessage'},id,'original rejected attempt');
        q=z.audit;fields(q,{'kind','passed','reasons','qualityPassed','exactZeroTime','coreCells','leftFrontCells'},id,'original rejected analytic audit');
        require(same_value(z.candidateIndex,k)&&no(z.passed)&&no(q.passed)&&isempty(z.errorIdentifier)&& ...
            isempty(z.errorMessage)&&strcmp(q.kind,'analytic_zero_time_axis_selection_audit')&& ...
            same_value(q.reasons,{'initial_axis_or_analytic_resolution_gate'})&&yes(q.qualityPassed)&& ...
            yes(q.exactZeroTime)&&pair(q.coreCells)&&number(q.leftFrontCells)&& ...
            (any(q.coreCells<p.transactionMinimumCoreCells)||q.leftFrontCells<p.minimumFrontCells)&& ...
            any(d.candidates(k).x==1)&&any(d.candidates(k).x==-1),id,'Original failure was not solely analytic core/front resolution.');
    end
end
o=e.observation;
fields(o,{'x','y','nodeCount','xQuality','yQuality','features','datum','derivative', ...
    'canonicalTime','physicalTime','normalizedTime','solverBudgetChanged'},id,'fixed initial observation');
ox=observer_axis(config.grid.xlim(2),rule,id);oy=observer_axis(config.grid.ymax,rule,id)';ox=[-fliplr(ox(2:end)),ox];
require(same_value(o.x,ox)&&same_value(o.y,oy)&&same_value(o.nodeCount,[numel(ox),numel(oy)])&& ...
    prod(o.nodeCount)<=rule.maximumObservationNodes&&same_value(o.canonicalTime,0)&& ...
    same_value(o.physicalTime,0)&&same_value(o.normalizedTime,0)&&no(o.solverBudgetChanged)&& ...
    strcmp(o.datum,'original_k8_initialDensity')&&strcmp(o.derivative,'paired_seven_point_fdMatrix'), ...
    id,'Saved observation axes must exactly reproduce the fixed one-dimensional rule.');
qx=ipm.mesh.quality(ox,7,ipm.mesh.quadrature(ox));qy=ipm.mesh.quality(oy,7,ipm.mesh.quadrature(oy));
require(same_value(qx,o.xQuality)&&same_value(qy,o.yQuality),id,'Observation quality must match the real saved axes.');
quality(qx,p.qualityLimits,id);quality(qy,p.qualityLimits,id);
rows=e.candidateAudits;require(isstruct(rows)&&numel(rows)>=1&&numel(rows)<=p.search.maximumPairCandidates, ...
    id,'Observed candidate stage must be a bounded nonempty top-three set.');
admitted=[];
for k=1:numel(rows)
    z=rows(k);fields(z,{'originalPairIndex','x','y','features','xQuality','yQuality', ...
        'geometryPassed','familyPassed','familyFailure','resourceMask','passed','nativeFlowEvaluated'},id,'observed pure candidate');
    require(same_value(z.originalPairIndex,k)&&no(z.nativeFlowEvaluated)&& ...
        islogical(z.passed)&&isscalar(z.passed),id,'Pure candidate rows must preserve original pair order and cannot claim a flow.');
    if z.passed
        require(yes(z.geometryPassed)&&yes(z.familyPassed)&&isempty(fieldnames(z.familyFailure)),id,'Admitted pure candidate gates disagree.');
        axis_pair(z.x,z.y,config.grid.nx,config.grid.ny,id);
        require(same_value(z.x([1,end]),config.grid.xlim)&&z.y(end)==config.grid.ymax&& ...
            any(z.x==1)&&any(z.x==-1)&&pair(z.features.actualCoreCells)&& ...
            all(z.features.actualCoreCells>=p.transactionMinimumCoreCells)&& ...
            number(z.features.leftFrontCells)&&z.features.leftFrontCells>=p.minimumFrontCells,id,'Observed candidate core/front/anchor failed.');
        xq=ipm.mesh.quality(z.x,7,ipm.mesh.quadrature(z.x));yq=ipm.mesh.quality(z.y,7,ipm.mesh.quadrature(z.y));
        require(same_value(xq,z.xQuality)&&same_value(yq,z.yQuality),id,'Observed candidate quality does not match its axes.');
        quality(xq,p.qualityLimits,id);quality(yq,p.qualityLimits,id);
        family=ipm.remesh.referenceAxisFamily(z.x,z.y,1,p);
        require(all([family.members.qualityPassed])&&family.members(1).resourceAdmitted&& ...
            same_value(z.resourceMask,[family.members.resourceAdmitted]),id,'Observed candidate family failed.');
        admitted(end+1)=k; %#ok<AGROW>
    end
end
require(same_value(admitted,e.admittedOriginalPairIndices)&& ...
    same_value(e.selectedLocalCandidateIndex,init.selectedCandidateIndex)&& ...
    e.selectedLocalCandidateIndex>=1&&e.selectedLocalCandidateIndex<=numel(admitted)&& ...
    same_value(e.selectedOriginalPairIndex,admitted(e.selectedLocalCandidateIndex)),id,'Selected local candidate must map to its admitted original pair.');
z=rows(e.selectedOriginalPairIndex);
require(same_value(init.selectedBaseX,z.x)&&same_value(init.selectedBaseY,z.y),id,'Native selected base must equal the actual observed candidate.');
end

function axis=observer_axis(endpoint,p,id)
inner=p.innerHalfWidth;spacing=p.innerSpacing;ramp=p.tailLogSlopeRampCells;
span=endpoint-inner;maximumLogRatio=log(p.maximumTailRatio);
n=ceil(log1p((p.maximumTailRatio-1)*span/spacing)/maximumLogRatio)+ramp;
require(isfinite(n)&&n>0&&n<=p.maximumObservationNodes,id,'Invalid saved observation tail budget.');
indices=(0:n-1)';exponent=indices-ramp/2;mask=indices<ramp;exponent(mask)=indices(mask).^2/(2*ramp);
require(spacing*n<span&&spacing*sum(exp(maximumLogRatio*exponent))>=span,id,'Saved observation tail has no fixed bracket.');
rate=fzero(@(v)spacing*sum(exp(v*exponent))-span,[0,maximumLogRatio]);
tail=inner+cumsum(spacing*exp(rate*exponent));tail(end)=endpoint;axis=[linspace(0,inner,161),tail'];
end
