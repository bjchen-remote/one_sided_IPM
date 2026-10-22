function value=searchEvidence(action,varargin)
%IPM.REMESH.SEARCHEVIDENCE Pure bounded v4 search identity and prefix contract.
% No state sampling, geometry search, LU, flow, history rewriting or I/O.
value=[];id='ipm:AutonomousMeshSearchEvidence';
switch action
 case 'registration'
  maximumMembers=8;
  if ~isempty(varargin),maximumMembers=varargin{1};end
  validateattributes(maximumMembers,{'numeric'},{'scalar','integer','>=',1,'<=',192});
  rankingRule='quality_margin_then_schedule_index_then_y_index';
  if numel(varargin)>=2 && varargin{2}==5
   rankingRule='minimax_adjacent_ratio_then_quality_margin_then_schedule_and_y';
  end
  value=struct('version',1,'algorithm','primary_then_lower_octave_dyadic_v1', ...
   'activation','evolved_remesh_only','lowerExtensionOctaves',1,'refinementDepth',2, ...
   'maximumAxisTrialsPerMember',500,'maximumPairsPerMember',3,'maximumMembers',maximumMembers, ...
   'maximumTotalPairDescriptors',3*maximumMembers,'stageCounts',[70,260,170], ...
   'stopRule','first_completed_stage_with_qualified_pair', ...
   'rankingRule',rankingRule);
 case 'empty'
  p=varargin{1};initial=varargin{2};assert(any(p.version==[4,5])&&islogical(initial)&&isscalar(initial),id);
  phase='not_requested';if initial,phase='initial_original_planner';end
  value=struct('version',1,'rule',p.axisSearchPolicy,'phase',phase, ...
   'members',struct([]),'candidates',struct([]),'filteredCandidateIndices',zeros(1,0));
 case 'validate'
  validate_evidence(varargin{:});
 case 'attempts'
  validate_attempts(varargin{:});
 otherwise
  error(id,'Unknown search evidence operation.');
end
end

function validate_evidence(e,p,family,sourceLevel,initial,requested)
id='ipm:AutonomousMeshSearchEvidence';
if p.version==5
 expectedRule=ipm.remesh.searchEvidence('registration',size(p.nodeFamily.cellFactors,1),5);
else
 expectedRule=ipm.remesh.searchEvidence('registration');
end
need(any(p.version==[4,5])&&same(p.axisSearchPolicy,expectedRule),id);
fields(e,{'version','rule','phase','members','candidates','filteredCandidateIndices'},id);
need(same(e.version,1)&&same(e.rule,p.axisSearchPolicy)&&ischar(e.phase)&&isrow(e.phase),id);
need(islogical(initial)&&isscalar(initial)&&islogical(requested)&&isscalar(requested),id);
if initial||~requested
 expected=ipm.remesh.searchEvidence('empty',p,initial);need(same(e,expected),id);return
end
need(strcmp(e.phase,'evolved_requested')&&isstruct(e.members)&&isrow(e.members)&& ...
 numel(e.members)<=p.axisSearchPolicy.maximumMembers&&isstruct(e.candidates)&& ...
 numel(e.candidates)<=p.axisSearchPolicy.maximumTotalPairDescriptors,id);
need(scalar_integer(sourceLevel)&&sourceLevel>=1&&sourceLevel<=numel(family.members),id);
factors=vertcat(family.members.cellFactors);counts=vertcat(family.members.nodeCount);
ids=find(all(factors>=family.members(sourceLevel).cellFactors,2)& ...
 [family.members.resourceAdmitted]'&[family.members.qualityPassed]');
[~,order]=sortrows([prod(counts(ids,:),2),ids],[1,2]);ids=ids(order);
need(numel(e.members)==numel(ids)&&ids(1)==sourceLevel,id);
[schedule,~]=ipm.remesh.hierarchicalAxisSchedule(p.search);
need(size(schedule,1)==500,id);offset=0;
for j=1:numel(ids)
 m=e.members(j);member=family.members(ids(j));
 fields(m,{'levelId','nodeCount','phase','evaluatedAxisTrials','primaryPairCount', ...
  'primaryYIndices','completedRefinementLevel','returnedPairCount','status'},id);
 need(same(m.levelId,ids(j))&&same(m.nodeCount,member.nodeCount)&& ...
  ischar(m.phase)&&isrow(m.phase)&&ischar(m.status)&&isrow(m.status),id);
 for name={'evaluatedAxisTrials','primaryPairCount','completedRefinementLevel','returnedPairCount'}
  need(scalar_integer(m.(name{1}))&&m.(name{1})>=0,id);
 end
 need(m.primaryPairCount<=3&&m.returnedPairCount<=3&&isnumeric(m.primaryYIndices)&& ...
  isa(m.primaryYIndices,'double')&&isrow(m.primaryYIndices)&&isreal(m.primaryYIndices)&& ...
  all(isfinite(m.primaryYIndices))&&all(mod(m.primaryYIndices,1)==0)&& ...
  all(m.primaryYIndices>=1&m.primaryYIndices<=4)&&all(diff(m.primaryYIndices)>0),id);
 if strcmp(m.phase,'primary')
  need(m.evaluatedAxisTrials==70&&m.completedRefinementLevel==0&& ...
   m.returnedPairCount==m.primaryPairCount,id);
  if m.primaryPairCount>0
   need(strcmp(m.status,'primary_qualified'),id);
  else
   need(isempty(m.primaryYIndices)&&strcmp(m.status,'primary_y_family_exhausted_no_x_refinement'),id);
  end
 else
  need(strcmp(m.phase,'refined')&&m.primaryPairCount==0&&~isempty(m.primaryYIndices)&& ...
   any(m.completedRefinementLevel==[1,2]),id);
  if m.completedRefinementLevel==1
   need(m.evaluatedAxisTrials==330&&m.returnedPairCount>0,id);
  else
   need(m.evaluatedAxisTrials==500,id);
  end
  if m.returnedPairCount>0
   need(strcmp(m.status,'refined_axes_proposed_require_actual_transfer'),id);
  else
   need(strcmp(m.status,'registered_search_budget_exhausted_not_global_infeasibility'),id);
  end
 end
 for local=1:m.returnedPairCount
  offset=offset+1;need(offset<=numel(e.candidates),id);d=e.candidates(offset);
  fields(d,{'targetLevelId','localPairIndex','unchanged','searchPhase','refinementLevel','xScheduleIndex','yTrialIndex'},id);
  need(same(d.targetLevelId,m.levelId)&&same(d.localPairIndex,local)&& ...
   islogical(d.unchanged)&&isscalar(d.unchanged)&&same(d.searchPhase,m.phase)&& ...
   same(d.refinementLevel,m.completedRefinementLevel)&& ...
   scalar_integer(d.xScheduleIndex)&&scalar_integer(d.yTrialIndex),id);
  if d.unchanged
   need(m.levelId==sourceLevel&&local==1&&strcmp(m.phase,'primary')&& ...
    same(d.xScheduleIndex,0)&&same(d.yTrialIndex,0),id);
  else
   need(d.xScheduleIndex>=1&&d.xScheduleIndex<=m.evaluatedAxisTrials&& ...
    schedule(d.xScheduleIndex,1)==d.refinementLevel&& ...
    any(m.primaryYIndices==d.yTrialIndex),id);
  end
 end
end
need(offset==numel(e.candidates),id);
if isempty(e.candidates),expected=zeros(1,0);else,expected=find(~[e.candidates.unchanged]);end
need(same(e.filteredCandidateIndices,expected),id);
% Identity duplication cannot be disguised by different local rank labels.
if numel(e.candidates)>1
 keys=[[e.candidates.targetLevelId]',[e.candidates.xScheduleIndex]',[e.candidates.yTrialIndex]'];
 need(size(unique(keys,'rows'),1)==size(keys,1),id);
end
end

function validate_attempts(attempts,e,acceptedTarget)
id='ipm:AutonomousMeshSearchEvidence';
need(isstruct(attempts)&&isrow(attempts)&&~isempty(attempts)&& ...
 numel(attempts)<=numel(e.filteredCandidateIndices),id);
accepted=~isempty(acceptedTarget);
for k=1:numel(attempts)
 a=attempts(k);need(isstruct(a)&&isscalar(a)&& ...
  all(isfield(a,{'candidateIndex','passed','targetLevelId','localPairIndex','searchDescriptor'})),id);
 d=e.candidates(e.filteredCandidateIndices(k));
 need(same(a.candidateIndex,k)&&same(a.targetLevelId,d.targetLevelId)&& ...
  same(a.localPairIndex,d.localPairIndex)&&same(a.searchDescriptor,d)&& ...
  islogical(a.passed)&&isscalar(a.passed)&&a.passed==(accepted&&k==numel(attempts)),id);
end
if accepted,need(same(attempts(end).targetLevelId,acceptedTarget),id);end
end
function fields(v,names,id)
need(isstruct(v)&&isscalar(v)&&isempty(setxor(fieldnames(v),names)),id);
end
function yes=scalar_integer(v)
yes=isa(v,'double')&&isreal(v)&&isscalar(v)&&isfinite(v)&&v==fix(v);
end
function yes=same(a,b)
yes=strcmp(class(a),class(b))&&isequal(size(a),size(b));if ~yes,return;end
if isstruct(a)
 names=fieldnames(a);if ~isempty(setxor(names,fieldnames(b))),yes=false;return;end
 for j=1:numel(a)
  for k=1:numel(names)
   if ~same(a(j).(names{k}),b(j).(names{k})),yes=false;return;end
  end
 end
else
 yes=isequaln(a,b);
end
end
function need(ok,id)
assert(ok,id,'Invalid v4 bounded search identity, stage, ordering, or actual attempt prefix.');
end
