function [schedule,registration]=hierarchicalAxisSchedule(search)
%IPM.REMESH.HIERARCHICALAXISSCHEDULE Finite, state-independent integer refinement.
% First preserve the registered coarse enumeration. Then extend each lower
% endpoint by one octave, and halve parameter intervals at two depths.
f=search.fineCells;r=search.roundingCells;a=search.coreFineCellFractions;
validateattributes(f,{'numeric'},{'row','real','finite','integer','>=',8});
validateattributes(r,{'numeric'},{'row','real','finite','integer','>=',2});
validateattributes(a,{'numeric'},{'row','real','finite','>',0,'<',1});
assert(numel(f)*numel(r)*numel(a)<=500&&numel(unique(a))==numel(a));
assert(isrow(f)&&isrow(r)&&numel(f)>=2&&numel(r)>=2);
df=diff(f);dr=diff(r);assert(df(1)>0&&dr(1)>0&&all(df==df(1))&&all(dr==dr(1)));
assert(mod(df(1),4)==0&&mod(dr(1),4)==0&&mod(f(1),2)==0&&mod(r(1),2)==0);
minimumFine=max(8,f(1)/2);minimumRounding=max(2,r(1)/2);
assert(mod(f(end)-minimumFine,df(1)/2)==0&&mod(r(end)-minimumRounding,dr(1)/2)==0);
refinementMaximumRounding=r(end);
if isequal(r,[8,16,24,32,40])
    % The opt-in balanced-density primary grid replaces 48 by 8. Extend
    % only its refined range one half-step beyond 40 so stage one still
    % contains the registered 260 unique supplementary trials.
    refinementMaximumRounding=r(end)+dr(1)/2;
end
registration=struct('version',1,'lowerExtensionOctaves',1,'refinementDepth',2, ...
    'maximumAxisTrials',500,'maximumPairCandidates',3,'constructorMinimumFine',8, ...
    'constructorMinimumRounding',2,'fineBounds',[minimumFine,f(end)], ...
    'roundingBounds',[minimumRounding,refinementMaximumRounding],'fineBaseStep',df(1), ...
    'roundingBaseStep',dr(1),'fractions',a,'primaryEnumerationUnchanged',true, ...
    'stageOrder','complete primary; complete first refinement; budgeted second refinement', ...
    'refinedOrder','fine+2*(rounding+2), fine, rounding, original fraction order', ...
    'stopAfterFirstStageWithQualifiedPair',true,'stepOrClockInput',false, ...
    'globalFeasibilityClaim',false);
schedule=zeros(0,6);seen=zeros(0,3);
for level=0:2
    if level==0,ff=f;rr=r;
    else,ff=minimumFine:df(1)/2^level:f(end);rr=minimumRounding:dr(1)/2^level:refinementMaximumRounding;end
    rows=zeros(numel(ff)*numel(rr)*numel(a),6);n=0;
    for fine=ff
        for rounding=rr
            for j=1:numel(a)
                n=n+1;rows(n,:)=[level,fine,rounding,a(j),fine+2*(rounding+2),j];
            end
        end
    end
    if level>0
        rows=rows(~ismember(rows(:,2:4),seen,'rows'),:);
        rows=sortrows(rows,[5,2,3,6]);
    end
    remaining=registration.maximumAxisTrials-size(schedule,1);
    rows=rows(1:min(size(rows,1),remaining),:);
    schedule=[schedule;rows];seen=schedule(:,2:4); %#ok<AGROW>
    if size(schedule,1)==registration.maximumAxisTrials,break;end
end
assert(size(unique(schedule(:,2:4),'rows'),1)==size(schedule,1));
end
