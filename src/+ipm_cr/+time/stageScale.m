function [s,clockInfo]=stageScale(previous,z,physicalIncrement)
%STAGESCALE One scale representation for RHS evaluation and accepted commit.
% The RK physical-time increment is supplied before it is added to the large
% clock. Recovering it from z(3)-previous.physicalTime would lose small steps.
s=previous;s.logC_l=z(1);s.logC_omega=z(2);s.X_shift=z(4);
s.canonicalTime=z(5);s.densityOrigin=z(6);
high=previous.physicalTime;low=previous.physicalTimeLow;
tail=zeros(1,2);
if isfield(previous,'physicalTimeTail'),tail=previous.physicalTimeTail;end
before=[high,low,reshape(tail,1,[])];
if physicalIncrement==0
    % A zero stage offset preserves the incoming clock representation.
    clockInfo=struct('kind','physical_clock_four_words_v1','legacyCompatible',true, ...
        'strictlyIncreased',false,'wordCount',4,'discardedAbsolute',0, ...
        'discardAccounting','Zero increment preserves every incoming clock word.');
else
    [s.physicalTime,s.physicalTimeLow,s.physicalTimeTail,clockInfo]= ...
        ipm_cr.time.clock(high,low,tail,physicalIncrement);
end
afterTail=zeros(1,2);
if isfield(s,'physicalTimeTail'),afterTail=s.physicalTimeTail;end
clockInfo.increment=physicalIncrement;
clockInfo.beforeWords=before;
clockInfo.afterWords=[s.physicalTime,s.physicalTimeLow,reshape(afterTail,1,[])];
end
