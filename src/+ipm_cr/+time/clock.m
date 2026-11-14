classdef clock
    %CLOCK Compensated physical time and RK stage scales.
    methods (Static)
        function [high,low,tail,info]=add(high,low,tail,increment)
            %IPM_CLOCK_ADD Positive physical-clock increment with four binary64 words.
            % Experimental helper only. OLD physicalTime/physicalTimeLow imports use
            % tail=zeros(1,2). The physical clock is high+low+tail(1)+tail(2), with the
            % words retained separately. Summing them in binary64 loses the extension.
            %
            % An empty tail imports an old two-word clock. While that old arithmetic
            % increases, its high/low bits are retained exactly and tail stays zero.
            % At its first stagnation, error-free expansion arithmetic retains the
            % increment in the third/fourth words. The exact-sign expansion comparison is returned as a diagnostic;
            % numerical verification is performed independently by tests.
            if isempty(tail),tail=zeros(1,2);end
            tail=reshape(tail,1,2);old=[high,low,tail];

            % Preserve the pre-extension trajectory wherever its old clock can advance.
            [legacyHigh,legacyLow]=legacy_add(high,low,increment);
            legacyIncreased=isfinite(legacyHigh)&&isfinite(legacyLow)&& ...
                (legacyHigh>high||(legacyHigh==high&&legacyLow>low));
            if ~any(tail)&&legacyIncreased
                change=expansion([legacyHigh,legacyLow,-high,-low]);
                legacyIncreased=~isempty(change)&&change(end)>0;
            end
            if ~any(tail)&&legacyIncreased
                high=legacyHigh;low=legacyLow;
                if nargout>3
                    info=struct('kind','physical_clock_four_words_v1','legacyCompatible',true, ...
                        'strictlyIncreased',true,'wordCount',4,'discardedAbsolute',0, ...
                        'discardAccounting','Legacy-compatible path retains historical rounding unchanged.');
                end
                return;
            end

            e=expansion([old,increment]);
            discarded=0;
            if numel(e)>4
                discarded=sum(abs(e(1:end-4)));
                e=e(end-3:end);
            end
            words=zeros(1,4);words(1:numel(e))=fliplr(e);
            change=expansion([words,-old]);
            increased=~isempty(change)&&change(end)>0;
            high=words(1);low=words(2);tail=words(3:4);
            if nargout>3
                info=struct('kind','physical_clock_four_words_v1','legacyCompatible',false, ...
                    'strictlyIncreased',increased,'wordCount',4,'discardedAbsolute',discarded, ...
                    'discardAccounting','Absolute bound on omitted fifth and smaller expansion terms in this update.');
            end
        end

        function [s,clockInfo]=stage(previous,z,physicalIncrement)
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
                    ipm_cr.time.clock.add(high,low,tail,physicalIncrement);
            end
            afterTail=zeros(1,2);
            if isfield(s,'physicalTimeTail'),afterTail=s.physicalTimeTail;end
            clockInfo.increment=physicalIncrement;
            clockInfo.beforeWords=before;
            clockInfo.afterWords=[s.physicalTime,s.physicalTimeLow,reshape(afterTail,1,[])];
        end
    end
end

function [high,low]=legacy_add(high,low,increment)
% Identical operation order to the frozen solver's clock_add.
s=high+increment;b=s-high;e=(high-(s-b))+(increment-b);
t=low+e;high=s+t;low=(s-high)+t;
end

function e=expansion(values)
% Grow exact nonoverlapping expansions, then remove redundant components.
% Components are ordered from the smallest magnitude to the largest.
e=[];
for value=values
    q=value;out=zeros(1,numel(e)+1);n=0;
    for item=e
        [q,residual]=two_sum(q,item);
        if residual~=0,n=n+1;out(n)=residual;end
    end
    if q~=0||n==0,n=n+1;out(n)=q;end
    e=out(1:n);
end
if isempty(e),return;end
q=e(end);buffer=zeros(size(e));bottom=numel(e);
for j=numel(e)-1:-1:1
    [next,residual]=two_sum(q,e(j));
    if residual~=0
        buffer(bottom)=next;bottom=bottom-1;q=residual;
    else
        q=next;
    end
end
buffer(bottom)=q;q=buffer(bottom);out=zeros(size(e));n=0;
for j=bottom+1:numel(e)
    [q,residual]=two_sum(buffer(j),q);
    if residual~=0,n=n+1;out(n)=residual;end
end
if q~=0||n==0,n=n+1;out(n)=q;end
e=out(1:n);
end

function [sumValue,errorValue]=two_sum(a,b)
% Knuth's error-free TwoSum; no fast-math or reassociation is admissible.
sumValue=a+b;virtualB=sumValue-a;
errorValue=(a-(sumValue-virtualB))+(b-virtualB);
end
