function [s,r,h,event]=event(initial,s,r,h,wallKnots,options,trialProbe,event)
    %EVENT Locate original-peak crossings using complete RK4 trials.
    % Only event searching changes H; all trial states still originate
    % at INITIAL and are accepted by the ordinary outer-step gate.
    event.initialProposedDt=h;event.acceptedDt=h;
    probes=0;guardCuts=0;
    start=ipm_cr.time.peak(initial.rhs);
    if isempty(wallKnots)||~start.available||~isfinite(start.x),return;end
    while true
        crossing=find_crossing(r,h,start.x);
        if ~crossing.found,break;end
        event.knotX=crossing.knot;event.direction=crossing.direction;
        xtol=64*eps(max(1,abs(crossing.knot)));
        endpoint=ipm_cr.time.peak(s.rhs);
        fStart=crossing.direction*(start.x-crossing.knot);
        fEnd=crossing.direction*(endpoint.x-crossing.knot);
        if abs(fEnd)<=xtol
            event.located=true;event.rootResidual=endpoint.x-crossing.knot;
            event.status='located_peak_cell_boundary';break;
        end
        if fEnd<0
            % An internal RK prediction crossed, but its true substep
            % endpoint did not bracket that crossing. A half step is a
            % guard, not a falsely reported located physical event.
            h=h/2;guardCuts=guardCuts+1;
            if h<options.minDt||initial.scale.canonicalTime+h<=initial.scale.canonicalTime
                error('ipm_cr:PeakEventMinimumStep', ...
                    'A peak-cell stage guard requires dt below minDt.');
            end
            if guardCuts>24
                error('ipm_cr:PeakEventUnresolved','Peak-cell stage guards did not resolve a crossing.');
            end
            [s,r]=probe(h);event.status='stage_crossing_guard';
            event.rootResidual=s.rhs.measurement.peak.x-crossing.knot;
            continue;
        end
        % The endpoint brackets the first observed knot in this
        % direction. Locate endpoint.x(h)=knot using complete RK4
        % substeps, with safeguarded secant/bisection interpolation.
        lo=0;hi=h;flo=fStart;fhi=fEnd;
        bestState=s;bestRecords=r;bestH=h;bestResidual=fEnd;
        located=false;slowProgress=0;
        for iteration=1:64
            if abs(bestResidual)<=xtol,located=true;break;end
            timeTolerance=16*eps(max(1,abs(initial.scale.canonicalTime)+hi));
            if hi-lo<=timeTolerance
                if abs(bestResidual)<=xtol,located=true;end
                break;
            end
            width=hi-lo;probeH=lo-flo*width/(fhi-flo);
            % A nearly linear peak trajectory legitimately puts its
            % next secant very near a bracket endpoint. Preserve it
            % unless rounding leaves the bracket, or two consecutive
            % probes failed to reduce the best residual by 10 percent.
            if ~isfinite(probeH)||probeH<=lo||probeH>=hi||slowProgress>=2
                probeH=lo+width/2;slowProgress=0;
            end
            if probeH<options.minDt||initial.scale.canonicalTime+probeH<=initial.scale.canonicalTime
                error('ipm_cr:PeakEventMinimumStep', ...
                    'A located peak-cell event requires dt below minDt.');
            end
            [candidateState,candidateRecords]=probe(probeH);
            snapshot=ipm_cr.time.peak(candidateState.rhs);
            value=crossing.direction*(snapshot.x-crossing.knot);
            if ~isfinite(value)
                error('ipm_cr:PeakEventUnresolved','The endpoint peak is nonfinite during event localization.');
            end
            previousBest=abs(bestResidual);
            if abs(value)<previousBest
                bestState=candidateState;bestRecords=candidateRecords;
                bestH=probeH;bestResidual=value;
            end
            if abs(bestResidual)>.9*previousBest,slowProgress=slowProgress+1;
            else,slowProgress=0;end
            if value<0,lo=probeH;flo=value;
            else,hi=probeH;fhi=value;end
        end
        located=located||abs(bestResidual)<=xtol;
        if ~located
            error('ipm_cr:PeakEventUnresolved', ...
                ['The original maximum did not continuously reach knot %.17g; ' ...
                'best x residual %.3e. A nonunique peak or unresolved event cannot be projected.'], ...
                crossing.knot,bestResidual);
        end
        s=bestState;r=bestRecords;h=bestH;
        event.located=true;
        event.rootResidual=crossing.direction*bestResidual;
        event.rootTimeBracket=[lo,hi];
        event.status='located_peak_cell_boundary';break;
    end
    event.limited=h<event.initialProposedDt;
    event.acceptedDt=h;event.probeCount=probes;
    event.stageGuardCuts=guardCuts;

    function [probedState,probedRecords]=probe(step)
        probes=probes+1;
        [probedState,probedRecords]=trialProbe(step,probes);
    end
    function crossing=find_crossing(r,h,x0)
        crossing=struct('found',false,'knot',NaN,'direction',0);
        firstTime=Inf;
        % Exclude only a local ulp band around the accepted starting knot.
        % A step leaving a located event then seeks the NEXT knot, avoiding
        % repeated zero-length steps without changing the measured peak.
        if numel(wallKnots)<3,return;end
        knots=wallKnots(2:end-1);
        tolerance=64*eps(max(1,abs(knots)));
        for j=2:numel(r.peaks)
            p=r.peaks(j);
            if ~p.available||~isfinite(p.x)||r.abscissae(j)<=0,continue;end
            if p.x>x0
                index=find(knots>x0+tolerance&knots<=p.x+tolerance,1,'first');direction=1;
            elseif p.x<x0
                index=find(knots<x0-tolerance&knots>=p.x-tolerance,1,'last');direction=-1;
            else,continue;end
            if isempty(index),continue;end
            estimate=h*r.abscissae(j)*abs((knots(index)-x0)/(p.x-x0));
            if estimate<firstTime
                crossing=struct('found',true,'knot',knots(index),'direction',direction);
                firstTime=estimate;
            end
        end
    end
end
