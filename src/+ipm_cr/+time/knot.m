classdef knot
    %KNOT Peak crossings of the original wall spline.
    methods (Static)
        function p=peak(f)
            % Observe existing measurements only; missing diagnostics never reject a step.
            p=struct('x',NaN,'cellIndex',NaN,'cellFraction',NaN, ...
                'tieNonunique',false,'curvature',NaN,'available',false);
            if ~isstruct(f)||~isscalar(f)||~isfield(f,'measurement')|| ...
                    ~isstruct(f.measurement)||~isfield(f.measurement,'peak')|| ...
                    ~isstruct(f.measurement.peak)||~isscalar(f.measurement.peak)
                return;
            end
            source=f.measurement.peak;names={'x','cellIndex','cellFraction','tieNonunique','curvature'};
            for k=1:numel(names)
                if isfield(source,names{k}),p.(names{k})=source.(names{k});end
            end
            p.available=true;
        end

        function [out,records,h,p]=fit(initial,out,records,h,knots,evaluate,p)
            % Only a unique, locally concave peak crossing adjacent spline cells defines
            % this partition. A disconnected maximizing-branch change evolves directly.
            p.kind='RK4_original_spline_knot_fixed4_secant';
            a=ipm_cr.time.knot.peak(initial.rhs);b=ipm_cr.time.knot.peak(out.rhs);
            direction=sign(b.x-a.x);p.direction=direction;
            if ~a.available||~b.available||direction==0,return;end
            [ca,p.startClassification]=peak_cell(a.x,knots,direction,true);
            [cb,p.endClassification]=peak_cell(b.x,knots,direction,false);
            if cb-ca~=direction||any([records.peaks.tieNonunique])|| ...
                    a.curvature>=0||b.curvature>=0,return;end
            node=ca+(direction>0);[v0,c0]=knot_slope(initial.rhs,node,knots);
            [v1,c1]=knot_slope(out.rhs,node,knots);
            if ~(c0<0&&c1<0&&v0*v1<=0),return;end
            p.fitted=true;p.knotIndex=node;p.knotX=knots(node);
            previousH=0;previousValue=v0;value=v1;
            for iteration=1:4
                % Equal function values do not define a secant. Keep the last computed
                % trial; no new time interval can be obtained from this pair.
                if value==0||value==previousValue,break;end
                nextH=h-value*(h-previousH)/(value-previousValue);
                if nextH==h,break;end
                previousH=h;previousValue=value;h=nextH;
                [out,records]=evaluate(h);value=knot_slope(out.rhs,node,knots);
                p.secantIterations=iteration;
            end
            p.segmentDt=h;p.knotSlope=value;p.rawPeakKnotDefect=out.rhs.measurement.peak.x-knots(node);
        end
    end
end

function [cell,info]=peak_cell(value,knots,direction,isStart)
% Coordinate identity only: outgoing ownership at a start knot, incoming
% ownership at an end knot. Neither the measured peak nor state is changed.
knots=knots(:)';[distance,node]=min(abs(knots-value));
tolerance=256*eps(max([1,abs(value),abs(knots(node))]));onKnot=distance<=tolerance;
classified=value;
if onKnot
    classified=knots(node);cell=node-(direction<0);
    if ~isStart,cell=node-(direction>0);end
else,cell=find(knots<value,1,'last');end
cell=max(1,min(numel(knots)-1,cell));
info=struct('observedX',value,'classifiedX',classified,'coordinateDifference',classified-value, ...
    'onKnot',onKnot,'tolerance',tolerance,'cellIndex',cell);
end
function [slope,curvature]=knot_slope(f,node,knots)
q=f.measurement.q;peak=f.measurement.peak;
if isfield(peak,'evaluationAnchor'),q=q-peak.evaluationAnchor;end
M=peak.splineSecondDerivatives;width=knots(node+1)-knots(node);
slope=(q(node+1)-q(node))/width-width*(2*M(node)+M(node+1))/6;
curvature=M(node);
end
