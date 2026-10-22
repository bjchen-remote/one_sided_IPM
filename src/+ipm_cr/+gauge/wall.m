function [measurement,rateEvaluator]=wall(field,r,wallOps,options)
%WALLMEASURE Original I/P/C with a constant-shifted natural-spline peak.
% Subtracting max(q) changes neither the mathematical natural spline nor
% its maximizing location/functional rows. It prevents a P-scaled tie band
% from selecting a nearby node instead of the stationary maximum. In the
% shifted data every node is nonpositive; the local overshoot and its tie
% tolerance then have the same small scale near a knot crossing.
% I, C, all original references, wall sampling, and rate sampling remain
% those of ipm_cr.gauge.measure. The maintained packages are not modified.
if nargin<4,options=struct();end
[measurement,originalRate]=ipm_cr.gauge.measure(field,r,wallOps,options);
anchor=max(measurement.q);
peak=ipm_cr.gauge.peak(measurement.q-anchor,wallOps);
peak.evaluationAnchor=anchor;
peak.centeredMaximumCandidateValue=peak.maximumCandidateValue;
peak.value=peak.value+anchor;
peak.gridValue=peak.gridValue+anchor;
peak.maximumCandidateValue=peak.maximumCandidateValue+anchor;
peak.tieToleranceUsesCenteredValues=true;
assert(isfinite(peak.value)&&peak.value>0, ...
    'ipm_cr:ExplicitWallPeak','The original natural-spline maximum must remain finite and positive.');
measurement.peak=peak;measurement.P=peak.value;
measurement.peakEvaluation='constant_shifted_original_global_natural_cubic_spline';
if ~isempty(measurement.rate),measurement.rate=correct_rate(measurement.rate);end
rateEvaluator=[];
if nargout>1,rateEvaluator=@prepared_rate;end

    function observed=prepared_rate(rhs,rateOptions)
        if nargin<2,rateOptions=struct();end
        original=originalRate(rhs,rateOptions);
        observed=measurement;
        observed.rate=correct_rate(original.rate);
    end
    function rate=correct_rate(rate)
        rate.P_dot=peak.weights*rate.qDot(peak.indices)';
        rate.normalized(2)=rate.P_dot/peak.value;
    end
end
