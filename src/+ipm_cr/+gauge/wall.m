function [measurement,rateEvaluator]=wall(field,r,wallOps,options)
%WALLMEASURE Original I/P/C with a constant-shifted natural-spline peak.
% Subtracting max(q) changes neither the mathematical natural spline nor
% its maximizing location/functional rows. It prevents a P-scaled tie band
% from selecting a nearby node instead of the stationary maximum. In the
% shifted data every node is nonpositive; the local overshoot and its tie
% tolerance then have the same small scale near a knot crossing.
% I, C, all original references, wall sampling, and rate sampling remain
% those of ipm_cr.gauge.measure. Its bound rate closure uses the same single
% centered peak computation. Independent endpoint audits call wall afresh.
if nargin<4,options=struct();end
assert(isstruct(options)&&isscalar(options),'polar:WallOptions', ...
    'Measurement options must be a scalar structure.');
options.centerPeak=true;
[measurement,rateEvaluator]=ipm_cr.gauge.measure(field,r,wallOps,options);
end
