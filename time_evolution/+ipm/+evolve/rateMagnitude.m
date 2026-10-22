function [rateMagnitude,normalizedRates] = rateMagnitude( ...
    rateVector,translationScale)
%IPM.EVOLVE.RATEMAGNITUDE Measure scaling rates without modifying time.
%   The translation component is normalized by a fixed coordinate scale.
%   This is a diagnostic only: it never clips a rate or changes the
%   canonical-time vector field.

if nargin ~= 2
    error('ipm:RateMagnitudeArguments', ...
        'Use rateVector and translationScale.');
end
if ~isnumeric(rateVector) || isempty(rateVector) || ...
        any(~isfinite(rateVector),'all')
    error('ipm:RateMagnitudeRates', ...
        'rateVector must be a nonempty finite numeric vector.');
end
validateattributes(translationScale,{'numeric'}, ...
    {'scalar','real','finite','positive'});
normalizedRates = abs(rateVector(:));
normalizedRates(end) = normalizedRates(end)/translationScale;
rateMagnitude = max(normalizedRates);
end
