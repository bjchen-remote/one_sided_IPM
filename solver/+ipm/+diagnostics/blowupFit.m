function fit = blowupFit(t,gradInf)
%IPM.DIAGNOSTICS.BLOWUPFIT Compare exponential growth with A*(T-t)^(-p).
%   This is a diagnostic, not proof of singularity. The primary report uses
%   the last 40 valid points; 20/40/80/160-point fits expose window drift.

fit = empty_fit();
t = t(:);
gradInf = gradInf(:);
if numel(t) ~= numel(gradInf)
    error('ipm:BlowupFitSize', ...
        'Time and maximum-value vectors must have the same length.');
end
valid = isfinite(t) & isfinite(gradInf) & gradInf > 0;
t = t(valid);
g = gradInf(valid);
if numel(t) < 12 || t(end) <= t(1)
    return;
end

primary = fit_window(t,g,min(40,numel(t)));
fit.candidate = primary.candidate;
fit.T = primary.T;
fit.p = primary.p;
fit.r2Power = primary.r2Power;
fit.r2Exponential = primary.r2Exponential;
fit.lambda = primary.lambda;
fit.powerProfileBoundaryHit = primary.powerProfileBoundaryHit;
fit.points = primary.points;

requestedWindows = [20,40,80,160];
windows = requestedWindows(requestedWindows <= numel(t));
if isempty(windows) || windows(end) ~= numel(t)
    windows = unique([windows,min(numel(t),160)]);
end
fit.windowPoints = windows;
fit.windowT = NaN(size(windows));
fit.windowP = NaN(size(windows));
fit.windowR2Power = NaN(size(windows));
fit.windowR2Exponential = NaN(size(windows));
fit.windowLambda = NaN(size(windows));
fit.windowPowerProfileBoundaryHit = false(size(windows));
fit.windowCandidate = false(size(windows));
for index = 1:numel(windows)
    windowFit = fit_window(t,g,windows(index));
    fit.windowT(index) = windowFit.T;
    fit.windowP(index) = windowFit.p;
    fit.windowR2Power(index) = windowFit.r2Power;
    fit.windowR2Exponential(index) = windowFit.r2Exponential;
    fit.windowLambda(index) = windowFit.lambda;
    fit.windowPowerProfileBoundaryHit(index) = ...
        windowFit.powerProfileBoundaryHit;
    fit.windowCandidate(index) = windowFit.candidate;
end
end

function fit = fit_window(t,g,n)
t = t(end-n+1:end);
g = g(end-n+1:end);
fit = empty_fit();
fit.points = n;
logg = log(g(:));

exponentialMatrix = [ones(n,1),t(:)];
exponentialCoefficients = exponentialMatrix\logg;
fit.lambda = exponentialCoefficients(2);
exponentialError = sum((logg-exponentialMatrix*exponentialCoefficients).^2);
totalVariation = sum((logg-mean(logg)).^2);
fit.r2Exponential = 1-exponentialError/max(totalVariation,eps);

span = max(t(end)-t(1),eps);
minimumOffset = max(span/1e4,eps);
maximumOffset = 1e6*span;
lowerBound = log(minimumOffset);
upperBound = log(maximumOffset);
objective = @(logOffset) power_error( ...
    logOffset,t(:),logg,t(end));
logOffset = fminbnd(objective,lowerBound,upperBound, ...
    optimset('Display','off','TolX',1e-10));
offset = exp(logOffset);
fit.T = t(end)+offset;
singularCoordinate = -log(fit.T-t(:));
powerMatrix = [ones(n,1),singularCoordinate];
coefficients = powerMatrix\logg;
fit.p = coefficients(2);
bestError = sum((logg-powerMatrix*coefficients).^2);
fit.powerProfileBoundaryHit = ...
    logOffset-lowerBound <= 1e-3*(upperBound-lowerBound) || ...
    upperBound-logOffset <= 1e-3*(upperBound-lowerBound);
fit.r2Power = 1-bestError/max(totalVariation,eps);
fit.candidate = isfinite(fit.r2Power) && fit.r2Power > 0.98 && ...
    fit.r2Power > fit.r2Exponential+0.01 && ...
    fit.p > 0 && ~fit.powerProfileBoundaryHit;
end

function errorValue = power_error(logOffset,t,logg,minimumTerminalTime)
terminalTime = minimumTerminalTime+exp(logOffset);
coordinate = -log(terminalTime-t);
design = [ones(numel(t),1),coordinate];
coefficients = design\logg;
if coefficients(2) <= 0
    errorValue = inf;
else
    residual = logg-design*coefficients;
    errorValue = sum(residual.^2);
end
end

function fit = empty_fit()
fit = struct('candidate',false,'T',NaN,'p',NaN,'r2Power',NaN, ...
    'r2Exponential',NaN,'lambda',NaN, ...
    'powerProfileBoundaryHit',false,'points',0,'windowPoints',[], ...
    'windowT',[],'windowP',[],'windowR2Power',[], ...
    'windowR2Exponential',[],'windowLambda',[], ...
    'windowPowerProfileBoundaryHit',[],'windowCandidate',[]);
end
