function fit = canonicalMaximumGrowthFit(canonicalTau,maximumValues)
%IPM.DIAGNOSTICS.CANONICALMAXIMUMGROWTHFIT Classify M versus canonical tau.
%   FIT = IPM.DIAGNOSTICS.CANONICALMAXIMUMGROWTHFIT(TAU,M) compares two
%   positive-growth laws on the same log(M) error domain:
%
%     linear maximum:      M = A*(1+k*(tau-tau0)),
%     exponential maximum: M = A*exp(lambda*(tau-tau0)).
%
%   A quadratic model for log(M) is included only as a curvature guard.  A
%   significant, reproducible preference for that model downgrades the
%   window to CURVED_NONASYMPTOTIC instead of forcing a linear/exponential
%   choice.  It is not reported as a third asymptotic growth law.
%
%   Fits use at most 64 uniformly spaced canonical-tau scoring points,
%   canonical-tau trapezoidal weights, a 70/30 blocked validation split,
%   three robust local-quadratic estimates of d(log(M))/d(tau), and trailing
%   windows selected by 0.25/0.5/1/1.5 accumulated e-folds.  Decisions must
%   also survive uniform canonical-tau resampling of the full input.

canonicalTau = canonicalTau(:);
maximumValues = maximumValues(:);
if numel(canonicalTau) ~= numel(maximumValues)
    error('ipm:CanonicalMaximumGrowthSize', ...
        ['Canonical-tau and maximum-value vectors must have the same ' ...
        'length.']);
end
fit = empty_fit();
if any(~isfinite(canonicalTau)) || any(~isfinite(maximumValues)) || ...
        any(maximumValues <= 0)
    error('ipm:CanonicalMaximumGrowthValues', ...
        ['Canonical tau and maximum values must be finite, with positive ' ...
        'maxima.']);
end
if any(diff(canonicalTau) <= 0)
    error('ipm:CanonicalMaximumGrowthTau', ...
        'Canonical tau must be strictly increasing.');
end
if numel(canonicalTau) < 16
    return;
end

logMaximum = log(maximumValues);
tauSpan = canonicalTau(end)-canonicalTau(1);
bandwidthFactors = [0.75,1,1.5];
bandwidths = 0.03*tauSpan*bandwidthFactors;
[rates,rateDerivatives] = local_log_rates( ...
    canonicalTau,logMaximum,bandwidths);

requestedLogGrowth = [0.25,0.5,1,1.5];
windows = build_growth_windows(canonicalTau,logMaximum,rates, ...
    bandwidths,requestedLogGrowth);

samplingGridPoints = min(numel(canonicalTau),512);
samplingTau = linspace( ...
    canonicalTau(1),canonicalTau(end),samplingGridPoints)';
samplingLogMaximum = interp1( ...
    canonicalTau,logMaximum,samplingTau,'pchip');
[samplingRates,~] = local_log_rates( ...
    samplingTau,samplingLogMaximum,bandwidths);
samplingWindows = build_growth_windows(samplingTau,samplingLogMaximum, ...
    samplingRates,bandwidths,requestedLogGrowth);
for index = 1:numel(windows)
    match = find([samplingWindows.requestedLogGrowth] == ...
        windows(index).requestedLogGrowth,1);
    windows(index).preSamplingDecision = windows(index).modelDecision;
    if isempty(match)
        windows(index).uniformInputDecision = 'unavailable';
        windows(index).inputSamplingStable = false;
    else
        windows(index).uniformInputDecision = ...
            samplingWindows(match).modelDecision;
        windows(index).inputSamplingStable = strcmp( ...
            windows(index).modelDecision, ...
            samplingWindows(match).modelDecision);
    end
    if ~windows(index).inputSamplingStable
        windows(index).modelDecision = 'indistinguishable';
        windows(index).curvatureGuardTriggered = false;
    end
end

totalLogGrowth = logMaximum(end)-logMaximum(1);
fit.valid = true;
fit.points = numel(canonicalTau);
fit.coordinateName = 'canonical_tau';
fit.canonicalTau = canonicalTau;
fit.maximum = maximumValues;
fit.logMaximum = logMaximum;
fit.canonicalLogarithmicRate = rates;
fit.canonicalLogarithmicRateDerivative = rateDerivatives;
fit.bandwidthsCanonicalTau = bandwidths;
fit.bandwidthFactors = bandwidthFactors;
fit.rateEstimator = 'robust_local_quadratic_canonical_tau';
fit.modelErrorDomain = 'log_maximum';
fit.coordinateWeighting = 'trapezoidal_canonical_tau';
fit.maximumComparisonPoints = 64;
fit.validationTrainingFraction = 0.7;
fit.requestedLogGrowth = requestedLogGrowth;
fit.samplingSensitivityGridPoints = samplingGridPoints;
fit.minimumCredibleEFolds = 3;
fit.totalLogGrowth = totalLogGrowth;
fit.eFoldCount = totalLogGrowth;
fit.enoughDynamicRange = totalLogGrowth >= fit.minimumCredibleEFolds;
fit.finalCanonicalRate = rates(end,:);
fit.finalCanonicalRateDerivative = rateDerivatives(end,:);
fit.windows = windows;
fit.curvatureGuardTriggered = any([windows.curvatureGuardTriggered]);
[fit.modelConsensus,fit.consensusModel] = window_consensus(windows);
fit.credibleClassification = false;
fit.credibilityLimitations = { ...
    'Canonical-tau derivatives are correlated local-polynomial diagnostics.', ...
    'A curvature-guard decision is nonasymptotic, not a third growth law.', ...
    ['Linear/exponential preference must persist across e-fold windows, ' ...
    'bandwidths, input sampling, grids, time steps, and boxes.'], ...
    'At least three e-folds are required for a credible classification.'};
end

function windows = build_growth_windows( ...
        tau,logMaximum,rates,bandwidths,requestedLogGrowth)
totalLogGrowth = logMaximum(end)-logMaximum(1);
windows = repmat(empty_window(),0,1);
for index = 1:numel(requestedLogGrowth)
    target = requestedLogGrowth(index);
    if totalLogGrowth < target
        continue;
    end
    [windowTau,windowLogMaximum,windowRates] = growth_window( ...
        tau,logMaximum,rates,target);
    if isempty(windowTau)
        continue;
    end
    windows(end+1,1) = fit_window(windowTau,windowLogMaximum, ...
        windowRates,bandwidths,target); %#ok<AGROW>
end
end

function window = fit_window(tau,y,rates,bandwidths,requestedGrowth)
rateTau = tau;
originalPointCount = numel(tau);
[tau,y] = fixed_tau_comparison_series(tau,y,64);

linearMaximum = fit_linear_maximum(tau,y);
exponentialMaximum = fit_exponential_maximum(tau,y);
quadraticGuard = fit_quadratic_log_guard(tau,y);
models = {linearMaximum,exponentialMaximum,quadraticGuard};

rawAicc = cellfun(@(model)model.aicc,models);
[~,rawPreferredIndex] = min(rawAicc);
for index = 1:numel(models)
    models{index}.rawDeltaAicc = ...
        rawAicc(index)-rawAicc(rawPreferredIndex);
    predictedRate = model_canonical_rate(models{index},rateTau);
    models{index}.canonicalRateLogRmse = rate_log_rmse( ...
        rates,predictedRate,rateTau);
end

selectionAicc = rawAicc;
for index = 1:numel(models)
    if ~models{index}.valid || ~models{index}.identifiable
        selectionAicc(index) = inf;
    end
end
[preferredIndex,preferredAicc] = finite_minimum(selectionAicc);
for index = 1:numel(models)
    if isfinite(preferredAicc) && isfinite(selectionAicc(index))
        models{index}.deltaAicc = selectionAicc(index)-preferredAicc;
    else
        models{index}.deltaAicc = inf;
    end
end

trainingEndTau = tau(1)+0.7*(tau(end)-tau(1));
trainCount = find(tau <= trainingEndTau,1,'last');
trainCount = max(6,trainCount);
trainCount = min(trainCount,numel(tau)-3);
train = 1:trainCount;
validation = trainCount+1:numel(tau);
trainingModels = { ...
    fit_linear_maximum(tau(train),y(train)), ...
    fit_exponential_maximum(tau(train),y(train)), ...
    fit_quadratic_log_guard(tau(train),y(train))};
for index = 1:numel(models)
    prediction = model_log_value(trainingModels{index},tau(validation));
    models{index}.validationLogRmse = weighted_rmse( ...
        y(validation)-prediction,tau_quadrature_weights(tau(validation)));
end
validationRmse = cellfun(@(model)model.validationLogRmse,models);
[~,rawValidationIndex] = min(validationRmse);
selectionValidation = validationRmse;
for index = 1:numel(models)
    if ~models{index}.valid || ~models{index}.identifiable || ...
            ~trainingModels{index}.valid || ...
            ~trainingModels{index}.identifiable
        selectionValidation(index) = inf;
    end
end
[validationIndex,validationMinimum] = finite_minimum(selectionValidation);

ratePreferred = cell(1,numel(bandwidths));
for bandwidthIndex = 1:numel(bandwidths)
    errors = cellfun( ...
        @(model)model.canonicalRateLogRmse(bandwidthIndex),models);
    for modelIndex = 1:numel(models)
        if ~models{modelIndex}.valid || ~models{modelIndex}.identifiable
            errors(modelIndex) = inf;
        end
    end
    [rateIndex,rateMinimum] = finite_minimum(errors);
    if isfinite(rateMinimum)
        ratePreferred{bandwidthIndex} = models{rateIndex}.name;
    else
        ratePreferred{bandwidthIndex} = 'none';
    end
end

window = empty_window();
window.requestedLogGrowth = requestedGrowth;
window.actualLogGrowth = y(end)-y(1);
window.startCanonicalTau = tau(1);
window.endCanonicalTau = tau(end);
window.points = originalPointCount;
window.comparisonPoints = numel(tau);
window.linearMaximum = models{1};
window.exponentialMaximum = models{2};
window.quadraticLogGuard = models{3};
window.rawPreferredByAicc = models{rawPreferredIndex}.name;
if isfinite(preferredAicc)
    window.preferredByAicc = models{preferredIndex}.name;
end
window.rawPreferredByValidation = models{rawValidationIndex}.name;
if isfinite(validationMinimum)
    window.preferredByValidation = models{validationIndex}.name;
end
window.preferredByRateBandwidth = ratePreferred;
if isfinite(preferredAicc)
    window.rateBandwidthAgreementCount = nnz(strcmp( ...
        ratePreferred,models{preferredIndex}.name));
end
window.rateBandwidthConsensus = ...
    window.rateBandwidthAgreementCount >= ceil(numel(bandwidths)/2);
window.bandwidthsCanonicalTau = bandwidths;

finiteAicc = sort(selectionAicc(isfinite(selectionAicc)));
if numel(finiteAicc) >= 2
    window.aiccSeparation = finiteAicc(2)-finiteAicc(1);
elseif isscalar(finiteAicc)
    window.aiccSeparation = inf;
end

if isfinite(preferredAicc) && preferredIndex <= 2
    alternativeIndex = 3-preferredIndex;
    if isfinite(selectionValidation(alternativeIndex))
        window.validationImprovementOverAlternative = 1- ...
            selectionValidation(preferredIndex) / ...
            max(selectionValidation(alternativeIndex),eps);
    end
end

sameValidationChoice = preferredIndex == validationIndex && ...
    isfinite(preferredAicc) && isfinite(validationMinimum);
guardSelected = preferredIndex == 3 && sameValidationChoice && ...
    window.aiccSeparation >= 10 && ...
    models{3}.curvatureAcrossWindow >= 0.02 && ...
    window.rateBandwidthConsensus;
if guardSelected
    window.curvatureGuardTriggered = true;
    window.modelDecision = 'curved_nonasymptotic';
elseif preferredIndex <= 2 && sameValidationChoice && ...
        window.aiccSeparation >= 10 && ...
        window.rateBandwidthConsensus && ...
        window.validationImprovementOverAlternative >= 0.2
    window.modelDecision = models{preferredIndex}.name;
end
end

function model = fit_linear_maximum(tau,y)
referenceTau = tau(1);
x = tau-referenceTau;
span = max(x(end),eps);
lowerBound = log(1e-8);
upperBound = log(1e6);
weights = tau_quadrature_weights(tau);
objective = @(logTheta) linear_maximum_objective( ...
    logTheta,x,y,weights,span);
options = optimset('Display','off','TolX',1e-12);
[logTheta,~,exitFlag] = fminbnd( ...
    objective,lowerBound,upperBound,options);
k = exp(logTheta)/span;
shape = log1p(k*x);
logA = weighted_constant(y-shape,weights);
prediction = logA+shape;
model = base_model('linear_maximum',2,tau,y,prediction,weights);
model.referenceCanonicalTau = referenceTau;
model.parameters = struct( ...
    'logA',logA,'A',exp(logA),'k',k,'B',exp(logA)*k);
model.parameterBounds = [exp(lowerBound)/span,exp(upperBound)/span];
model.profileBoundaryHit = profile_boundary_hit( ...
    logTheta,lowerBound,upperBound);
model.optimizerExitFlag = exitFlag;
model.valid = k > 0 && isfinite(logA);
model.identifiable = model.valid && ~model.profileBoundaryHit;
if ~model.valid
    model.aicc = inf;
end
end

function value = linear_maximum_objective(logTheta,x,y,weights,span)
k = exp(logTheta)/span;
shape = log1p(k*x);
logA = weighted_constant(y-shape,weights);
residual = y-(logA+shape);
value = sum(weights.*residual.^2);
end

function model = fit_exponential_maximum(tau,y)
referenceTau = tau(1);
x = tau-referenceTau;
design = [ones(numel(tau),1),x];
weights = tau_quadrature_weights(tau);
coefficients = weighted_fit(design,y,weights);
model = base_model( ...
    'exponential_maximum',2,tau,y,design*coefficients,weights);
model.referenceCanonicalTau = referenceTau;
model.parameters = struct('logA',coefficients(1), ...
    'A',exp(coefficients(1)),'lambda',coefficients(2));
model.valid = coefficients(2) > 0;
model.identifiable = model.valid;
if ~model.valid
    model.aicc = inf;
end
end

function model = fit_quadratic_log_guard(tau,y)
referenceTau = tau(1);
x = tau-referenceTau;
design = [ones(numel(tau),1),x,x.^2];
weights = tau_quadrature_weights(tau);
coefficients = weighted_fit(design,y,weights);
prediction = design*coefficients;
model = base_model( ...
    'quadratic_log_guard',3,tau,y,prediction,weights);
model.referenceCanonicalTau = referenceTau;
model.parameters = struct('logA',coefficients(1), ...
    'A',exp(coefficients(1)),'lambda',coefficients(2), ...
    'quadraticCoefficient',coefficients(3));
span = x(end);
model.curvatureAcrossWindow = abs(coefficients(3))*span^2;
endpointRates = coefficients(2) + ...
    2*coefficients(3)*[x(1),x(end)];
model.valid = all(endpointRates > 0);
model.identifiable = model.valid && ...
    model.curvatureAcrossWindow >= 0.02;
if ~model.valid
    model.aicc = inf;
end
end

function model = base_model(name,parameterCount,tau,y,prediction,weights)
residual = y-prediction;
sse = sum(weights.*residual.^2);
points = numel(tau);
weightedScale = sum(weights.*max(abs(y),1).^2);
sseForInformationCriterion = max(sse,eps^2*weightedScale);
model = empty_model(name);
model.valid = true;
model.parameterCount = parameterCount;
model.sseLogValue = sse;
model.logValueRmse = weighted_rmse(residual,weights);
model.aicc = points*log(sseForInformationCriterion/points) + ...
    2*parameterCount + ...
    2*parameterCount*(parameterCount+1) / ...
    max(points-parameterCount-1,1);
end

function values = model_log_value(model,tau)
x = tau-model.referenceCanonicalTau;
switch model.name
    case 'linear_maximum'
        values = model.parameters.logA+log1p(model.parameters.k*x);
    case 'exponential_maximum'
        values = model.parameters.logA+model.parameters.lambda*x;
    case 'quadratic_log_guard'
        values = model.parameters.logA + model.parameters.lambda*x + ...
            model.parameters.quadraticCoefficient*x.^2;
    otherwise
        values = NaN(size(tau));
end
end

function values = model_canonical_rate(model,tau)
x = tau-model.referenceCanonicalTau;
switch model.name
    case 'linear_maximum'
        values = model.parameters.k./(1+model.parameters.k*x);
    case 'exponential_maximum'
        values = repmat(model.parameters.lambda,size(tau));
    case 'quadratic_log_guard'
        values = model.parameters.lambda + ...
            2*model.parameters.quadraticCoefficient*x;
    otherwise
        values = NaN(size(tau));
end
end

function errors = rate_log_rmse(rates,prediction,tau)
errors = NaN(1,size(rates,2));
for index = 1:size(rates,2)
    valid = isfinite(rates(:,index)) & rates(:,index) > 0 & ...
        isfinite(prediction) & prediction > 0;
    if any(valid)
        residual = log(rates(valid,index))-log(prediction(valid));
        errors(index) = weighted_rmse( ...
            residual,tau_quadrature_weights(tau(valid)));
    end
end
end

function [rates,rateDerivatives] = local_log_rates(tau,y,bandwidths)
rates = NaN(numel(tau),numel(bandwidths));
rateDerivatives = NaN(numel(tau),numel(bandwidths));
for index = 1:numel(bandwidths)
    [rates(:,index),rateDerivatives(:,index)] = local_log_rate( ...
        tau,y,bandwidths(index));
end
end

function [rate,rateDerivative] = local_log_rate(tau,y,bandwidth)
points = numel(tau);
minimumPoints = min(7,points);
rate = NaN(points,1);
rateDerivative = NaN(points,1);
for center = 1:points
    distance = abs(tau-tau(center));
    selection = find(distance <= bandwidth);
    if numel(selection) < minimumPoints
        [~,order] = sort(distance,'ascend');
        selection = sort(order(1:minimumPoints));
    end
    localScale = max(distance(selection));
    if localScale == 0
        continue;
    end
    coordinate = (tau(selection)-tau(center))/localScale;
    kernel = max(1-abs(coordinate).^3,0).^3;
    kernel = max(kernel,1e-6);
    design = [ones(numel(selection),1),coordinate,coordinate.^2];
    quadratureWeights = tau_quadrature_weights(tau(selection));
    weights = kernel.*quadratureWeights;
    coefficients = weighted_fit(design,y(selection),weights);
    for iteration = 1:3
        residual = y(selection)-design*coefficients;
        residualScale = 1.4826*median(abs( ...
            residual-median(residual)))+eps;
        huber = min(1,1.345*residualScale./max(abs(residual),eps));
        weights = kernel.*quadratureWeights.*huber;
        coefficients = weighted_fit(design,y(selection),weights);
    end
    rate(center) = coefficients(2)/localScale;
    rateDerivative(center) = 2*coefficients(3)/localScale^2;
end
end

function [windowTau,windowValues,windowRates] = growth_window( ...
        tau,values,rates,target)
threshold = values(end)-target;
first = find(values <= threshold,1,'last');
if isempty(first) || first >= numel(values)
    windowTau = [];
    windowValues = [];
    windowRates = [];
    return;
end
latestStart = numel(tau)-9;
if first >= latestStart
    first = latestStart;
    windowTau = tau(first:end);
    windowValues = values(first:end);
    windowRates = rates(first:end,:);
    return;
end
denominator = values(first+1)-values(first);
if denominator <= 0
    windowTau = [];
    windowValues = [];
    windowRates = [];
    return;
end
fraction = (threshold-values(first))/denominator;
crossingTau = tau(first)+fraction*(tau(first+1)-tau(first));
crossingRates = rates(first,:) + fraction* ...
    (rates(first+1,:)-rates(first,:));
windowTau = [crossingTau;tau(first+1:end)];
windowValues = [threshold;values(first+1:end)];
windowRates = [crossingRates;rates(first+1:end,:)];
end

function [comparisonTau,comparisonValues] = ...
        fixed_tau_comparison_series(tau,values,targetPoints)
comparisonPoints = min(numel(tau),targetPoints);
comparisonTau = linspace(tau(1),tau(end),comparisonPoints)';
if comparisonPoints == numel(tau) && ...
        all(abs(comparisonTau-tau) <= ...
        64*eps(max(1,max(abs(tau)))))
    comparisonValues = values;
else
    comparisonValues = interp1(tau,values,comparisonTau,'pchip');
end
end

function weights = tau_quadrature_weights(tau)
points = numel(tau);
if points <= 1
    weights = ones(size(tau));
    return;
end
spacing = diff(tau);
weights = [spacing(1)/2; ...
    (spacing(1:end-1)+spacing(2:end))/2;spacing(end)/2];
weights = weights*points/sum(weights);
end

function value = weighted_constant(values,weights)
value = sum(weights.*values)/sum(weights);
end

function coefficients = weighted_fit(design,values,weights)
weightedDesign = design.*sqrt(weights);
weightedValues = values.*sqrt(weights);
coefficients = weightedDesign\weightedValues;
end

function value = weighted_rmse(residual,weights)
value = sqrt(sum(weights.*residual.^2)/sum(weights));
end

function hit = profile_boundary_hit(value,lowerBound,upperBound)
tolerance = 1e-3*(upperBound-lowerBound);
hit = value-lowerBound <= tolerance || ...
    upperBound-value <= tolerance;
end

function [index,value] = finite_minimum(values)
finiteValues = values;
finiteValues(~isfinite(finiteValues)) = inf;
[value,index] = min(finiteValues);
if isempty(value) || isinf(value)
    index = 1;
    value = inf;
end
end

function [hasConsensus,modelName] = window_consensus(windows)
hasConsensus = false;
modelName = 'none';
if isempty(windows)
    return;
end
decisions = {windows.modelDecision};
if all(strcmp(decisions,'curved_nonasymptotic'))
    modelName = 'curved_nonasymptotic';
    return;
end
if numel(windows) >= 2 && ...
        any(strcmp(decisions{1}, ...
        {'linear_maximum','exponential_maximum'})) && ...
        all(strcmp(decisions,decisions{1}))
    hasConsensus = true;
    modelName = decisions{1};
else
    modelName = 'window_dependent';
end
end

function fit = empty_fit()
fit = struct('valid',false,'points',0, ...
    'coordinateName','canonical_tau','canonicalTau',[], ...
    'maximum',[],'logMaximum',[], ...
    'canonicalLogarithmicRate',[], ...
    'canonicalLogarithmicRateDerivative',[], ...
    'bandwidthsCanonicalTau',[],'bandwidthFactors',[], ...
    'rateEstimator','','modelErrorDomain','', ...
    'coordinateWeighting','','maximumComparisonPoints',64, ...
    'validationTrainingFraction',0.7, ...
    'requestedLogGrowth',[0.25,0.5,1,1.5], ...
    'samplingSensitivityGridPoints',0,'minimumCredibleEFolds',3, ...
    'totalLogGrowth',NaN,'eFoldCount',NaN, ...
    'enoughDynamicRange',false,'finalCanonicalRate',[], ...
    'finalCanonicalRateDerivative',[], ...
    'windows',repmat(empty_window(),0,1), ...
    'curvatureGuardTriggered',false, ...
    'modelConsensus',false,'consensusModel','none', ...
    'credibleClassification',false,'credibilityLimitations',{{}});
end

function window = empty_window()
window = struct('requestedLogGrowth',NaN,'actualLogGrowth',NaN, ...
    'startCanonicalTau',NaN,'endCanonicalTau',NaN, ...
    'points',0,'comparisonPoints',0, ...
    'linearMaximum',empty_model('linear_maximum'), ...
    'exponentialMaximum',empty_model('exponential_maximum'), ...
    'quadraticLogGuard',empty_model('quadratic_log_guard'), ...
    'rawPreferredByAicc','none','preferredByAicc','none', ...
    'rawPreferredByValidation','none','preferredByValidation','none', ...
    'preferredByRateBandwidth',{{}}, ...
    'rateBandwidthAgreementCount',0,'rateBandwidthConsensus',false, ...
    'bandwidthsCanonicalTau',[],'aiccSeparation',NaN, ...
    'validationImprovementOverAlternative',NaN, ...
    'curvatureGuardTriggered',false, ...
    'preSamplingDecision','indistinguishable', ...
    'uniformInputDecision','unavailable', ...
    'inputSamplingStable',false,'modelDecision','indistinguishable');
end

function model = empty_model(name)
model = struct('name',name,'valid',false,'identifiable',false, ...
    'parameterCount',0,'parameters',struct(), ...
    'referenceCanonicalTau',NaN,'sseLogValue',NaN, ...
    'logValueRmse',NaN,'aicc',inf,'rawDeltaAicc',NaN, ...
    'deltaAicc',NaN,'validationLogRmse',NaN, ...
    'canonicalRateLogRmse',[],'curvatureAcrossWindow',NaN, ...
    'parameterBounds',[],'profileBoundaryHit',false, ...
    'optimizerExitFlag',NaN);
end
