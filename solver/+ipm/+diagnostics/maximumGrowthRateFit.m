function fit = maximumGrowthRateFit(t,maximumValues)
%IPM.DIAGNOSTICS.MAXIMUMGROWTHRATEFIT Fit the logarithmic maximum growth rate.
%   For M(t)>0, the diagnostic estimates gamma=d(log(M))/dt on physical
%   time and compares three models under the same log(M) error model:
%
%     exponential:       gamma = lambda,
%     finite-time power: gamma = p/(T-t),
%     double exponential:gamma = gamma0*exp(kappa*(t-t0)).
%
%   Local rates are descriptive derivative diagnostics.  AICc and blocked
%   validation errors are computed from log(M), rather than from the noisy
%   differentiated series.  Fits use physical-time quadrature weights and
%   at most 64 equally spaced scoring times.  Trailing windows are selected
%   by accumulated log growth, not by solver-step count.  A decision is
%   downgraded when it changes after uniform resampling of the full input.

t = t(:);
maximumValues = maximumValues(:);
if numel(t) ~= numel(maximumValues)
    error('ipm:MaximumGrowthRateSize', ...
        'Time and maximum-value vectors must have the same length.');
end
fit = empty_fit();
if any(~isfinite(t)) || any(~isfinite(maximumValues)) || ...
        any(maximumValues <= 0)
    error('ipm:MaximumGrowthRateValues', ...
        'Times and maximum values must be finite, with positive maxima.');
end
if any(diff(t) <= 0)
    error('ipm:MaximumGrowthRateTime', ...
        'Physical times must be strictly increasing.');
end
if numel(t) < 16
    return;
end

logMaximum = log(maximumValues);
timeSpan = t(end)-t(1);
baseBandwidth = 0.03*timeSpan;
bandwidthFactors = [0.75,1,1.5];
bandwidths = baseBandwidth*bandwidthFactors;
rates = NaN(numel(t),numel(bandwidths));
rateDerivatives = NaN(numel(t),numel(bandwidths));
for index = 1:numel(bandwidths)
    [rates(:,index),rateDerivatives(:,index)] = local_log_rate( ...
        t,logMaximum,bandwidths(index));
end
instantaneousTerminalTime = NaN(size(rates));
instantaneousPower = NaN(size(rates));
accelerating = isfinite(rates) & rates > 0 & ...
    isfinite(rateDerivatives) & rateDerivatives > 0;
timeMatrix = repmat(t,1,numel(bandwidths));
instantaneousTerminalTime(accelerating) = timeMatrix(accelerating) + ...
    rates(accelerating)./rateDerivatives(accelerating);
instantaneousPower(accelerating) = rates(accelerating).^2 ./ ...
    rateDerivatives(accelerating);

requestedLogGrowth = [0.25,0.5,1,1.5];
totalLogGrowth = logMaximum(end)-logMaximum(1);
windows = build_growth_windows( ...
    t,logMaximum,rates,bandwidths,requestedLogGrowth);
samplingGridPoints = min(numel(t),512);
samplingTime = linspace(t(1),t(end),samplingGridPoints)';
samplingLogMaximum = interp1(t,logMaximum,samplingTime,'pchip');
samplingRates = interp1(t,rates,samplingTime,'pchip');
samplingWindows = build_growth_windows(samplingTime, ...
    samplingLogMaximum,samplingRates,bandwidths,requestedLogGrowth);
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
    end
end

fit.valid = true;
fit.points = numel(t);
fit.time = t;
fit.maximum = maximumValues;
fit.logMaximum = logMaximum;
fit.logarithmicRate = rates;
fit.logarithmicRateDerivative = rateDerivatives;
fit.instantaneousTerminalTime = instantaneousTerminalTime;
fit.instantaneousPower = instantaneousPower;
fit.bandwidths = bandwidths;
fit.bandwidthFactors = bandwidthFactors;
fit.rateEstimator = 'robust_local_quadratic';
fit.modelErrorDomain = 'log_maximum';
fit.timeWeighting = 'trapezoidal_physical_time';
fit.maximumComparisonPoints = 64;
fit.samplingSensitivityGridPoints = samplingGridPoints;
fit.minimumCredibleEFolds = 3;
fit.totalLogGrowth = totalLogGrowth;
fit.eFoldCount = totalLogGrowth;
fit.enoughDynamicRange = totalLogGrowth >= 3;
fit.finalRate = rates(end,:);
fit.finalRateDerivative = rateDerivatives(end,:);
fit.finalInstantaneousTerminalTime = instantaneousTerminalTime(end,:);
fit.finalInstantaneousPower = instantaneousPower(end,:);
fit.windows = windows;
if isempty(windows)
    fit.modelConsensus = false;
    fit.consensusModel = 'none';
else
    preferred = {windows.modelDecision};
    fit.modelConsensus = numel(windows) >= 2 && ...
        ~strcmp(preferred{1},'indistinguishable') && ...
        all(strcmp(preferred,preferred{1}));
    if fit.modelConsensus
        fit.consensusModel = preferred{1};
    else
        fit.consensusModel = 'window_dependent';
    end
end
fit.credibleClassification = false;
fit.credibilityLimitations = { ...
    'Derivative estimates are correlated local-polynomial diagnostics.', ...
    'Equally spaced scoring values remain correlated samples of one run.', ...
    ['Model preference must be stable across growth windows, bandwidths, ' ...
    'and raw versus uniform input sampling.'], ...
    'At least three e-folds and grid/time-step/box replication are required.'};
end

function window = fit_window(t,y,rates,bandwidths,requestedGrowth)
rateTime = t;
originalPointCount = numel(t);
[t,y] = fixed_time_comparison_series(t,y,64);
timeSpan = t(end)-t(1);
exponential = fit_exponential(t,y);
powerLaw = fit_power_law(t,y,t(end),timeSpan);
doubleExponential = fit_double_exponential(t,y,timeSpan);
models = {exponential,powerLaw,doubleExponential};

aicc = cellfun(@(model)model.aicc,models);
[~,rawPreferredIndex] = min(aicc);
for index = 1:numel(models)
    models{index}.rawDeltaAicc = aicc(index)-aicc(rawPreferredIndex);
    predictedRate = model_rate(models{index},rateTime);
    models{index}.rateLogRmse = rate_log_rmse( ...
        rates,predictedRate,rateTime);
end

selectionAicc = aicc;
for index = 1:numel(models)
    if ~models{index}.valid || ~models{index}.identifiable
        selectionAicc(index) = inf;
    end
end
[~,preferredIndex] = min(selectionAicc);
for index = 1:numel(models)
    if isfinite(selectionAicc(index))
        models{index}.deltaAicc = ...
            selectionAicc(index)-selectionAicc(preferredIndex);
    else
        models{index}.deltaAicc = inf;
    end
end

trainingEndTime = t(1)+0.7*(t(end)-t(1));
trainCount = find(t <= trainingEndTime,1,'last');
trainCount = max(6,trainCount);
trainCount = min(trainCount,numel(t)-3);
train = 1:trainCount;
validation = trainCount+1:numel(t);
trainingModels = { ...
    fit_exponential(t(train),y(train)), ...
    fit_power_law(t(train),y(train),t(end),timeSpan), ...
    fit_double_exponential(t(train),y(train),timeSpan)};
for index = 1:numel(models)
    prediction = model_log_value(trainingModels{index},t(validation));
    models{index}.validationLogRmse = weighted_rmse( ...
        y(validation)-prediction,time_quadrature_weights(t(validation)));
end
validationRmse = cellfun( ...
    @(model)model.validationLogRmse,models);
[~,rawValidationIndex] = min(validationRmse);
selectionValidation = validationRmse;
for index = 1:numel(models)
    if ~models{index}.valid || ~models{index}.identifiable || ...
            ~trainingModels{index}.valid || ...
            ~trainingModels{index}.identifiable
        selectionValidation(index) = inf;
    end
end
[~,validationIndex] = min(selectionValidation);

ratePreferred = cell(1,numel(bandwidths));
for bandwidthIndex = 1:numel(bandwidths)
    errors = cellfun( ...
        @(model)model.rateLogRmse(bandwidthIndex),models);
    for modelIndex = 1:numel(models)
        if ~models{modelIndex}.valid || ~models{modelIndex}.identifiable
            errors(modelIndex) = inf;
        end
    end
    [~,rateIndex] = min(errors);
    ratePreferred{bandwidthIndex} = models{rateIndex}.name;
end

window = empty_window();
window.requestedLogGrowth = requestedGrowth;
window.actualLogGrowth = y(end)-y(1);
window.startTime = t(1);
window.endTime = t(end);
window.points = originalPointCount;
window.comparisonPoints = numel(t);
window.exponential = models{1};
window.powerLaw = models{2};
window.doubleExponential = models{3};
window.rawPreferredByAicc = models{rawPreferredIndex}.name;
window.preferredByAicc = models{preferredIndex}.name;
window.rawPreferredByValidation = models{rawValidationIndex}.name;
window.preferredByValidation = models{validationIndex}.name;
window.preferredByRateBandwidth = ratePreferred;
window.rateBandwidthAgreementCount = nnz(strcmp( ...
    ratePreferred,models{preferredIndex}.name));
window.rateBandwidthConsensus = ...
    window.rateBandwidthAgreementCount >= ceil(numel(bandwidths)/2);
window.bandwidths = bandwidths;
finiteAicc = sort(selectionAicc(isfinite(selectionAicc)));
if numel(finiteAicc) >= 2
    window.aiccSeparation = finiteAicc(2)-finiteAicc(1);
elseif isscalar(finiteAicc)
    window.aiccSeparation = inf;
end
sameValidationChoice = validationIndex == preferredIndex;
validationImprovement = 0;
if preferredIndex ~= 1
    validationImprovement = 1- ...
        selectionValidation(preferredIndex)/ ...
        max(selectionValidation(1),eps);
end
window.validationImprovementOverExponential = validationImprovement;
if isfinite(selectionAicc(preferredIndex)) && ...
        window.aiccSeparation >= 10 && sameValidationChoice && ...
        window.rateBandwidthConsensus && ...
        (preferredIndex == 1 || validationImprovement >= 0.2)
    window.modelDecision = models{preferredIndex}.name;
else
    window.modelDecision = 'indistinguishable';
end
end

function windows = build_growth_windows( ...
        t,logMaximum,rates,bandwidths,requestedLogGrowth)
totalLogGrowth = logMaximum(end)-logMaximum(1);
windows = repmat(empty_window(),0,1);
for index = 1:numel(requestedLogGrowth)
    target = requestedLogGrowth(index);
    if totalLogGrowth < target
        continue;
    end
    [windowTime,windowLogMaximum,windowRates] = growth_window( ...
        t,logMaximum,rates,target);
    if isempty(windowTime)
        continue;
    end
    window = fit_window(windowTime,windowLogMaximum, ...
        windowRates,bandwidths,target);
    windows(end+1,1) = window; %#ok<AGROW>
end
end

function model = fit_exponential(t,y)
referenceTime = t(1);
design = [ones(numel(t),1),t-referenceTime];
weights = time_quadrature_weights(t);
coefficients = weighted_fit(design,y,weights);
model = base_model('exponential',2,t,y,design*coefficients,weights);
model.referenceTime = referenceTime;
model.parameters = struct('logM0',coefficients(1), ...
    'lambda',coefficients(2));
model.valid = coefficients(2) > 0;
if ~model.valid
    model.aicc = inf;
end
end

function model = fit_power_law(t,y,minimumTerminalTime,searchSpan)
minimumOffset = max(searchSpan/1e4, ...
    100*eps(max(abs(minimumTerminalTime),1)));
maximumOffset = 1e6*searchSpan;
objective = @(logOffset) power_objective( ...
    logOffset,t,y,minimumTerminalTime);
lowerBound = log(minimumOffset);
upperBound = log(maximumOffset);
options = optimset('Display','off','TolX',1e-10);
[logOffset,~,exitFlag] = fminbnd( ...
    objective,lowerBound,upperBound,options);
terminalTime = minimumTerminalTime+exp(logOffset);
coordinate = -log(terminalTime-t);
design = [ones(numel(t),1),coordinate];
weights = time_quadrature_weights(t);
coefficients = weighted_fit(design,y,weights);
model = base_model( ...
    'finite_time_power',3,t,y,design*coefficients,weights);
model.referenceTime = t(1);
model.parameters = struct('logAmplitude',coefficients(1), ...
    'p',coefficients(2),'T',terminalTime);
model.valid = coefficients(2) > 0 && terminalTime > t(end);
model.curvatureAcrossWindow = log( ...
    (terminalTime-t(1))/(terminalTime-t(end)));
model.parameterBounds = [minimumOffset,maximumOffset];
model.profileBoundaryHit = profile_boundary_hit( ...
    logOffset,lowerBound,upperBound);
model.optimizerExitFlag = exitFlag;
model.identifiable = model.valid && ...
    model.curvatureAcrossWindow >= 0.05 && ...
    ~model.profileBoundaryHit;
if ~model.valid
    model.aicc = inf;
end
end

function value = power_objective(logOffset,t,y,minimumTerminalTime)
terminalTime = minimumTerminalTime+exp(logOffset);
coordinate = -log(terminalTime-t);
design = [ones(numel(t),1),coordinate];
weights = time_quadrature_weights(t);
coefficients = weighted_fit(design,y,weights);
if coefficients(2) <= 0
    value = inf;
else
    residual = y-design*coefficients;
    value = sum(weights.*residual.^2);
end
end

function model = fit_double_exponential(t,y,searchSpan)
referenceTime = t(1);
objective = @(logTheta) double_exponential_objective( ...
    logTheta,t,y,referenceTime,searchSpan);
lowerBound = log(1e-6);
upperBound = log(20);
options = optimset('Display','off','TolX',1e-10);
[logTheta,~,exitFlag] = fminbnd( ...
    objective,lowerBound,upperBound,options);
kappa = exp(logTheta)/searchSpan;
coordinate = expm1(kappa*(t-referenceTime))/kappa;
design = [ones(numel(t),1),coordinate];
weights = time_quadrature_weights(t);
coefficients = weighted_fit(design,y,weights);
model = base_model( ...
    'double_exponential',3,t,y,design*coefficients,weights);
model.referenceTime = referenceTime;
model.parameters = struct('logM0',coefficients(1), ...
    'gamma0',coefficients(2),'kappa',kappa);
model.valid = coefficients(2) > 0 && kappa > 0;
model.curvatureAcrossWindow = kappa*(t(end)-t(1));
model.parameterBounds = [exp(lowerBound)/searchSpan, ...
    exp(upperBound)/searchSpan];
model.profileBoundaryHit = profile_boundary_hit( ...
    logTheta,lowerBound,upperBound);
model.optimizerExitFlag = exitFlag;
model.identifiable = model.valid && ...
    model.curvatureAcrossWindow >= 0.05 && ...
    ~model.profileBoundaryHit;
if ~model.valid
    model.aicc = inf;
end
end

function value = double_exponential_objective( ...
        logTheta,t,y,referenceTime,searchSpan)
kappa = exp(logTheta)/searchSpan;
coordinate = expm1(kappa*(t-referenceTime))/kappa;
design = [ones(numel(t),1),coordinate];
weights = time_quadrature_weights(t);
coefficients = weighted_fit(design,y,weights);
if coefficients(2) <= 0
    value = inf;
else
    residual = y-design*coefficients;
    value = sum(weights.*residual.^2);
end
end

function model = base_model(name,parameterCount,t,y,prediction,weights)
residual = y-prediction;
sse = sum(weights.*residual.^2);
points = numel(t);
model = empty_model(name);
model.valid = true;
model.parameterCount = parameterCount;
model.sseLogValue = sse;
model.logValueRmse = weighted_rmse(residual,weights);
model.aicc = points*log(max(sse/points,realmin)) + ...
    2*parameterCount + ...
    2*parameterCount*(parameterCount+1) / ...
    max(points-parameterCount-1,1);
end

function values = model_log_value(model,t)
switch model.name
    case 'exponential'
        values = model.parameters.logM0 + ...
            model.parameters.lambda*(t-model.referenceTime);
    case 'finite_time_power'
        values = model.parameters.logAmplitude - ...
            model.parameters.p*log(model.parameters.T-t);
    case 'double_exponential'
        kappa = model.parameters.kappa;
        coordinate = expm1(kappa*(t-model.referenceTime))/kappa;
        values = model.parameters.logM0 + ...
            model.parameters.gamma0*coordinate;
end
end

function values = model_rate(model,t)
switch model.name
    case 'exponential'
        values = repmat(model.parameters.lambda,size(t));
    case 'finite_time_power'
        values = model.parameters.p./(model.parameters.T-t);
    case 'double_exponential'
        values = model.parameters.gamma0*exp( ...
            model.parameters.kappa*(t-model.referenceTime));
end
end

function errors = rate_log_rmse(rates,prediction,t)
errors = NaN(1,size(rates,2));
for index = 1:size(rates,2)
    valid = isfinite(rates(:,index)) & rates(:,index) > 0 & ...
        isfinite(prediction) & prediction > 0;
    if any(valid)
        residual = log(rates(valid,index))-log(prediction(valid));
        errors(index) = weighted_rmse( ...
            residual,time_quadrature_weights(t(valid)));
    end
end
end

function [rate,rateDerivative] = local_log_rate(t,y,bandwidth)
points = numel(t);
minimumPoints = min(7,points);
rate = NaN(points,1);
rateDerivative = NaN(points,1);
for center = 1:points
    distance = abs(t-t(center));
    selection = find(distance <= bandwidth);
    if numel(selection) < minimumPoints
        [~,order] = sort(distance,'ascend');
        selection = sort(order(1:minimumPoints));
    end
    localScale = max(distance(selection));
    if localScale == 0
        continue;
    end
    coordinate = (t(selection)-t(center))/localScale;
    kernel = max(1-abs(coordinate).^3,0).^3;
    kernel = max(kernel,1e-6);
    design = [ones(numel(selection),1),coordinate,coordinate.^2];
    quadratureWeights = time_quadrature_weights(t(selection));
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

function [windowTime,windowValues,windowRates] = growth_window( ...
        t,values,rates,target)
threshold = values(end)-target;
first = find(values <= threshold,1,'last');
if isempty(first)
    windowTime = [];
    windowValues = [];
    windowRates = [];
    return;
end
latestStart = numel(t)-9;
if first >= latestStart
    first = latestStart;
    windowTime = t(first:end);
    windowValues = values(first:end);
    windowRates = rates(first:end,:);
    return;
end
fraction = (threshold-values(first)) / ...
    (values(first+1)-values(first));
crossingTime = t(first)+fraction*(t(first+1)-t(first));
crossingRates = rates(first,:) + fraction* ...
    (rates(first+1,:)-rates(first,:));
windowTime = [crossingTime;t(first+1:end)];
windowValues = [threshold;values(first+1:end)];
windowRates = [crossingRates;rates(first+1:end,:)];
end

function [comparisonTime,comparisonValues] = ...
        fixed_time_comparison_series(t,values,targetPoints)
if numel(t) <= targetPoints
    comparisonTime = t;
    comparisonValues = values;
    return;
end
comparisonTime = linspace(t(1),t(end),targetPoints)';
comparisonValues = interp1(t,values,comparisonTime,'pchip');
end

function weights = time_quadrature_weights(t)
points = numel(t);
if points <= 1
    weights = ones(size(t));
    return;
end
spacing = diff(t);
weights = [spacing(1)/2; ...
    (spacing(1:end-1)+spacing(2:end))/2;spacing(end)/2];
weights = weights*points/sum(weights);
end

function value = weighted_rmse(residual,weights)
value = sqrt(sum(weights.*residual.^2)/sum(weights));
end

function hit = profile_boundary_hit(value,lowerBound,upperBound)
tolerance = 1e-3*(upperBound-lowerBound);
hit = value-lowerBound <= tolerance || ...
    upperBound-value <= tolerance;
end

function coefficients = weighted_fit(design,values,weights)
weightedDesign = design.*sqrt(weights);
weightedValues = values.*sqrt(weights);
coefficients = weightedDesign\weightedValues;
end

function fit = empty_fit()
fit = struct('valid',false,'points',0,'time',[],'maximum',[], ...
    'logMaximum',[],'logarithmicRate',[], ...
    'logarithmicRateDerivative',[], ...
    'instantaneousTerminalTime',[],'instantaneousPower',[], ...
    'bandwidths',[], ...
    'bandwidthFactors',[],'rateEstimator','', ...
    'modelErrorDomain','','timeWeighting','', ...
    'maximumComparisonPoints',64, ...
    'samplingSensitivityGridPoints',0,'minimumCredibleEFolds',3, ...
    'totalLogGrowth',NaN,'eFoldCount',NaN, ...
    'enoughDynamicRange',false,'finalRate',[], ...
    'finalRateDerivative',[], ...
    'finalInstantaneousTerminalTime',[], ...
    'finalInstantaneousPower',[], ...
    'windows',repmat(empty_window(),0,1), ...
    'modelConsensus',false,'consensusModel','none', ...
    'credibleClassification',false,'credibilityLimitations',{{}});
end

function window = empty_window()
window = struct('requestedLogGrowth',NaN,'actualLogGrowth',NaN, ...
    'startTime',NaN,'endTime',NaN,'points',0,'comparisonPoints',0, ...
    'exponential',empty_model('exponential'), ...
    'powerLaw',empty_model('finite_time_power'), ...
    'doubleExponential',empty_model('double_exponential'), ...
    'rawPreferredByAicc','none','preferredByAicc','none', ...
    'rawPreferredByValidation','none','preferredByValidation','none', ...
    'preferredByRateBandwidth',{{}}, ...
    'rateBandwidthAgreementCount',0, ...
    'rateBandwidthConsensus',false,'bandwidths',[], ...
    'aiccSeparation',NaN, ...
    'validationImprovementOverExponential',NaN, ...
    'preSamplingDecision','indistinguishable', ...
    'uniformInputDecision','unavailable', ...
    'inputSamplingStable',false, ...
    'modelDecision','indistinguishable');
end

function model = empty_model(name)
model = struct('name',name,'valid',false,'parameterCount',0, ...
    'parameters',struct(),'referenceTime',NaN,'sseLogValue',NaN, ...
    'logValueRmse',NaN,'aicc',inf,'rawDeltaAicc',NaN, ...
    'deltaAicc',NaN, ...
    'validationLogRmse',NaN,'rateLogRmse',[], ...
    'identifiable',true,'curvatureAcrossWindow',NaN, ...
    'parameterBounds',[],'profileBoundaryHit',false, ...
    'optimizerExitFlag',NaN);
end
