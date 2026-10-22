function window = outerWallDensityWindow(wallDensity,ops,settings)
%IPM.EVOLVE.OUTERWALLDENSITYWINDOW Fixed X=2 wall-density L2 observable.
% The compact window excludes the tracked singular core near X=1. Its
% normalized quadrature is recomputed on each mesh, while its physical
% coordinate support stays fixed through adaptive remeshing.

if ~isnumeric(wallDensity) || ~isreal(wallDensity) || ...
        ~isvector(wallDensity) || numel(wallDensity) ~= numel(ops.x) || ...
        any(~isfinite(wallDensity),'all')
    error('ipm:OmegaGaugeWindowDensity', ...
        'The outer wall gauge requires finite wall-density samples.');
end
wallDensity = reshape(wallDensity,1,[]);
center = 2;
radius = settings.omegaGaugeWindowRadius;
if ~isfinite(radius) || radius <= 0 || ...
        center-radius <= 0 || center+radius >= ops.x(end)
    error('ipm:OmegaGaugeWindowDomain', ...
        'The X=2 wall-density window must lie inside the positive domain.');
end

coordinate = (ops.x-center)/radius;
support = abs(coordinate) < 1;
bump = zeros(size(ops.x));
bump(support) = cos(pi*coordinate(support)/2).^2;
xQuadrature = sum(ops.integrationWeights,1);
weightedBump = bump.*xQuadrature;
supportPoints = nnz(weightedBump > 0);
if supportPoints < 3 || any(~isfinite(weightedBump),'all') || ...
        any(weightedBump < 0) || sum(weightedBump) <= 0
    error('ipm:OmegaGaugeWindowSupport', ...
        'The X=2 wall-density window needs three positive-weight nodes.');
end
weights = weightedBump/sum(weightedBump);
moment = sum(weights.*wallDensity.^2);
value = sqrt(moment);
condition = value/max(abs(wallDensity));
if ~isfinite(value) || value < settings.cOmegaGaugeFloor || ...
        condition < settings.pointStrainConditionFloor
    error('ipm:DegenerateNormalization', ...
        'The X=2 wall-density window is too small or poorly conditioned.');
end
window = struct('center',center,'radius',radius, ...
    'normalizedWeights',weights,'supportPoints',supportPoints, ...
    'order',2,'moment',moment,'energy',moment, ...
    'value',value,'condition',condition);
end
