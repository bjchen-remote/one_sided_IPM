function window = wallWindowMoment(wallSlope,ops,settings,order)
%IPM.EVOLVE.WALLWINDOWMOMENT Fixed compact wall-slope even moment.

if ~isnumeric(wallSlope) || ~isreal(wallSlope) || ...
        ~isvector(wallSlope) || numel(wallSlope) ~= numel(ops.x) || ...
        any(~isfinite(wallSlope),'all')
    error('ipm:OmegaGaugeWindowSlope', ...
        'The wall-window gauge requires one finite wall-slope sample per x node.');
end
if ~isnumeric(order) || ~isscalar(order) || ~isfinite(order) || ...
        ~any(order == [2,4])
    error('ipm:OmegaGaugeWindowOrder', ...
        'The supported wall-window moment orders are 2 and 4.');
end
wallSlope = reshape(wallSlope,1,[]);
geometry = ipm.evolve.wallWindowGeometry(ops,settings);
weights = geometry.normalizedWeights;
moment = sum(weights.*wallSlope.^order);
energy = sum(weights.*wallSlope.^2);
value = moment^(1/order);
wallScale = max(abs(wallSlope));
condition = value/max(wallScale,realmin);
if ~isfinite(value) || value < settings.cOmegaGaugeFloor || ...
        condition < settings.pointStrainConditionFloor
    error('ipm:DegenerateNormalization', ...
        ['The anchor-wall window gauge has negligible wall-slope ' ...
        'moment relative to the full wall profile.']);
end

window = geometry;
window.order = order;
window.value = value;
window.energy = energy;
window.moment = moment;
window.condition = condition;
end
