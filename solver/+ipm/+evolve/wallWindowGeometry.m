function window = wallWindowGeometry(ops,settings)
%IPM.EVOLVE.WALLWINDOWGEOMETRY Fixed compact wall-window quadrature.

center = settings.transportAnchorX;
radius = settings.omegaGaugeWindowRadius;
anchorTolerance = 100*eps(max(1,abs(center)));
if ~isfinite(center) || abs(center-1) > anchorTolerance
    error('ipm:OmegaGaugeWindowAnchor', ...
        'The fixed anchor-wall window must be centered at X=1.');
end
if ~isfinite(radius) || radius <= 0 || ...
        center-radius <= 0 || center+radius >= ops.x(end)
    error('ipm:OmegaGaugeWindowDomain', ...
        ['The open support of the fixed wall window must lie strictly ' ...
        'in x>0 and inside the positive x boundary.']);
end

coordinate = (ops.x-center)/radius;
support = abs(coordinate) < 1;
bump = zeros(size(ops.x));
bump(support) = exp(1-1./(1-coordinate(support).^2));
xQuadrature = sum(ops.integrationWeights,1);
weightedBump = bump.*xQuadrature;
supportPoints = nnz(weightedBump > 0);
if supportPoints < 3 || any(~isfinite(weightedBump),'all') || ...
        sum(weightedBump) <= 0
    error('ipm:OmegaGaugeWindowSupport', ...
        ['The fixed anchor-wall window must contain at least three ' ...
        'positive-weight x quadrature nodes.']);
end

window = struct('center',center,'radius',radius, ...
    'support',support,'supportPoints',supportPoints, ...
    'normalizedWeights',weightedBump/sum(weightedBump));
end
