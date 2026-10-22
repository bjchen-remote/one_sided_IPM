function ops = build(config,gridOverride)
%IPM.MESH.BUILD Grid, derivatives, and Dirichlet -Laplacian.
%   Arrays use rows for y and columns for x. Stretched mode clusters points
%   near x=0,y=0 while retaining a much more distant artificial boundary.

grid = config.grid;
if nargin >= 2 && ~isempty(gridOverride)
    grid = gridOverride;
end
quadrantOnly = isfield(grid,'quadrantOnly') && grid.quadrantOnly;
physics = config.physics;
elliptic = config.elliptic;
transport = config.transport;
scaling = config.scaling;

stretch = grid.gridStretch(:)';
if isscalar(stretch)
    stretch = [stretch,stretch];
end
useCustomGrid = isfield(grid,'customX') && ~isempty(grid.customX);
% A dynamic half-plane profile needs an exact origin node even when the
% left and right far-field distances differ. Keep the historical symmetric
% and even-node mappings byte-for-byte compatible.
centerAsymmetricHalfPlane = ~quadrantOnly && ~useCustomGrid && ...
    strcmp(physics.symmetryMode,'half_plane') && ...
    grid.xlim(1) < 0 && grid.xlim(2) > 0 && mod(grid.nx,2) == 1 && ...
    abs(sum(grid.xlim)) > 100*eps(max(abs(grid.xlim)));
if useCustomGrid
    x = grid.customX(:)';
    if isfield(grid,'customY') && ~isempty(grid.customY)
        y = grid.customY(:);
    else
        error('ipm:CustomGridPair','customX and customY must be supplied together.');
    end
    if quadrantOnly
        sx = linspace(0,1,grid.nx);
    else
        sx = linspace(-1,1,grid.nx);
    end
    sy = linspace(0,1,grid.ny)';
    metricX = ones(size(x));
    metricY = ones(size(y));
elseif strcmpi(string(grid.gridMode),'stretched') && any(stretch > 0)
    if quadrantOnly
        sx = linspace(0,1,grid.nx);
        midpoint = 0;
        halfWidth = grid.xlim(2);
    else
        sx = linspace(-1,1,grid.nx);
        midpoint = mean(grid.xlim);
        halfWidth = diff(grid.xlim)/2;
    end
    [normalizedX,normalizedMetricX] = sinh_axis(sx,stretch(1));
    x = midpoint+halfWidth*normalizedX;
    if centerAsymmetricHalfPlane
        [x,normalizedMetricX] = asymmetric_centered_axis( ...
            sx,grid.xlim,stretch(1));
        halfWidth = 1;
    end
    sy = linspace(0,1,grid.ny)';
    [normalizedY,normalizedMetricY] = sinh_axis(sy,stretch(2));
    y = grid.ymax*normalizedY;
    metricX = halfWidth*normalizedMetricX;
    metricY = grid.ymax*normalizedMetricY;
else
    if quadrantOnly
        sx = linspace(0,1,grid.nx);
    else
        sx = linspace(-1,1,grid.nx);
    end
    sy = linspace(0,1,grid.ny)';
    x = linspace(grid.xlim(1),grid.xlim(2),grid.nx);
    if centerAsymmetricHalfPlane
        [x,metricX] = asymmetric_centered_axis(sx,grid.xlim,0);
    end
    y = linspace(0,grid.ymax,grid.ny)';
    if ~centerAsymmetricHalfPlane
        if quadrantOnly
            metricX = diff(grid.xlim)*ones(size(sx));
        else
            metricX = diff(grid.xlim)/2*ones(size(sx));
        end
    end
    metricY = grid.ymax*ones(size(sy));
end
if quadrantOnly
    assert(x(1)==0 && x(end)>0,'ipm:QuadrantAxis', ...
        'A quadrant mesh must contain x=0 and only nonnegative x.');
end

dxFaces = diff(x);
dyFaces = diff(y);
highOrder = strcmp(transport.spatialDiscretization,'high_order');
sixthOrder = strcmp(transport.spatialDiscretization,'sixth_order');
if sixthOrder
    Dx = ipm.mesh.fdMatrix(x,1,9);
    Dy = ipm.mesh.fdMatrix(y,1,9);
elseif highOrder
    Dx = ipm.mesh.fdMatrix(x,1,7);
    Dy = ipm.mesh.fdMatrix(y,1,7);
elseif useCustomGrid
    Dx = nonuniform_first_derivative(x);
    Dy = nonuniform_first_derivative(y);
else
    Dx = mapped_first_derivative(sx,metricX);
    Dy = mapped_first_derivative(sy,metricY);
end
if quadrantOnly
    if sixthOrder
        stencilWidth = 9;
    elseif highOrder
        stencilWidth = 7;
    else
        stencilWidth = 3;
    end
    Dx = ipm.mesh.quadrantDerivative(x,1,'even',stencilWidth);
    DxOdd = ipm.mesh.quadrantDerivative(x,1,'odd',stencilWidth);
else
    DxOdd = Dx;
end
if quadrantOnly && highOrder && useCustomGrid
    % The logical X map is odd. Its metric is even; a one-sided derivative
    % at the symmetry axis would disagree with the paired full-grid map.
    metricX = (ipm.mesh.quadrantDerivative(sx,1,'odd',7)*x')';
    metricY = ipm.mesh.fdMatrix(sy,1,7)*y;
    if any(metricX <= 0) || any(metricY <= 0)
        error('ipm:HighOrderGridMetric', ...
            'The quadrant high-order mapped metric must remain positive.');
    end
end

% Reject an inadmissible sixth-order axis before assembling or factoring
% the two-dimensional Poisson operator.  Adaptive-remesh retries can visit
% several rejected candidates, so every quality gate intentionally depends
% only on the inexpensive one-dimensional geometry.
if sixthOrder
    if useCustomGrid
        metricX = (ipm.mesh.fdMatrix(sx,1,9)*x')';
        metricY = ipm.mesh.fdMatrix(sy,1,9)*y;
    end
    if any(metricX <= 0) || any(metricY <= 0)
        error('ipm:SixthOrderGridMetric', ...
            'The sixth-order mapped-grid metric must remain positive.');
    end
    integrationHx = ipm.mesh.quadrature(sx,8).*metricX;
    integrationHy = ipm.mesh.quadrature(sy,8).*metricY;
    integrationHx = integrationHx* ...
        ((x(end)-x(1))/sum(integrationHx));
    integrationHy = integrationHy* ...
        ((y(end)-y(1))/sum(integrationHy));
    if any(integrationHx <= 0) || any(integrationHy <= 0)
        error('ipm:SixthOrderQuadratureWeights', ...
            ['The sixth-order grid is too irregular to retain a positive ' ...
            'mapped quadrature norm.']);
    end

    qualityX = ipm.mesh.quality(x,9,integrationHx);
    qualityY = ipm.mesh.quality(y,9,integrationHy);
    if useCustomGrid
        alternateMetricX = (ipm.mesh.fdMatrix(sx,1,7)*x')';
        alternateMetricY = ipm.mesh.fdMatrix(sy,1,7)*y;
        qualityX.metricConsistencyError = max( ...
            abs(metricX-alternateMetricX))/max(metricX);
        qualityY.metricConsistencyError = max( ...
            abs(metricY-alternateMetricY))/max(metricY);
        metricCheck = 'nine_vs_seven_point';
    else
        discreteMetricX = (ipm.mesh.fdMatrix(sx,1,9)*x')';
        discreteMetricY = ipm.mesh.fdMatrix(sy,1,9)*y;
        qualityX.metricConsistencyError = max( ...
            abs(metricX-discreteMetricX))/max(metricX);
        qualityY.metricConsistencyError = max( ...
            abs(metricY-discreteMetricY))/max(metricY);
        metricCheck = 'discrete_vs_exact_map';
    end
    policy = ipm.mesh.sixthOrderPolicy();
    validate_sixth_order_grid(qualityX,'x',policy);
    validate_sixth_order_grid(qualityY,'y',policy);
end

nxi = grid.nx-2;
nyi = grid.ny-2;
if sixthOrder
    if quadrantOnly
        negativeDxx = -ipm.mesh.quadrantDerivative(x,2,'odd',9);
    else
        negativeDxx = -ipm.mesh.fdMatrix(x,2,9);
    end
    negativeDyy = -ipm.mesh.fdMatrix(y,2,9);
    Tx = negativeDxx(2:end-1,2:end-1);
    Ty = negativeDyy(2:end-1,2:end-1);
    xBoundaryCoefficients = negativeDxx(2:end-1,[1,end]);
    yBoundaryCoefficients = negativeDyy(2:end-1,[1,end]);
elseif highOrder
    if quadrantOnly
        negativeDxx = -ipm.mesh.quadrantDerivative(x,2,'odd',7);
    else
        negativeDxx = -ipm.mesh.fdMatrix(x,2,7);
    end
    negativeDyy = -ipm.mesh.fdMatrix(y,2,7);
    Tx = negativeDxx(2:end-1,2:end-1);
    Ty = negativeDyy(2:end-1,2:end-1);
    xBoundaryCoefficients = negativeDxx(2:end-1,[1,end]);
    yBoundaryCoefficients = negativeDyy(2:end-1,[1,end]);
else
    [Tx,xBoundaryCoefficients] = negative_laplacian_1d(x);
    [Ty,yBoundaryCoefficients] = negative_laplacian_1d(y);
end
separablePoisson = (highOrder || sixthOrder) && ...
    strcmp(elliptic.anisotropicPoissonSolver,'separable');
if separablePoisson
    % The two-dimensional operator is exactly a Kronecker sum.  Retain its
    % one-dimensional factors, but never assemble or factor the large sum.
    A = [];
else
    A = kron(speye(nxi),Ty)+kron(Tx,speye(nyi));
end

hx = zeros(1,grid.nx);
hx([1,end]) = [dxFaces(1),dxFaces(end)]/2;
hx(2:end-1) = (x(3:end)-x(1:end-2))/2;
hy = zeros(grid.ny,1);
hy([1,end]) = [dyFaces(1),dyFaces(end)]/2;
hy(2:end-1) = (y(3:end)-y(1:end-2))/2;

ops.x = x;
ops.y = y;
ops.baseX = x;
ops.baseY = y;
ops.dx = min(dxFaces);
ops.dy = min(dyFaces);
ops.dxFaces = dxFaces;
ops.dyFaces = dyFaces;
ops.X = repmat(x,grid.ny,1);
ops.Y = repmat(y,1,grid.nx);
ops.Dx = Dx;
if quadrantOnly
    ops.DxOdd = DxOdd;
end
ops.Dy = Dy;
ops.A = A;
ops.Tx = Tx;
ops.Ty = Ty;
ops.hx = hx;
ops.hy = hy;
ops.weights = hy*hx;
ops.integrationWeights = ops.weights;
if highOrder
    if quadrantOnly
        integrationHx = ipm.mesh.quadrantQuadrature(x);
    else
        integrationHx = ipm.mesh.quadrature(x);
    end
    integrationHy = ipm.mesh.quadrature(y);
    if any(integrationHx <= 0) || any(integrationHy <= 0)
        error('ipm:HighOrderQuadratureWeights', ...
            ['The high-order grid is too irregular to retain a positive ' ...
            'quadrature norm.']);
    end
    ops.integrationWeights = integrationHy*integrationHx;
    if quadrantOnly
        ops.integrationHx = integrationHx;
        ops.integrationHy = integrationHy;
    end
    if separablePoisson
        ops.poisson = [];
        [Qy,Sy] = schur(full(Ty),'real');
        [Qx,Sx] = schur(full(Tx.'),'real');
        ops.separablePoisson = struct('Qy',Qy,'Sy',Sy, ...
            'Qx',Qx,'Sx',Sx);
    else
        ops.poisson = decomposition(A,'lu');
    end
elseif sixthOrder
    ops.integrationWeights = integrationHy*integrationHx;
    ops.integrationHx = integrationHx;
    ops.integrationHy = integrationHy;
    if separablePoisson
        ops.poisson = [];
        [Qy,Sy] = schur(full(Ty),'real');
        [Qx,Sx] = schur(full(Tx.'),'real');
        ops.separablePoisson = struct('Qy',Qy,'Sy',Sy, ...
            'Qx',Qx,'Sx',Sx);
    else
        ops.poisson = decomposition(A,'lu');
    end
else
    poissonWeights = kron(hx(2:end-1)',hy(2:end-1));
    weightedA = spdiags(poissonWeights,0,numel(poissonWeights), ...
        numel(poissonWeights))*A;
    weightedA = (weightedA+weightedA')/2;
    ops.poisson = decomposition(weightedA,'chol');
    ops.poissonWeights = poissonWeights;
end
ops.poissonBoundary = struct('x',xBoundaryCoefficients, ...
    'y',yBoundaryCoefficients);
ops.farBoundaryMode = lower(char(elliptic.farBoundaryMode));
ops.transportBoundaryMode = lower(char(transport.transportBoundaryMode));
ops.transportScheme = lower(char(transport.transportScheme));
ops.spatialDiscretization = lower(char(transport.spatialDiscretization));
ops.wallTransportMode = lower(char(transport.wallTransportMode));
ops.wenoEpsilon = transport.wenoEpsilon;
if isfield(transport,'wenoAlphaMode')
    ops.wenoAlphaMode = lower(char(transport.wenoAlphaMode));
else
    ops.wenoAlphaMode = 'global';
end
if isfield(transport,'wenoInflowDensityDegree')
    ops.wenoInflowDensityDegree = transport.wenoInflowDensityDegree;
end
if strcmp(ops.transportScheme,'weno5_nonuniform')
    ops.wenoX = ipm.field.weno5Geometry(x);
    ops.wenoY = ipm.field.weno5Geometry(y');
end
if highOrder
    if useCustomGrid && ~quadrantOnly
        metricX = (ipm.mesh.fdMatrix(sx,1,7)*x')';
        metricY = ipm.mesh.fdMatrix(sy,1,7)*y;
    end
    if any(metricX <= 0) || any(metricY <= 0)
        error('ipm:HighOrderGridMetric', ...
            'The high-order mapped-grid metric must remain positive.');
    end
    ops.computationalSpacingX = sx(2)-sx(1);
    ops.computationalSpacingY = sy(2)-sy(1);
    ops.metricX = metricX;
    ops.metricY = metricY;
end
if sixthOrder
    computationalSpacingX = sx(2)-sx(1);
    computationalSpacingY = sy(2)-sy(1);
    nearestFaceX = [dxFaces(1), ...
        min(dxFaces(1:end-1),dxFaces(2:end)),dxFaces(end)];
    nearestFaceY = [dyFaces(1); ...
        min(dyFaces(1:end-1),dyFaces(2:end));dyFaces(end)];
    transportSpacingX = min( ...
        metricX*computationalSpacingX,nearestFaceX);
    transportSpacingY = min( ...
        metricY*computationalSpacingY,nearestFaceY);

    ops.computationalSpacingX = computationalSpacingX;
    ops.computationalSpacingY = computationalSpacingY;
    ops.metricX = metricX;
    ops.metricY = metricY;
    ops.transportSpacingX = transportSpacingX;
    ops.transportSpacingY = transportSpacingY;
    ops.derivativeStencilWidth = 9;
    ops.quadratureStencilWidth = 8;
    ops.gridQuality = struct('x',qualityX,'y',qualityY, ...
        'metricCheck',metricCheck, ...
        'maximumCellRatio',policy.maximumAdjacentCellRatio, ...
        'minimumStencilRcond',policy.minimumStencilRcond, ...
        'minimumQuadratureWeightRatio', ...
        policy.minimumQuadratureWeightRatio, ...
        'maximumMetricError',policy.maximumMetricError);
end
ops.symmetryMode = lower(char(physics.symmetryMode));
if quadrantOnly
    ops.quadrantOnly = true;
end
ops.greenSourceTolerance = elliptic.greenSourceTolerance;
ops.greenMaxSources = elliptic.greenMaxSources;
ops.nx = grid.nx;
ops.ny = grid.ny;
ops.remeshCount = 0;
ops.rescalingMode = lower(char(scaling.rescalingMode));
ops.dynamicScaleGeometry = lower(char(scaling.dynamicScaleGeometry));
ops.anisotropicPoisson = struct( ...
    'solver',lower(char(elliptic.anisotropicPoissonSolver)), ...
    'tolerance',elliptic.anisotropicPoissonTolerance, ...
    'maxIterations',elliptic.anisotropicPoissonMaxIterations);
ops.rescaling = struct('enabled',false, ...
    'scalingContract',lower(char(scaling.scalingContract)), ...
    'anisotropicGaugeMode',lower(char(scaling.anisotropicGaugeMode)), ...
    'fixedCX',scaling.anisotropicFixedCX, ...
    'fixedCY',scaling.anisotropicFixedCY, ...
    'fixedCOmega',scaling.anisotropicFixedCOmega, ...
    'fixedCR',scaling.anisotropicFixedCR, ...
    'cOmegaGauge',lower(char(scaling.cOmegaGauge)), ...
    'cOmegaGaugeFloor',scaling.cOmegaGaugeFloor, ...
    'pointStrainConditionFloor',scaling.pointStrainConditionFloor, ...
    'omegaGaugeWindowRadius',scaling.omegaGaugeWindowRadius, ...
    'lengthGauge',lower(char(scaling.lengthGauge)), ...
    'transportAnchorX',scaling.transportAnchorX, ...
    'lengthGaugeLevels',scaling.lengthGaugeLevels(:)', ...
    'targetPeakPoints',scaling.targetPeakPoints, ...
    'adaptiveLevels',scaling.adaptiveLevels(:)', ...
    'targetLevelPoints',scaling.targetLevelPoints(:)', ...
    'targetPeakExpansion',scaling.targetPeakExpansion, ...
    'minimumPeakPoints',scaling.minimumPeakPoints, ...
    'minimumLevelPoints',scaling.minimumLevelPoints(:)', ...
    'peakTrackingWidthFactor',scaling.peakTrackingWidthFactor);
if isfield(scaling,'peakPositionRestoringRate')
    ops.rescaling.peakPositionRestoringRate = ...
        scaling.peakPositionRestoringRate;
else
    ops.rescaling.peakPositionRestoringRate = 0;
end
if isfield(scaling,'peakCurvatureRadius')
    ops.rescaling.peakCurvatureRadius = ...
        scaling.peakCurvatureRadius;
else
    ops.rescaling.peakCurvatureRadius = 0.75;
end
ops.rescaling.travelingWaveGauge = lower(char(scaling.travelingWaveGauge));
end

function validate_sixth_order_grid(quality,axisName,policy)
if quality.maximumAdjacentCellRatio > ...
        policy.maximumAdjacentCellRatio
    error('ipm:SixthOrderGridCellRatio', ...
        ['The sixth-order %s grid has adjacent-cell ratio %.6g, ' ...
        'exceeding the policy limit %.6g.'],axisName, ...
        quality.maximumAdjacentCellRatio, ...
        policy.maximumAdjacentCellRatio);
end
if quality.minimumStencilRcond < policy.minimumStencilRcond
    error('ipm:SixthOrderGridStencilCondition', ...
        ['The sixth-order %s grid has minimum normalized stencil ' ...
        'rcond %.3e, below the policy floor %.3e.'],axisName, ...
        quality.minimumStencilRcond,policy.minimumStencilRcond);
end
if quality.minimumQuadratureWeightRatio < ...
        policy.minimumQuadratureWeightRatio
    error('ipm:SixthOrderGridQuadratureMargin', ...
        ['The sixth-order %s grid has normalized quadrature-weight ' ...
        'margin %.3e, below the policy floor %.3e.'],axisName, ...
        quality.minimumQuadratureWeightRatio, ...
        policy.minimumQuadratureWeightRatio);
end
if quality.metricConsistencyError > policy.maximumMetricError
    error('ipm:SixthOrderGridMetricConsistency', ...
        ['The sixth-order %s grid has metric consistency error %.3e, ' ...
        'exceeding the policy limit %.3e.'],axisName, ...
        quality.metricConsistencyError,policy.maximumMetricError);
end
end

function [mapped,metric] = sinh_axis(reference,stretch)
if stretch == 0
    mapped = reference;
    metric = ones(size(reference));
else
    mapped = sinh(stretch*reference)/sinh(stretch);
    metric = stretch*cosh(stretch*reference)/sinh(stretch);
end
end

function [axis,metric] = asymmetric_centered_axis(reference,limits,stretch)
% Smooth asymmetric deformation of the old centered sinh map. The even
% term changes the two outer distances while preserving x(0)=0 and one
% common core spacing on the negative and positive sides.
left = -limits(1);
right = limits(2);
halfSum = (left+right)/2;
halfDifference = (right-left)/2;
if stretch == 0
    oddMap = reference;
    oddMetric = ones(size(reference));
    evenMap = reference.^2;
    evenMetric = 2*reference;
else
    [oddMap,oddMetric] = sinh_axis(reference,stretch);
    halfSinh = sinh(stretch/2);
    evenMap = (sinh(stretch*reference/2)/halfSinh).^2;
    evenMetric = stretch*sinh(stretch*reference)/(2*halfSinh^2);
end
axis = halfSum*oddMap+halfDifference*evenMap;
metric = halfSum*oddMetric+halfDifference*evenMetric;
axis([1,end]) = limits;
axis((numel(axis)+1)/2) = 0;
if any(~isfinite(metric)) || any(metric <= 0) || ...
        any(~isfinite(axis)) || any(diff(axis) <= 0)
    error('ipm:HalfPlaneGridMetric', ...
        ['The asymmetric half-plane x grid is not monotone. Increase ' ...
        'gridStretch or reduce the left/right domain imbalance.']);
end
end

function D = mapped_first_derivative(reference,metric)
reference = reference(:);
metric = metric(:);
n = numel(reference);
h = reference(2)-reference(1);
e = ones(n,1);
Dref = spdiags([-e,e],[-1,1],n,n)/(2*h);
Dref(1,:) = 0;
Dref(1,1:3) = [-3,4,-1]/(2*h);
Dref(end,:) = 0;
Dref(end,end-2:end) = [1,-4,3]/(2*h);
D = spdiags(1./metric,0,n,n)*Dref;
end

function [L,boundaryCoefficients] = negative_laplacian_1d(z)
z = z(:);
nInterior = numel(z)-2;
interior = (2:numel(z)-1)';
hLeft = z(interior)-z(interior-1);
hRight = z(interior+1)-z(interior);
left = -2./(hLeft.*(hLeft+hRight));
center = 2./(hLeft.*hRight);
right = -2./(hRight.*(hLeft+hRight));
rows = [(2:nInterior)';(1:nInterior)';(1:nInterior-1)'];
cols = [(1:nInterior-1)';(1:nInterior)';(2:nInterior)'];
values = [left(2:end);center;right(1:end-1)];
L = sparse(rows,cols,values,nInterior,nInterior);
boundaryCoefficients = [left(1),right(end)];
end

function D = nonuniform_first_derivative(z)
z = z(:);
n = numel(z);
rows = zeros(3*n,1);
cols = zeros(3*n,1);
values = zeros(3*n,1);
cursor = 0;
for index = 1:n
    if index == 1
        stencil = 1:3;
    elseif index == n
        stencil = n-2:n;
    else
        stencil = index-1:index+1;
    end
    nodes = z(stencil)-z(index);
    vandermonde = [ones(1,3);nodes';(nodes.^2)'];
    coefficients = vandermonde\[0;1;0];
    slots = cursor+(1:3);
    rows(slots) = index;
    cols(slots) = stencil;
    values(slots) = coefficients;
    cursor = cursor+3;
end
D = sparse(rows,cols,values,n,n);
end
