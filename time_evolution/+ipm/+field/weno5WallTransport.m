function [wallRhs,vertical] = weno5WallTransport( ...
        rho,u1,u2,cL,cR,ops,columns,vertical)
%IPM.FIELD.WENO5WALLTRANSPORT Exact open WENO5-FD RHS on selected wall nodes.
%   The horizontal LF split uses the complete wall row, so its alpha
%   (global or local), split-flux scales, and artificial-boundary ghosts
%   are identical to the full two-dimensional transport. The vertical pass
%   uses the full column batch once, matching the full transport's ghost-extrapolation
%   matrix multiplication order; VERTICAL caches the requested wall values
%   for repeated traveling-wave probes.
%   The explicit density-gradient assembly reuses the full transport on its
%   wall row with an exact zero normal velocity and no vertical cache.

if ~strcmp(ops.transportScheme,'weno5_fd') || ...
        ipm.mesh.isQuadrant(ops)
    error('ipm:Weno5WallTransportMode', ...
        'Local wall transport requires full-half-plane WENO5-FD.');
end
if isfield(ops,'rescaling') && isfield(ops.rescaling,'gaugeTransportAssembly') && ...
        ismember(ops.rescaling.gaugeTransportAssembly, ...
        {'net_velocity_density_centered','net_velocity_density_centered_advective'})
    % Match the complete assembler's density constant-mode convention.
    [~,anchorColumn]=min(abs(ops.x));
    rho=rho-rho(1,anchorColumn);
end
columns = reshape(columns,1,[]);
if isempty(columns) || any(columns < 1) || ...
        any(columns > ops.nx) || any(columns ~= floor(columns))
    error('ipm:Weno5WallTransportColumns', ...
        'Wall columns must be valid integer node indices.');
end
if isfield(ops,'rescaling') && isfield(ops.rescaling,'gaugeTransportAssembly') && ...
        strcmp(ops.rescaling.gaugeTransportAssembly,'net_velocity_density_centered_advective')
    wallVelocity=u1(1,:)+cL*ops.x+cR;
    wallAll=ipm.field.transport(rho(1,:),wallVelocity, ...
        zeros(size(wallVelocity),'like',rho),ops,'open');
    wallRhs=wallAll(columns);vertical=[];
    return
end
options = struct('lowerBoundary','extrapolate', ...
    'upperBoundary','extrapolate','extrapolationDegree',5, ...
    'epsilon',ops.wenoEpsilon);
options.alphaMode = ops.wenoAlphaMode;
if isfield(ops,'wenoInflowDensityDegree')
    options.inflowDensityDegree = ops.wenoInflowDensityDegree;
end

if nargin < 8 || isempty(vertical)
    % At the wall, local-alpha WENO uses only faces 0 and 1. Their
    % lower-boundary ghosts depend on the first six physical y nodes.
    % Retain every x line to preserve the original extrapolation order.
    if strcmp(ops.wenoAlphaMode,'local')
        verticalDensity = rho(1:6,:)';
        verticalVelocity = u2(1:6,:)+cL*ops.Y(1:6,:);
        verticalMetric = ops.metricY(1:6)';
    else
        verticalDensity = rho';
        verticalVelocity = u2+cL*ops.Y;
        verticalMetric = ops.metricY';
    end
    densityY = ipm.field.weno5FluxDerivative( ...
        verticalDensity,verticalVelocity',verticalMetric, ...
        ops.computationalSpacingY,options);
    unitY = ones(size(verticalVelocity'),'like',rho);
    divergenceY = ipm.field.weno5FluxDerivative( ...
        unitY,verticalVelocity',verticalMetric, ...
        ops.computationalSpacingY,options);
    vertical = struct('columns',columns, ...
        'density',densityY(columns,1)', ...
        'divergence',divergenceY(columns,1)');
elseif ~isequal(vertical.columns,columns)
    error('ipm:Weno5WallTransportCache', ...
        'The cached vertical wall derivatives belong to other columns.');
end

wallDensity = rho(1,:);
wallVelocity = u1(1,:)+cL*ops.x+cR;
densityX = ipm.field.weno5FluxDerivative( ...
    wallDensity,wallVelocity,ops.metricX, ...
    ops.computationalSpacingX,options);
unitX = ones(size(wallDensity),'like',rho);
divergenceX = ipm.field.weno5FluxDerivative( ...
    unitX,wallVelocity,ops.metricX, ...
    ops.computationalSpacingX,options);

switch ops.wallTransportMode
    case 'conservative_flux'
        wallRhs = -(densityX(columns)+vertical.density- ...
            wallDensity(columns).*(divergenceX(columns)+ ...
            vertical.divergence));
    case 'advective_upwind'
        wallRhs = -(densityX(columns)- ...
            wallDensity(columns).*divergenceX(columns));
    otherwise
        error('ipm:WallTransportMode', ...
            'Unknown wall transport mode ''%s''.',ops.wallTransportMode);
end
end
