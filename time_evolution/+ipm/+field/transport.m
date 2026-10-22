function [rhs,form] = transport( ...
        rho,u1,u2,ops,boundaryMode,closedMassRateTarget)
%IPM.FIELD.TRANSPORT Selectable transport on a nonuniform grid.
%   MUSCL_MINMOD is the robust second-order fallback. WENO5_NONUNIFORM uses
%   fifth-order essentially non-oscillatory point reconstruction on faces with a
%   five-point stencil and MUSCL on the two boundary layers. The physical
%   wall is always impermeable. In the finite-volume paths, OPEN uses
%   boundary nodal states on the three artificial sides and CLOSED sets those
%   fluxes to zero. WENO5_FD is the laboratory mapped finite-difference
%   WENO-Z5 path with a free-stream correction and, for the conservative
%   wall mode, a closed-boundary quadrature projection. WENO7_FD is the
%   paired sixth-order feature path on an admitted smooth mapped grid. Its
%   optional CLOSEDMASSRATETARGET lets the rescaled assembler supply a
%   nonzero weighted target; ordinary transport calls default to zero.
%   The explicit net_velocity_density_centered_advective assembly instead
%   applies one-sided WENO5 density gradients to the actual net velocity.
%   It retains WENO5_FD as the kernel selector and its saved boundary options.

if nargin < 5
    boundaryMode = ops.transportBoundaryMode;
end
if nargin < 6
    closedMassRateTarget = 0;
end
validateattributes(closedMassRateTarget,{'numeric'}, ...
    {'scalar','real','finite'},mfilename,'closedMassRateTarget');
form = struct( ...
    'bulkConservative', ...
    ~ismember(ops.transportScheme,{'weno5_fd','weno7_fd'}), ...
    'advectiveWall',strcmp(ops.wallTransportMode,'advective_upwind'));
if closedMassRateTarget ~= 0 && ...
        ~strcmp(ops.transportScheme,'weno7_fd')
    error('ipm:ClosedMassRateTarget', ...
        'A nonzero closed mass-rate target is supported only by WENO7.');
end
if strcmp(ops.transportScheme,'weno7_fd')
    rhs = weno7_fd_transport( ...
        rho,u1,u2,ops,boundaryMode,closedMassRateTarget);
    return;
elseif strcmp(ops.transportScheme,'weno5_fd')
    if isfield(ops,'rescaling') && isfield(ops.rescaling,'gaugeTransportAssembly') && ...
            strcmp(ops.rescaling.gaugeTransportAssembly,'net_velocity_density_centered_advective')
        rhs = density_advective_weno5(rho,u1,u2,ops,boundaryMode);
    else
        rhs = weno5_fd_transport(rho,u1,u2,ops,boundaryMode);
    end
    return;
end
nx = ops.nx;
ny = ops.ny;

left = zeros(ny,nx);
right = zeros(ny,nx);
deltaX = (rho(:,2:end)-rho(:,1:end-1))./ops.dxFaces;
left(:,2:end) = deltaX;
right(:,1:end-1) = deltaX;
slopeX = minmod(left,right);

down = zeros(ny,nx);
up = zeros(ny,nx);
deltaY = (rho(2:end,:)-rho(1:end-1,:))./ops.dyFaces;
down(2:end,:) = deltaY;
up(1:end-1,:) = deltaY;
slopeY = minmod(down,up);

faceU1 = 0.5*(u1(:,1:end-1)+u1(:,2:end));
rhoLeft = rho(:,1:end-1)+0.5*ops.dxFaces.*slopeX(:,1:end-1);
rhoRight = rho(:,2:end)-0.5*ops.dxFaces.*slopeX(:,2:end);
if strcmp(ops.transportScheme,'weno5_nonuniform')
    faces = ops.wenoX.faces;
    rhoLeft(:,faces) = ipm.field.weno5Reconstruct( ...
        rho,ops.wenoX.left,ops.wenoEpsilon);
    rhoRight(:,faces) = ipm.field.weno5Reconstruct( ...
        rho,ops.wenoX.right,ops.wenoEpsilon);
end
fluxX = zeros(ny,nx+1);
fluxX(:,2:nx) = faceU1 .* ((faceU1 >= 0).*rhoLeft + ...
    (faceU1 < 0).*rhoRight);

faceU2 = 0.5*(u2(1:end-1,:)+u2(2:end,:));
rhoDown = rho(1:end-1,:)+0.5*ops.dyFaces.*slopeY(1:end-1,:);
rhoUp = rho(2:end,:)-0.5*ops.dyFaces.*slopeY(2:end,:);
if strcmp(ops.transportScheme,'weno5_nonuniform')
    faces = ops.wenoY.faces;
    rhoDown(faces,:) = ipm.field.weno5Reconstruct( ...
        rho',ops.wenoY.left,ops.wenoEpsilon)';
    rhoUp(faces,:) = ipm.field.weno5Reconstruct( ...
        rho',ops.wenoY.right,ops.wenoEpsilon)';
end
fluxY = zeros(ny+1,nx);
fluxY(2:ny,:) = faceU2 .* ((faceU2 >= 0).*rhoDown + ...
    (faceU2 < 0).*rhoUp);

if strcmpi(string(boundaryMode),'open')
    fluxX(:,1) = u1(:,1).*rho(:,1);
    fluxX(:,end) = u1(:,end).*rho(:,end);
    fluxY(end,:) = u2(end,:).*rho(end,:);
elseif ~strcmpi(string(boundaryMode),'closed')
    error('ipm:TransportBoundaryMode', ...
        'Transport boundary mode must be ''open'' or ''closed''.');
end
% fluxY(1,:)=0 is the physical no-penetration wall.

divergence = (fluxX(:,2:end)-fluxX(:,1:end-1))./ops.hx + ...
    (fluxY(2:end,:)-fluxY(1:end-1,:))./ops.hy;
rhs = -divergence;

wallMode = lower(char(ops.wallTransportMode));
switch wallMode
    case 'conservative_flux'
        % Keep the bottom half-control-volume update from the common flux
        % divergence. This is the exactly conservative default.
    case 'advective_upwind'
        % On Y=0, impermeability reduces transport to the trace equation
        % rho_t+u_1 rho_x=0. Use the already selected MUSCL/WENO upwind
        % reconstruction on both faces of each bottom-row control volume.
        rhs(1,:) = advective_wall_rhs( ...
            rho(1,:),u1(1,:),rhoLeft(1,:),rhoRight(1,:),ops.hx);
    otherwise
        error('ipm:WallTransportMode', ...
            ['Wall transport mode must be ''conservative_flux'' or ' ...
            '''advective_upwind''.']);
end
end

function rhs = weno7_fd_transport( ...
        rho,u1,u2,ops,boundaryMode,closedMassRateTarget)
% Mapped WENO-Z7 with a discrete free-stream correction. Degree-six ghost
% continuation supplies the sixth-order boundary closure without the larger
% amplification of the unused degree-seven alternative.
options = struct('lowerBoundary','extrapolate', ...
    'upperBoundary','extrapolate','extrapolationDegree',6, ...
    'epsilon',ops.wenoEpsilon);
optionsX = options;
if ipm.mesh.isQuadrant(ops)
    optionsX.lowerBoundary = 'reflect';
end
drhoX = ipm.field.weno7FluxDerivative(rho,u1,ops.metricX, ...
    ops.computationalSpacingX,optionsX);
drhoY = ipm.field.weno7FluxDerivative(rho',u2',ops.metricY', ...
    ops.computationalSpacingY,options)';
unit = ones(size(rho),'like',rho);
divergenceX = ipm.field.weno7FluxDerivative(unit,u1,ops.metricX, ...
    ops.computationalSpacingX,optionsX);
divergenceY = ipm.field.weno7FluxDerivative(unit',u2',ops.metricY', ...
    ops.computationalSpacingY,options)';
rhs = -(drhoX+drhoY-rho.*(divergenceX+divergenceY));

wallMode = lower(char(ops.wallTransportMode));
if strcmp(wallMode,'advective_upwind')
    rhs(1,:) = -(drhoX(1,:)-rho(1,:).*divergenceX(1,:));
elseif ~strcmp(wallMode,'conservative_flux')
    error('ipm:WallTransportMode', ...
        ['Wall transport mode must be ''conservative_flux'' or ' ...
        '''advective_upwind''.']);
end

if strcmpi(string(boundaryMode),'closed') && ...
        strcmp(wallMode,'conservative_flux')
    weights = ops.integrationWeights;
    massRate = sum(rhs.*weights,'all');
    rhs = rhs- ...
        (massRate-closedMassRateTarget)/sum(weights,'all');
elseif ~strcmpi(string(boundaryMode),'open') && ...
        ~strcmpi(string(boundaryMode),'closed')
    error('ipm:TransportBoundaryMode', ...
        'Transport boundary mode must be ''open'' or ''closed''.');
end
end

function rhs=density_advective_weno5(rho,u1,u2,ops,boundaryMode)
% Reconstruct density gradients; the physical velocity never enters the
% nonlinear WENO weights. This explicit assembly is advective throughout.
if ~strcmpi(string(boundaryMode),'open')
    error('ipm:AdvectiveWenoBoundary','The selected density-gradient assembly requires open transport.');
end
rhs=density_advective_dimension(rho,u1,ops.metricX,ops.computationalSpacingX,ops)+ ...
    density_advective_dimension(rho',u2',ops.metricY',ops.computationalSpacingY,ops)';
end

function rhs=density_advective_dimension(rho,velocity,metric,spacing,ops)
% The zero-velocity shortcut also admits the one-row impermeable wall.
if ~isequal(size(velocity),size(rho))
    if isvector(velocity) && numel(velocity)==size(rho,2) && ...
            ~isequal(size(velocity),[size(rho,1),1])
        velocity=reshape(velocity,1,[]);
    end
    velocity=ones(size(rho),'like',rho).*velocity;
    if ~isequal(size(velocity),size(rho))
        error('ipm:AdvectiveWenoVelocitySize','Velocity must broadcast to the density line batch.');
    end
end
rhs=zeros(size(rho),'like',rho);
if ~any(velocity(:)),return;end
lowerInflow=velocity(:,1)>0;upperInflow=velocity(:,end)<0;
options=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate', ...
    'extrapolationDegree',5,'epsilon',ops.wenoEpsilon,'alphaMode',ops.wenoAlphaMode);
if isfield(ops,'wenoInflowDensityDegree')
    options.inflowDensityDegree=ops.wenoInflowDensityDegree;
end
plus=zeros(size(rho),'like',rho);minus=plus;
% Real outflow lines retain the existing five-argument MEX. Only actual
% inflow lines need the MATLAB density-ghost continuation. Both derivative
% directions receive exactly the same physical boundary masks.
groups={find(~lowerInflow & ~upperInflow),find(lowerInflow | upperInflow)};
for group=1:numel(groups)
    lines=groups{group};if isempty(lines),continue;end
    batch=rho(lines,:);unit=ones(size(batch),'like',rho);
    options.physicalLowerInflow=lowerInflow(lines);
    options.physicalUpperInflow=upperInflow(lines);
    plus(lines,:)=ipm.field.weno5FluxDerivative(batch,unit,1,spacing,options);
    minus(lines,:)=-ipm.field.weno5FluxDerivative(batch,-unit,1,spacing,options);
end
rhs=-(max(velocity,0).*plus+min(velocity,0).*minus)./metric;
end

function rhs = weno5_fd_transport(rho,u1,u2,ops,boundaryMode)
% Standard mapped finite-difference WENO-Z fluxes are corrected by their
% constant-state residual.  This converts the discretely imperfect
% divergence form into the equivalent advective form while preserving a
% uniform scalar to roundoff on stretched grids.
options = struct('lowerBoundary','extrapolate', ...
    'upperBoundary','extrapolate','extrapolationDegree',5, ...
    'epsilon',ops.wenoEpsilon);
options.alphaMode = ops.wenoAlphaMode;
if isfield(ops,'wenoInflowDensityDegree')
    options.inflowDensityDegree = ops.wenoInflowDensityDegree;
end
optionsX = options;
if ipm.mesh.isQuadrant(ops)
    optionsX.lowerBoundary = 'reflect';
    optionsX.normalizationNodeCount = 2*ops.nx-1;
    optionsX.symmetricFluxScale = true;
end
drhoX = ipm.field.weno5FluxDerivative(rho,u1,ops.metricX, ...
    ops.computationalSpacingX,optionsX);
drhoY = ipm.field.weno5FluxDerivative(rho',u2',ops.metricY', ...
    ops.computationalSpacingY,options)';
unit = ones(size(rho),'like',rho);
divergenceX = ipm.field.weno5FluxDerivative(unit,u1,ops.metricX, ...
    ops.computationalSpacingX,optionsX);
divergenceY = ipm.field.weno5FluxDerivative(unit',u2',ops.metricY', ...
    ops.computationalSpacingY,options)';
rhs = -(drhoX+drhoY-rho.*(divergenceX+divergenceY));

wallMode = lower(char(ops.wallTransportMode));
if strcmp(wallMode,'advective_upwind')
    rhs(1,:) = -(drhoX(1,:)-rho(1,:).*divergenceX(1,:));
elseif ~strcmp(wallMode,'conservative_flux')
    error('ipm:WallTransportMode', ...
        ['Wall transport mode must be ''conservative_flux'' or ' ...
        '''advective_upwind''.']);
end

if strcmpi(string(boundaryMode),'closed') && ...
        strcmp(wallMode,'conservative_flux')
    weights = ops.integrationWeights;
    massRate = sum(rhs.*weights,'all');
    rhs = rhs-massRate/sum(weights,'all');
elseif ~strcmpi(string(boundaryMode),'open') && ...
        ~strcmpi(string(boundaryMode),'closed')
    error('ipm:TransportBoundaryMode', ...
        'Transport boundary mode must be ''open'' or ''closed''.');
end
end

function slope = minmod(a,b)
slope = 0.5*(sign(a)+sign(b)).*min(abs(a),abs(b));
end

function rhs = advective_wall_rhs(rho,u,rhoLeft,rhoRight,hx)
positiveFaces = [rho(1),rhoLeft,rho(end)];
negativeFaces = [rho(1),rhoRight,rho(end)];
positiveDerivative = (positiveFaces(2:end)-positiveFaces(1:end-1))./hx;
negativeDerivative = (negativeFaces(2:end)-negativeFaces(1:end-1))./hx;
rhoX = (u >= 0).*positiveDerivative+(u < 0).*negativeDerivative;
rhs = -u.*rhoX;
end
