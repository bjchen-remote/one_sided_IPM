function [rhs,baseRhs,transportU1,transportU2,conservativeSource] = ...
    assembleRhs(rho,u1,u2,cX,cY,cOmega,cR,ops,boundaryMode)
%IPM.EVOLVE.ASSEMBLERHS Common rescaled transport/source assembly.
%   BASERHS excludes the amplitude source so gauges can use it directly.
%   Explicit closed sixth-order assembly gives WENO7 the weighted spatial
%   source target required by its advective-form mass identity.

if nargin < 9
    boundaryMode = 'open';
end
if cX == 0 && cR == 0
    transportU1 = u1;
else
    transportU1 = u1+cX*ops.X+cR;
end
if cY == 0
    transportU2 = u2;
else
    transportU2 = u2+cY*ops.Y;
end
if cX == cY
    % Preserve the maintained isotropic multiplication and grouping.
    spatialSource = 2*cX;
else
    spatialSource = cX+cY;
end
closedMassRateTarget = 0;
if strcmp(ops.spatialDiscretization,'sixth_order') && ...
        strcmpi(string(boundaryMode),'closed') && ...
        strcmp(ops.wallTransportMode,'conservative_flux')
    closedMassRateTarget = spatialSource* ...
        sum(rho.*ops.integrationWeights,'all');
end
transportDensity=rho;
if isfield(ops,'rescaling') && isfield(ops.rescaling,'gaugeTransportAssembly') && ...
        ismember(ops.rescaling.gaugeTransportAssembly, ...
        {'net_velocity_density_centered','net_velocity_density_centered_advective'})
    % Advective transport is invariant to a spatially constant density.
    % Center only its nonlinear WENO input; the amplitude source below
    % and the velocity/normalization observables retain the original rho.
    [~,anchorColumn]=min(abs(ops.x));
    transportDensity=rho-rho(1,anchorColumn);
end
[transportRhs,transportForm] = ipm.field.transport( ...
    transportDensity,transportU1,transportU2,ops,boundaryMode, ...
    closedMassRateTarget);
if spatialSource == 0 || ~transportForm.bulkConservative
    % Mapped WENO-FD subtracts its own discrete velocity divergence to
    % preserve free streams. It therefore already represents the advective
    % form and must not receive this conversion a second time.
    baseRhs = transportRhs;
else
    baseRhs = transportRhs+spatialSource*rho;
end
if transportForm.advectiveWall
    % The wall trace is already in advective form, so it does not receive
    % the conservative-to-advective spatial-divergence correction.
    baseRhs(1,:) = transportRhs(1,:);
end
if cOmega == 0
    rhs = baseRhs;
else
    rhs = baseRhs+cOmega*rho;
end
conservativeSource = spatialSource+cOmega;
end
