function [u1,u2,psi,source,poissonResidual,poissonSolveInfo] = ...
    velocity(rho,ops,kappa)
%IPM.FIELD.VELOCITY Darcy/Biot-Savart law for IPM on the upper half-plane.
%   -Delta psi = d_x rho and u = grad^perp psi = (-psi_y, psi_x).

if nargin < 3
    kappa = 1;
end
source = rho * ops.Dx';
[u1,u2,psi,poissonResidual,poissonSolveInfo] = ...
    ipm.field.biotSavart(source,ops,kappa);
end
