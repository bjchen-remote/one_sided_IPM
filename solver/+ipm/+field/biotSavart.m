function [u1,u2,psi,poissonResidual,poissonSolveInfo] = ...
    biotSavart(omega,ops,kappa)
%IPM.FIELD.BIOTSAVART Map omega=d_x rho to u=grad^perp(-Delta)^(-1)omega.

if nargin < 3
    kappa = 1;
end
[psi,interiorRhs,~,poissonSolveInfo,poissonOperator] = ...
    ipm.field.poisson(omega,ops,kappa);
u1 = -(ops.Dy*psi);
u2 = psi*ipm.mesh.oddDx(ops)';

% The physical wall has constant psi=0, hence exact no penetration there.
% Green data on artificial far boundaries must be allowed to carry flux.
u2(1,:) = 0;

highOrder = strcmp(ops.spatialDiscretization,'high_order') || ...
    strcmp(ops.spatialDiscretization,'sixth_order');
if highOrder && isfinite(poissonSolveInfo.relativeResidual)
    % poisson already formed this same unweighted algebraic residual after
    % the direct solve.  Reuse it instead of repeating a large sparse
    % matrix-vector product at every Runge--Kutta stage.
    poissonResidual = poissonSolveInfo.relativeResidual;
else
    psiInterior = psi(2:end-1,2:end-1);
    r = poissonOperator*psiInterior(:)-interiorRhs(:);
    poissonResidual = norm(r,inf)/max(norm(interiorRhs(:),inf),eps);
end
poissonSolveInfo.algebraicResidual = poissonResidual;
if isnan(poissonSolveInfo.relativeResidual)
    poissonSolveInfo.relativeResidual = poissonResidual;
end
end
