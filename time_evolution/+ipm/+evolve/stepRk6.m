function [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
    stepRk6(rho,dt,ops,scale,initialRhsCache)
%IPM.EVOLVE.STEPRK6 Eight-stage sixth-order time-step entry point.

if nargin < 4
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'rk6');
elseif nargin < 5
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'rk6',scale);
else
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'rk6',scale,initialRhsCache);
end
end
