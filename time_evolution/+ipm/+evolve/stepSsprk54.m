function [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
    stepSsprk54(rho,dt,ops,scale,initialRhsCache)
%IPM.EVOLVE.STEPSSPRK54 Five-stage fourth-order time-step entry point.

if nargin < 4
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'ssprk54');
elseif nargin < 5
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'ssprk54',scale);
else
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'ssprk54',scale,initialRhsCache);
end
end
