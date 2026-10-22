function [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
    stepSsprk3(rho,dt,ops,scale,initialRhsCache)
%IPM.EVOLVE.STEPSSPRK3 Maintained third-order time-step entry point.

if nargin < 4
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'ssprk3');
elseif nargin < 5
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'ssprk3',scale);
else
    [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
        ipm.evolve.stepRk(rho,dt,ops,'ssprk3',scale,initialRhsCache);
end
end
