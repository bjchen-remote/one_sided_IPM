function [rhoNew,flow,scaleNew,rhsCache,evaluationInfo] = ...
    stepRk(rho,dt,ops,method,scale,initialRhsCache)
%IPM.EVOLVE.STEPRK Shared stage-state contract for the supported RK methods.
%   Each method retains its explicit arithmetic expression. Only scale
%   packing, stage RHS evaluation, and final reconstruction are shared. An
%   optional exact cache may supply only the first stage at the identical
%   accepted field, grid, and scale.

isAnisotropic = strcmp(ops.dynamicScaleGeometry,'anisotropic');
if nargin < 5
    if isAnisotropic
        error('ipm:MissingAnisotropicScale', ...
            'Anisotropic %s requires the current scale state.', ...
            upper(method));
    end
    scale = ipm.evolve.initialScale(ops);
end
if nargin < 6
    initialRhsCache = [];
end
if ~any(strcmp(method,{'ssprk3','ssprk54','rk6'}))
    error('ipm:TimeIntegrator','Unknown RK method ''%s''.',method);
end

z = pack_scale(scale,isAnisotropic);
[cacheHit,k1,s1] = cached_first_stage( ...
    initialRhsCache,rho,z,ops,isAnisotropic);
if ~cacheHit
    [k1,s1] = stage_rhs(rho,z,ops,isAnisotropic);
end
switch method
    case 'ssprk3'
        stageCount = 3;
        q1 = rho+dt*k1;
        z1 = z+dt*s1;

        [k2,s2] = stage_rhs(q1,z1,ops,isAnisotropic);
        q2 = 0.75*rho+0.25*(q1+dt*k2);
        z2 = 0.75*z+0.25*(z1+dt*s2);

        [k3,s3] = stage_rhs(q2,z2,ops,isAnisotropic);
        rhoNew = (rho+2*(q2+dt*k3))/3;
        zNew = (z+2*(z2+dt*s3))/3;
    case 'ssprk54'
        stageCount = 5;
        tableau = ipm.evolve.ssprk54Tableau();
        A = tableau.A;
        b = tableau.b;

        q2 = rho+dt*A(2,1)*k1;
        z2 = z+dt*A(2,1)*s1;
        [k2,s2] = stage_rhs(q2,z2,ops,isAnisotropic);

        q3 = rho+dt*(A(3,1)*k1+A(3,2)*k2);
        z3 = z+dt*(A(3,1)*s1+A(3,2)*s2);
        [k3,s3] = stage_rhs(q3,z3,ops,isAnisotropic);

        q4 = rho+dt*(A(4,1)*k1+A(4,2)*k2+A(4,3)*k3);
        z4 = z+dt*(A(4,1)*s1+A(4,2)*s2+A(4,3)*s3);
        [k4,s4] = stage_rhs(q4,z4,ops,isAnisotropic);

        q5 = rho+dt*(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3+A(5,4)*k4);
        z5 = z+dt*(A(5,1)*s1+A(5,2)*s2+A(5,3)*s3+A(5,4)*s4);
        [k5,s5] = stage_rhs(q5,z5,ops,isAnisotropic);

        rhoNew = rho+dt*( ...
            b(1)*k1+b(2)*k2+b(3)*k3+b(4)*k4+b(5)*k5);
        zNew = z+dt*( ...
            b(1)*s1+b(2)*s2+b(3)*s3+b(4)*s4+b(5)*s5);
    case 'rk6'
        tableau = ipm.evolve.rk6Tableau();
        stageCount = tableau.stages;
        densityRates = cell(tableau.stages,1);
        scaleRates = cell(tableau.stages,1);
        densityRates{1} = k1;
        scaleRates{1} = s1;
        for stage = 2:tableau.stages
            stageDensity = rho;
            stageScale = z;
            for previous = 1:stage-1
                coefficient = tableau.A(stage,previous);
                if coefficient ~= 0
                    coefficient = dt*coefficient;
                    stageDensity = stageDensity + ...
                        coefficient*densityRates{previous};
                    stageScale = stageScale + ...
                        coefficient*scaleRates{previous};
                end
            end
            [densityRates{stage},scaleRates{stage}] = stage_rhs( ...
                stageDensity,stageScale,ops,isAnisotropic);
        end
        rhoNew = rho;
        zNew = z;
        for stage = 1:tableau.stages
            coefficient = tableau.b(stage);
            if coefficient ~= 0
                coefficient = dt*coefficient;
                rhoNew = rhoNew+coefficient*densityRates{stage};
                zNew = zNew+coefficient*scaleRates{stage};
            end
        end
end
scaleNew = unpack_scale(zNew,isAnisotropic);
% Canonical tau is the authoritative integration clock.  Reconstruct its
% constant-rate update directly instead of accumulating roundoff through
% the Runge--Kutta stage recombinations.
scaleNew.canonicalTime = scale.canonicalTime+dt;
[terminalRhoRate,flow] = ipm.evolve.flow(rhoNew,ops,scaleNew);
rhsCache = ipm.evolve.makeRhsCache( ...
    rhoNew,terminalRhoRate,flow,ops,scaleNew);
evaluationInfo = struct('stageCount',stageCount, ...
    'rhsEvaluations',stageCount+1-double(cacheHit), ...
    'initialStageCacheHit',cacheHit,'terminalRhsEvaluated',true);
end

function [hit,rhoRate,scaleRate] = ...
        cached_first_stage(cache,rho,z,ops,isAnisotropic)
hit = false;
rhoRate = [];
scaleRate = [];
required = {'schemaVersion','geometry','remeshCount','x','y','rho', ...
    'scale','rhoRate','scaleRate'};
if ~isstruct(cache) || ~isscalar(cache) || ...
        ~all(isfield(cache,required)) || ...
        ~(isnumeric(cache.schemaVersion) && isscalar(cache.schemaVersion) && ...
        isfinite(cache.schemaVersion) && cache.schemaVersion == 1) || ...
        ~strcmp(cache.geometry,ops.dynamicScaleGeometry) || ...
        ~isequal(cache.remeshCount,ops.remeshCount) || ...
        ~isequal(cache.x,ops.x) || ~isequal(cache.y,ops.y) || ...
        ~isequaln(cache.rho,rho) || ...
        ~isequaln(cache.scale,unpack_scale(z,isAnisotropic)) || ...
        ~isequal(size(cache.rhoRate),size(rho)) || ...
        ~isnumeric(cache.rhoRate) || ~isreal(cache.rhoRate) || ...
        any(~isfinite(cache.rhoRate),'all') || ...
        ~isnumeric(cache.scaleRate) || ~isreal(cache.scaleRate) || ...
        ~isequal(size(cache.scaleRate),size(z)) || ...
        any(~isfinite(cache.scaleRate),'all')
    return;
end
hit = true;
rhoRate = cache.rhoRate;
scaleRate = cache.scaleRate;
end

function [rhoRate,scaleRate] = stage_rhs(rho,z,ops,isAnisotropic)
scale = unpack_scale(z,isAnisotropic);
[rhoRate,flow] = ipm.evolve.flow(rho,ops,scale);
if isAnisotropic
    scaleRate = [flow.c_x;flow.c_y;flow.c_omega; ...
        exp(z(3)-z(1)); ...
        flow.c_x*z(5)+flow.c_r;1];
else
    scaleRate = [flow.c_l;flow.c_omega; ...
        exp(z(2)-z(1)); ...
        flow.c_l*z(4)+flow.c_r;1];
end
end

function z = pack_scale(scale,isAnisotropic)
if isAnisotropic
    z = [scale.logC_l;scale.logC_y;scale.logC_omega; ...
        scale.physicalTime;scale.X_shift;scale.canonicalTime];
else
    z = [scale.logC_l;scale.logC_omega;scale.physicalTime; ...
        scale.X_shift;scale.canonicalTime];
end
end

function scale = unpack_scale(z,isAnisotropic)
if isAnisotropic
    scale = struct('logC_l',z(1),'logC_x',z(1),'logC_y',z(2), ...
        'logC_omega',z(3),'physicalTime',z(4),'X_shift',z(5), ...
        'canonicalTime',z(6));
else
    scale = struct('logC_l',z(1),'logC_omega',z(2), ...
        'physicalTime',z(3),'X_shift',z(4),'canonicalTime',z(5));
end
end
