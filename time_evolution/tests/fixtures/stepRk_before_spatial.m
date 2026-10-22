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
if ~any(strcmp(method,{'ssprk3','ssprk54','rk6','sdirk3_alexander_practical_v1', ...
        'sdirk3_alexander_practical_v2','sdirk3_alexander_practical_v3','sdirk2_alexander_practical_v1','sdirk2_alexander_practical_v2'}))
    error('ipm:TimeIntegrator','Unknown RK method ''%s''.',method);
end

z = pack_scale(scale,isAnisotropic);
[cacheHit,k1,s1] = cached_first_stage( ...
    initialRhsCache,rho,z,ops,isAnisotropic);
if ~cacheHit
    [k1,s1,ops] = stage_rhs(rho,z,ops,isAnisotropic);
elseif ~isAnisotropic && ops.rescaling.enabled && ...
        strcmp(ops.rescaling.lengthGauge,'wall_gradient_integral')
    if isfield(initialRhsCache,'gaugeRates') && ...
            isnumeric(initialRhsCache.gaugeRates) && ...
            isreal(initialRhsCache.gaugeRates) && ...
            numel(initialRhsCache.gaugeRates) == 2 && ...
            all(isfinite(initialRhsCache.gaugeRates))
        ops.rescaling.gaugeRateGuess = initialRhsCache.gaugeRates(:)';
    else
        ops.rescaling.gaugeRateGuess = [s1(1),s1(4)-s1(1)*z(4)];
    end
end
switch method
    case 'ssprk3'
        stageCount = 3;
        q1 = rho+dt*k1;
        z1 = z+dt*s1;

        [k2,s2,ops] = stage_rhs(q1,z1,ops,isAnisotropic);
        q2 = 0.75*rho+0.25*(q1+dt*k2);
        z2 = 0.75*z+0.25*(z1+dt*s2);

        [k3,s3,ops] = stage_rhs(q2,z2,ops,isAnisotropic);
        rhoNew = (rho+2*(q2+dt*k3))/3;
        zNew = (z+2*(z2+dt*s3))/3;
    case 'ssprk54'
        stageCount = 5;
        tableau = ipm.evolve.ssprk54Tableau();
        A = tableau.A;
        b = tableau.b;

        q2 = rho+dt*A(2,1)*k1;
        z2 = z+dt*A(2,1)*s1;
        [k2,s2,ops] = stage_rhs(q2,z2,ops,isAnisotropic);

        q3 = rho+dt*(A(3,1)*k1+A(3,2)*k2);
        z3 = z+dt*(A(3,1)*s1+A(3,2)*s2);
        [k3,s3,ops] = stage_rhs(q3,z3,ops,isAnisotropic);

        q4 = rho+dt*(A(4,1)*k1+A(4,2)*k2+A(4,3)*k3);
        z4 = z+dt*(A(4,1)*s1+A(4,2)*s2+A(4,3)*s3);
        [k4,s4,ops] = stage_rhs(q4,z4,ops,isAnisotropic);

        q5 = rho+dt*(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3+A(5,4)*k4);
        z5 = z+dt*(A(5,1)*s1+A(5,2)*s2+A(5,3)*s3+A(5,4)*s4);
        [k5,s5,ops] = stage_rhs(q5,z5,ops,isAnisotropic);

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
            [densityRates{stage},scaleRates{stage},ops] = stage_rhs( ...
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
    case {'sdirk3_alexander_practical_v1','sdirk3_alexander_practical_v2','sdirk3_alexander_practical_v3','sdirk2_alexander_practical_v1','sdirk2_alexander_practical_v2'}
        assert(~isAnisotropic && isfield(ops,'timeIntegrationContract'), ...
            'ipm:Sdirk3NumericalTuple','The implicit method requires its versioned isotropic contract.');
        [rhoNew,zNew,implicitInfo,ops]=sdirk_step(rho,z,dt,ops,k1);
        stageCount=numel(ops.timeIntegrationContract.c);
end
scaleNew = unpack_scale(zNew,isAnisotropic);
% Canonical tau is the authoritative integration clock.  Reconstruct its
% constant-rate update directly instead of accumulating roundoff through
% the Runge--Kutta stage recombinations.
scaleNew.canonicalTime = scale.canonicalTime+dt;
if strcmp(ops.rescaling.cOmegaGauge,'wall_omega_tracked_spline_peak')
    ops.rescaling.peakEvaluationTime=scaleNew.canonicalTime;
end
[terminalRhoRate,flow] = ipm.evolve.flow(rhoNew,ops,scaleNew);
rhsCache = ipm.evolve.makeRhsCache( ...
    rhoNew,terminalRhoRate,flow,ops,scaleNew);
evaluationInfo = struct('stageCount',stageCount, ...
    'rhsEvaluations',stageCount+1-double(cacheHit), ...
    'initialStageCacheHit',cacheHit,'terminalRhsEvaluated',true);
if exist('implicitInfo','var')
    evaluationInfo.rhsEvaluations=implicitInfo.rhsEvaluations+1+double(~cacheHit);
    evaluationInfo.implicit=implicitInfo;
end
end

function [rhoNew,zNew,info,ops]=sdirk_step(rho,z,dt,ops,initialRate)
% Full native stage Newton; approximate transport appears only in the ILU.
% Source20 isotropic F is independent of scale. The triangular scale ODE is
% eliminated exactly with the same tableau once each density stage resolves.
c=ops.timeIntegrationContract;A=c.A;gamma=c.gamma;
stageCount=numel(c.c);residualPower=4;
if isfield(c,'stageResidualPower'),residualPower=c.stageResidualPower;end
densityTolerance=max(c.densityToleranceFloor, ...
    min(c.densityToleranceCeiling,c.densityToleranceCoefficient*dt^residualPower));
rxTolerance=max(c.rxToleranceFloor, ...
    min(c.rxToleranceCeiling,c.rxToleranceCoefficient*dt^residualPower));
rates=cell(1,stageCount);scaleRates=cell(1,stageCount);rhoStages=cell(1,stageCount);zStages=cell(1,stageCount);
evaluations=0;records=struct([]);stageIndex=0;stageOps=ops;
for stageIndex=1:stageCount
    before=evaluations;known=rho;knownScale=z;
    for previous=1:stageIndex-1
        known=known+dt*A(stageIndex,previous)*rates{previous};
        knownScale=knownScale+dt*A(stageIndex,previous)*scaleRates{previous};
    end
    if stageIndex==1,predictorRate=initialRate;else,predictorRate=rates{stageIndex-1};end
    trialRho=known+gamma*dt*predictorRate;
    variableScale=max(1,max(abs(rho),abs(known)));
    stageOps=ops;found=false;lastFlag=NaN;lastRelres=NaN;lastIterations=[0,0];lastLength=NaN;
    [rate,flow]=native(trialRho);residual=trialRho-known-gamma*dt*rate;
    [merit,densityNorm,rxNorm]=measure(residual);
    updates=0;cheapTrials=0;cheapAccepted=0;
    for newton=1:c.maximumNewtonIterations
        if merit<=1,found=true;break;end
        [L,U,rowBalance]=preconditioner(flow,variableScale);
        baseRho=trialRho;
        right=@(v)U\(L\v);
        accepted=false;
        if isfield(c,'preconditionerDefectTrials') && cheapTrials<c.preconditionerDefectTrials
            cheapTrials=cheapTrials+1;
            cheapScaledStep=right(-rowBalance.*residual(:));
            cheapStep=reshape(variableScale(:).*cheapScaledStep,size(rho));
            candidate=trialRho+cheapStep;candidateValid=all(isfinite(candidate),'all');
            if candidateValid
                try
                    [candidateRate,candidateFlow]=native(candidate);
                catch exception
                    if any(strcmp(exception.identifier,{'ipm:UnresolvedWallIntegralGauge', ...
                            'ipm:DegenerateWallIntegralGauge'}))
                        candidateValid=false;
                    else
                        rethrow(exception);
                    end
                end
            end
            if candidateValid
                candidateResidual=candidate-known-gamma*dt*candidateRate;
                [candidateMerit,candidateDensityNorm,candidateRxNorm]=measure(candidateResidual);
                if candidateMerit<=1 || candidateMerit<=c.preconditionerDefectRequiredReduction*merit
                    trialRho=candidate;rate=candidateRate;flow=candidateFlow;
                    residual=candidateResidual;merit=candidateMerit;
                    densityNorm=candidateDensityNorm;rxNorm=candidateRxNorm;
                    lastLength=1;updates=updates+1;cheapAccepted=cheapAccepted+1;accepted=true;
                    stageOps.rescaling.gaugeRateGuess=[flow.c_l,flow.c_r];
                end
            end
        end
        if ~accepted
        operator=@(v)jacobian(right(v),baseRho,variableScale,rowBalance);
        [linearStep,lastFlag,lastRelres,lastIterations]=gmres(operator, ...
            -rowBalance.*residual(:),c.gmresRestart,c.linearTolerance,c.gmresCycles);
        scaledStep=right(linearStep);step=reshape(variableScale(:).*scaledStep,size(rho));
        if any(~isfinite(step),'all') || ~any(step(:)),break;end
        for backtrack=0:c.maximumLineSearchBacktracks
            length=2^-backtrack;candidate=trialRho+length*step;
            if any(~isfinite(candidate),'all'),continue;end
            try
                [candidateRate,candidateFlow]=native(candidate);
            catch exception
                if any(strcmp(exception.identifier,{'ipm:UnresolvedWallIntegralGauge', ...
                        'ipm:DegenerateWallIntegralGauge'})),continue;end
                rethrow(exception);
            end
            candidateResidual=candidate-known-gamma*dt*candidateRate;
            [candidateMerit,candidateDensityNorm,candidateRxNorm]=measure(candidateResidual);
            if candidateMerit<=1 || candidateMerit<merit*(1-1e-3*length)
                trialRho=candidate;rate=candidateRate;flow=candidateFlow;
                residual=candidateResidual;merit=candidateMerit;
                densityNorm=candidateDensityNorm;rxNorm=candidateRxNorm;
                lastLength=length;updates=updates+1;accepted=true;
                stageOps.rescaling.gaugeRateGuess=[flow.c_l,flow.c_r];
                break;
            end
        end
        end
        if ~accepted,break;end
        if merit<=1,found=true;break;end
    end
    if ~found
        error('ipm:ImplicitStageConvergence', ...
            'SDIRK%d stage %d unresolved after %d updates: density %.3e/%.3e, full DxR %.3e/%.3e.', ...
            c.order,stageIndex,updates,densityNorm,densityTolerance,rxNorm,rxTolerance);
    end
    zi=knownScale;zi(1)=knownScale(1)+gamma*dt*flow.c_l;
    zi(2)=knownScale(2)+gamma*dt*flow.c_omega;
    zi(3)=knownScale(3)+gamma*dt*exp(zi(2)-zi(1));
    denominator=1-gamma*dt*flow.c_l;
    if ~isfinite(denominator) || denominator==0
        error('ipm:ImplicitStageClock','The scale-stage shift equation is singular.');
    end
    zi(4)=(knownScale(4)+gamma*dt*flow.c_r)/denominator;
    zi(5)=z(5)+c.c(stageIndex)*dt;
    if any(~isfinite(zi)),error('ipm:ImplicitStageClock','The scale stage is nonfinite.');end
    rates{stageIndex}=rate;rhoStages{stageIndex}=trialRho;zStages{stageIndex}=zi;
    scaleRates{stageIndex}=[flow.c_l;flow.c_omega;exp(zi(2)-zi(1)); ...
        flow.c_l*zi(4)+flow.c_r;1];
    ops.rescaling.gaugeRateGuess=[flow.c_l,flow.c_r];
    record=struct('stage',stageIndex,'newtonUpdates',updates,'rhsEvaluations',evaluations-before, ...
        'densityResidual',densityNorm,'densityTolerance',densityTolerance, ...
        'fullRxResidual',rxNorm,'rxTolerance',rxTolerance, ...
        'linearFlag',lastFlag,'reportedLinearRelativeResidual',lastRelres, ...
        'gmresIterations',lastIterations,'lastStepLength',lastLength, ...
        'fullTransportCflDiagnostic',dt*flow.generatorCflRate);
    if isfield(c,'preconditionerDefectTrials')
        record.preconditionerDefectTrials=cheapTrials;
        record.preconditionerDefectAccepted=cheapAccepted;
        record.completeJvGmresUpdates=updates-cheapAccepted;
    end
    if isempty(records),records=record;else,records(end+1)=record;end %#ok<AGROW>
end
rhoNew=rhoStages{stageCount};zNew=zStages{stageCount};
info=struct('contract',c.id,'stageCount',stageCount,'stages',records,'rhsEvaluations',evaluations, ...
    'completeNativeRhsUsed',true,'fullRxResidualChecked',true, ...
    'artificialTimeAdvanced',0,'realCanonicalStep',dt, ...
    'explicitCflLimitApplied',false,'nonlinearConvergenceEstablishesOrder',false);

    function [F,flow]=native(field)
        evaluations=evaluations+1;
        stageScale=unpack_scale(z,false);
        stageScale.canonicalTime=z(5)+c.c(stageIndex)*dt;
        [F,flow]=ipm.evolve.flow(field,stageOps,stageScale);
        if any(~isfinite(F),'all')
            error('ipm:ImplicitStageConvergence','The complete implicit native RHS is nonfinite.');
        end
    end
    function [merit,densityNorm,rxNorm]=measure(R)
        densityNorm=max(abs(R)./variableScale,[],'all');
        rxNorm=max(abs(R*ops.Dx'),[],'all');
        merit=max(densityNorm/densityTolerance,rxNorm/rxTolerance);
    end
    function result=jacobian(v,base,scales,rowBalance)
        magnitude=max(abs(v));
        if magnitude==0,result=zeros(size(v));return;end
        h=c.jacobianDifferenceStep/magnitude;
        direction=reshape(scales(:).*v,size(rho));
        plus=native(base+h*direction);minus=native(base-h*direction);
        action=direction-gamma*dt*(plus-minus)/(2*h);
        result=rowBalance.*action(:);
    end
    function [L,U,rowBalance]=preconditioner(flow,scales)
        % M=I-gamma*h*B has positive physical mass. This is not the
        % large artificial PTC block and receives no reaction shift.
        n=numel(rho);ny=size(rho,1);nx=size(rho,2);index=reshape(1:n,ny,nx);
        diagonal=flow.c_omega*ones(ny,nx);rows=cell(1,4);columns=rows;values=rows;
        dx=diff(ops.x);dy=diff(ops.y);vx=flow.transportU1;vy=flow.transportU2;
        weight=max(vx(:,2:end),0)./dx;
        rows{1}=index(:,2:end);columns{1}=index(:,1:end-1);values{1}=weight;
        diagonal(:,2:end)=diagonal(:,2:end)-weight;
        weight=max(-vx(:,1:end-1),0)./dx;
        rows{2}=index(:,1:end-1);columns{2}=index(:,2:end);values{2}=weight;
        diagonal(:,1:end-1)=diagonal(:,1:end-1)-weight;
        weight=max(vy(2:end,:),0)./dy;
        rows{3}=index(2:end,:);columns{3}=index(1:end-1,:);values{3}=weight;
        diagonal(2:end,:)=diagonal(2:end,:)-weight;
        weight=max(-vy(1:end-1,:),0)./dy;
        rows{4}=index(1:end-1,:);columns{4}=index(2:end,:);values{4}=weight;
        diagonal(1:end-1,:)=diagonal(1:end-1,:)-weight;
        row=index(:);column=row;value=diagonal(:);
        for axis=1:4
            row=[row;rows{axis}(:)];column=[column;columns{axis}(:)];value=[value;values{axis}(:)]; %#ok<AGROW>
        end
        B=sparse(row,column,value,n,n);M=speye(n)-gamma*dt*B;
        M=M*spdiags(scales(:),0,n,n);rowBalance=1./max(sum(abs(M),2),realmin);
        balanced=spdiags(rowBalance,0,n,n)*M;
        try
            [L,U]=ilu(balanced,struct('type','ilutp','droptol',c.iluDropTolerance));
        catch exception
            error('ipm:ImplicitStageConvergence','The implicit stage ILU failed: %s',exception.message);
        end
    end
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
if strcmp(ops.rescaling.cOmegaGauge,'wall_omega_tracked_spline_peak') && ...
        (~isfield(cache,'peakTrackingAnchor') || ...
         ~isfield(ops.rescaling,'trackedPeakAnchor') || ...
         ~isequaln(cache.peakTrackingAnchor,ops.rescaling.trackedPeakAnchor))
    return;
end
hit = true;
rhoRate = cache.rhoRate;
scaleRate = cache.scaleRate;
end

function [rhoRate,scaleRate,ops] = stage_rhs(rho,z,ops,isAnisotropic)
scale = unpack_scale(z,isAnisotropic);
if strcmp(ops.rescaling.cOmegaGauge,'wall_omega_tracked_spline_peak')
    ops.rescaling.peakEvaluationTime=scale.canonicalTime;
end
[rhoRate,flow] = ipm.evolve.flow(rho,ops,scale);
if isfield(flow,'peakTrackingAnchor')
    ops.rescaling.trackedPeakAnchor=flow.peakTrackingAnchor;
end
if isAnisotropic
    scaleRate = [flow.c_x;flow.c_y;flow.c_omega; ...
        exp(z(3)-z(1)); ...
        flow.c_x*z(5)+flow.c_r;1];
else
    scaleRate = [flow.c_l;flow.c_omega; ...
        exp(z(2)-z(1)); ...
        flow.c_l*z(4)+flow.c_r;1];
    if ops.rescaling.enabled && ...
            strcmp(ops.rescaling.lengthGauge,'wall_gradient_integral')
        % Continue the same gauge root through the RK stages. This is only
        % an initial guess; each stage rechecks the original exact equations.
        ops.rescaling.gaugeRateGuess = [flow.c_l,flow.c_r];
    end
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
