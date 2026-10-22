function [cL,cR,context] = generatorWallIntegralGauge( ...
        rho,flow,ops,peak,firstWeights,integral,context,cLGuess)
%IPM.EVOLVE.GENERATORWALLINTEGRALGAUGE Candidate fixed signed WENO split.
%   V3 subtracts a field-frozen wall-peak reference velocity from physical
%   transport and solves delta=c_r+Uref. V2 retains its uncentered operator.
%   This is a new discrete operator, not a reformulation of the net-velocity
%   WENO operator. Each physical/signed-unit transport is assembled once at
%   this stage's rho. In each sign quadrant its full RHS is affine in
%   |c_l|, |c_r|, and c_omega. Solve the original integral, peak, and phase
%   equations, then audit the actual combined full-array RHS before use.

r=ops.rescaling;
practical=strcmp(r.scalingContract,'wall_integral_centroid_practical_v1');
centered=practical || strcmp(r.scalingContract,'wall_integral_centered_generators_v3');
rootTolerance=1e-9;phasePoint=peak.x;referencePosition=peak.x;
if practical
    rootTolerance=1e-6;phasePoint=context.centroidDisplacement;
    referencePosition=context.centroidValue;
end
if ~any(strcmp(r.scalingContract, ...
        {'wall_integral_generators_v2','wall_integral_centered_generators_v3', ...
        'wall_integral_centroid_practical_v1'})) || ...
        ~strcmp(ops.spatialDiscretization,'high_order') || ...
        ~strcmp(ops.symmetryMode,'half_plane') || ...
        ~strcmp(ops.dynamicScaleGeometry,'isotropic') || ...
        ~strcmp(ops.transportScheme,'weno5_fd') || ...
        ~strcmp(ops.transportBoundaryMode,'open')
    error('ipm:GeneratorGaugeMode','The fixed-generator gauge requires its registered mode.');
end
if practical && isfield(r,'gaugeTransportAssembly') && ...
        ismember(r.gaugeTransportAssembly,{'net_velocity','net_velocity_density_centered', ...
        'net_velocity_density_centered_advective'})
    [cL,cR,context]=net_velocity_gauge( ...
        rho,flow,ops,context,cLGuess,@observables);
    return
end
referenceVelocity=0;physicalU1=flow.u1;
if centered
    referenceVelocity=interp1(ops.x,flow.u1(1,:),referencePosition,'linear');
    physicalU1=flow.u1-referenceVelocity;
end

zero=zeros(size(rho));
[~,physical]=ipm.evolve.assembleRhs( ...
    rho,physicalU1,flow.u2,0,0,0,0,ops,'open');
basis=cell(2,2); % row: dilation/translation; column: negative/positive.
signs=[-1,1];
for index=1:2
    signValue=signs(index);
    [~,basis{1,index}]=ipm.evolve.assembleRhs( ...
        rho,zero,zero,signValue,signValue,0,0,ops,'open');
    [~,basis{2,index}]=ipm.evolve.assembleRhs( ...
        rho,zero,zero,0,0,0,signValue,ops,'open');
end

physicalRates=observables(physical);
amplitudeRates=observables(rho);
target=[0;0;r.peakPositionRestoringRate*phasePoint];
rhs=-physicalRates-target;
dilationRates=zeros(3,2);translationRates=zeros(3,2);
for index=1:2
    dilationRates(:,index)=observables(basis{1,index});
    translationRates(:,index)=observables(basis{2,index});
end
guess=[cLGuess,-interp1(ops.x,flow.u1(1,:),referencePosition,'linear')- ...
    cLGuess*referencePosition-r.peakPositionRestoringRate*phasePoint];
if isfield(r,'gaugeRateGuess') && isnumeric(r.gaugeRateGuess) && ...
        isreal(r.gaugeRateGuess) && numel(r.gaugeRateGuess)==2 && ...
        all(isfinite(r.gaugeRateGuess))
    guess=reshape(r.gaugeRateGuess,1,2);
end
bestScore=Inf;found=false;bestBase=[];bestRates=[NaN,NaN,NaN];
bestResidual=[NaN;NaN;NaN];bestCondition=NaN;bestSigns=[NaN,NaN];
bestTranslationRate=NaN;
candidates=struct([]);
for leftIndex=1:2
    for rightIndex=1:2
        matrix=[dilationRates(:,leftIndex), ...
            translationRates(:,rightIndex),amplitudeRates];
        rowMaximum=max(abs(matrix),[],2);
        condition=0;coefficients=[NaN;NaN;NaN];
        rootResidual=[NaN;NaN;NaN];signConsistent=false;
        if all(isfinite(matrix),'all') && (practical || all(rowMaximum>0))
            if practical,rowMaximum(rowMaximum==0)=1;end
            balanced=matrix./rowMaximum;
            condition=rcond(balanced);
            if isfinite(condition) && (practical || condition>=1e-8)
                if practical && condition<1e-10
                    % Flat global maxima can make I/P rows dependent.
                    % Keep the nearest previous rates in the free subspace.
                    reference=[max(0,signs(leftIndex)*guess(1)); ...
                        max(0,signs(rightIndex)*(guess(2)+referenceVelocity));0];
                    coefficients=reference+lsqminnorm(balanced, ...
                        rhs./rowMaximum-balanced*reference,1e-12);
                else
                    coefficients=balanced\(rhs./rowMaximum);
                end
                if practical && all(coefficients(1:2)>=-1e-10*max(1,max(abs(coefficients))))
                    coefficients(1:2)=max(coefficients(1:2),0);
                end
                signConsistent=all(isfinite(coefficients)) && ...
                    coefficients(1)>=0 && coefficients(2)>=0;
                if signConsistent
                    trialBase=physical+coefficients(1)*basis{1,leftIndex}+ ...
                        coefficients(2)*basis{2,rightIndex};
                    trialFull=trialBase+coefficients(3)*rho;
                    rootResidual=observables(trialFull)+target;
                    if all(isfinite(rootResidual)) && max(abs(rootResidual))<=rootTolerance
                        trialRates=[signs(leftIndex)*coefficients(1), ...
                            signs(rightIndex)*coefficients(2)-referenceVelocity,coefficients(3)];
                        score=norm((trialRates(1:2)-guess)./max(1,abs(guess)),2);
                        if score<bestScore
                            found=true;bestScore=score;bestBase=trialBase;
                            bestRates=trialRates;bestResidual=rootResidual;
                            bestTranslationRate=signs(rightIndex)*coefficients(2);
                            bestCondition=condition;
                            bestSigns=[signs(leftIndex),signs(rightIndex)];
                        end
                    end
                end
            end
        end
        row=struct('signs',[signs(leftIndex),signs(rightIndex)], ...
            'condition',condition,'matrix',matrix,'coefficients',coefficients, ...
            'signConsistent',signConsistent,'residual',rootResidual);
        if isempty(candidates),candidates=row;else,candidates(end+1)=row;end %#ok<AGROW>
    end
end
if ~found
    error('ipm:UnresolvedWallIntegralGauge', ...
        ['The fixed-generator discrete gauge has no sign-consistent, ' ...
        'resolved root at the declared %.3g gate. Quadrants: %s'],rootTolerance,jsonencode(candidates));
end
cL=bestRates(1);cR=bestRates(2);
% Mapped WENO differentiates on the uniform computational axes and divides
% by J. Physical centered node gaps are not its characteristic spacing on
% a remeshed axis; use the same J*dXi as the actual transport operator.
if strcmp(ops.spatialDiscretization,'sixth_order')
    spacingX=ops.transportSpacingX;spacingY=ops.transportSpacingY;
elseif any(strcmp(ops.transportScheme,{'weno5_fd','weno7_fd'}))
    spacingX=ops.metricX*ops.computationalSpacingX;
    spacingY=ops.metricY*ops.computationalSpacingY;
else
    spacingX=ops.hx;spacingY=ops.hy;
end
physicalCfl=max(abs(physicalU1)./spacingX+abs(flow.u2)./spacingY,[],'all');
% Retain the solved delta before actual c_r=delta-Uref is rounded.
translationRate=bestTranslationRate;
dilationCfl=max(abs(ops.X)./spacingX+abs(ops.Y)./spacingY,[],'all');
translationCfl=max(1./spacingX);
context.generatorBaseRhs=bestBase;
context.generatorReferenceVelocity=referenceVelocity;
context.generatorTranslationRate=translationRate;
context.generatorAmplitudeRate=bestRates(3);
context.generatorCflRate=physicalCfl+abs(cL)*dilationCfl+abs(translationRate)*translationCfl;
context.generatorCflComponents=[physicalCfl,abs(cL)*dilationCfl,abs(translationRate)*translationCfl];
context.generatorRootResidual=bestResidual;
context.generatorRootTolerance=rootTolerance;
context.generatorRootCondition=bestCondition;
context.generatorSigns=bestSigns;
context.generatorCandidates=candidates;
context.integralRatioRate=bestResidual(1)-bestResidual(2);
context.positionRateResidual=bestResidual(3);
context.iterations=1;context.evaluations=5;
context.bracketFallbackUsed=false;
context.warmStartUsed=isfield(r,'gaugeRateGuess');
context.coldRestarted=false;
context.localWallTransport=false;

    function values=observables(field)
        % Evaluate the same original wall derivative coefficients accurately;
        % neither the signed basis nor the complete-RHS gate gets a new
        % tolerance. Only the arithmetic of these linear functionals differs.
        [peakHigh,peakLow]=ipm.field.compensatedMultiply( ...
            ops.Dx(peak.indices,:),field(1,:)');
        [pHigh,pLow]=ipm.field.compensatedMultiply( ...
            sparse(peak.weights),peakHigh,peakLow);
        if practical
            [phaseHigh,phaseLow]=ipm.field.compensatedMultiply( ...
                ops.Dx(context.centroidNodes,:),field(1,:)');
            [aHigh,aLow]=ipm.field.compensatedMultiply( ...
                sparse(firstWeights),phaseHigh,phaseLow);
            phaseDenominator=1;
        else
            [aHigh,aLow]=ipm.field.compensatedMultiply( ...
                sparse(firstWeights),peakHigh,peakLow);
            phaseDenominator=-2*peak.curvature;
        end
        selected=find(integral.support);
        [integralHigh,integralLow]=ipm.field.compensatedMultiply( ...
            ops.Dx(selected,:),field(1,:)');
        [iHigh,iLow]=ipm.field.compensatedMultiply( ...
            sparse(integral.weights(selected)),integralHigh,integralLow);
        values=[(iHigh+iLow)/integral.value; ...
            (pHigh+pLow)/peak.value; ...
            (aHigh+aLow)/phaseDenominator];
    end
end

function [cL,cR,context]=net_velocity_gauge(rho,flow,ops,context,cLGuess,observe)
% The I/P/C observables are unchanged. This explicitly selected assembly
% solves the original nonlinear net-velocity WENO wall equations, then
% assembles and audits the complete two-dimensional RHS once at that root.
r=ops.rescaling;rootTolerance=1e-6;
Uref=interp1(ops.x,flow.u1(1,:),context.centroidValue,'linear');
physicalU1=flow.u1-Uref;phase=r.peakPositionRestoringRate*context.centroidDisplacement;
target=[0;0;phase];evaluations=0;maximumEvaluations=160;found=false;polished=false;
verticalCache=[];verticalCL=NaN;
% A geometry-based cold seed avoids inheriting a pole from a different
% discrete assembly. The normal stage-to-stage warm seed is tried first.
[~,physicalWall]=evaluate([0;0]);amplitude=observe(rho);
continuous=[observe(-ops.x.*flow.source(1,:)),observe(-flow.source(1,:)),amplitude];
cold=lsqminnorm(continuous,-observe(physicalWall)-target);
seeds=[cold(1:2)';cLGuess,-cLGuess*context.centroidValue-phase;1,0;0,0;-1,0];
warm=false;
if isfield(r,'gaugeRateGuess') && numel(r.gaugeRateGuess)==2 && all(isfinite(r.gaugeRateGuess))
    guess=reshape(r.gaugeRateGuess,1,2);warm=true;
    if max(abs([guess(1),guess(2)+Uref]))<=8*max(1,max(abs(cold(1:2))))
        seeds=[[guess(1),guess(2)+Uref];seeds];
    else,seeds=[seeds;guess(1),guess(2)+Uref];end
end
last=[NaN;NaN];lastX=[NaN;NaN];condition=NaN;iterations=0;
for seedIndex=1:size(seeds,1)
    x=seeds(seedIndex,:)';
    for iteration=1:14
        iterations=iterations+1;
        if evaluations>=maximumEvaluations,break;end
        [g,wallBase]=evaluate(x);merit=max(abs(g));last=g;lastX=x;
        if merit<=rootTolerance
            % One optional Newton polish avoids algebraic truncation error
            % obscuring SSPRK3 order. Its failure never tightens the gate.
            if merit>1e-10 && ~polished
                J=jacobian(x);trial=x-lsqminnorm(J,g);
                [trialG,trialWall]=evaluate(trial);polished=true;
                if max(abs(trialG))<merit,x=trial;g=trialG;wallBase=trialWall;end
            end
            found=true;break;
        end
        J=jacobian(x);condition=rcond(J);step=-lsqminnorm(J,g);
        stepLimit=4*max(1,norm(x,inf));
        if norm(step,inf)>stepLimit,step=step*(stepLimit/norm(step,inf));end
        accepted=false;
        for damping=2.^(-(0:7))
            trial=x+damping*step;[trialG,~]=evaluate(trial);
            if max(abs(trialG))<merit*(1-1e-3*damping) || max(abs(trialG))<=rootTolerance
                x=trial;accepted=true;break;
            end
        end
        if ~accepted
            for lambda=max(1,norm(J,'fro')^2)*10.^(-4:1)
                trial=x-(J'*J+lambda*eye(2))\(J'*g);[trialG,~]=evaluate(trial);
                if max(abs(trialG))<merit,x=trial;accepted=true;break;end
            end
        end
        if ~accepted,break;end
    end
    if found || evaluations>=maximumEvaluations,break;end
end
if ~found
    error('ipm:UnresolvedWallIntegralGauge', ...
        'Net-velocity I/P/C root unresolved after %d wall evaluations: cl=%g delta=%g residual=[%g %g].', ...
        evaluations,lastX(1),lastX(2),last(1),last(2));
end
if ~isfinite(condition),condition=rcond(jacobian(x));end
baseValues=observe(wallBase);cOmega=-baseValues(2);cL=x(1);delta=x(2);cR=delta-Uref;
[~,base]=ipm.evolve.assembleRhs(rho,physicalU1,flow.u2,cL,cL,0,delta,ops,'open');
wallDifference=max(abs(base(1,:)-wallBase),[],'all')/max(1,max(abs(base(1,:)),[],'all'));
full=base+cOmega*rho;residual=observe(full)+target;
if wallDifference>1e-10 || max(abs(residual))>1e-5 || any(~isfinite(residual))
    error('ipm:UnresolvedWallIntegralGauge', ...
        'Net-velocity complete RHS audit failed: wall difference=%g, I/P/C=[%g %g %g].', ...
        wallDifference,residual(1),residual(2),residual(3));
end
spacingX=ops.metricX*ops.computationalSpacingX;spacingY=ops.metricY*ops.computationalSpacingY;
vx=physicalU1+cL*ops.X+delta;vy=flow.u2+cL*ops.Y;
rate=max(abs(vx)./spacingX+abs(vy)./spacingY,[],'all');
context.generatorBaseRhs=base;context.generatorReferenceVelocity=Uref;
context.generatorTranslationRate=delta;context.generatorAmplitudeRate=cOmega;
% Legacy common column name: in net mode this is the actual net CFL rate.
context.generatorCflRate=rate;context.generatorCflComponents=[rate,0,0];
context.generatorRootResidual=residual;context.generatorRootTolerance=rootTolerance;
context.generatorRootCondition=condition;context.generatorSigns=sign([cL,delta]);
context.generatorCandidates=struct([]);context.integralRatioRate=residual(1)-residual(2);
context.positionRateResidual=residual(3);context.iterations=iterations;
context.evaluations=evaluations+1;context.bracketFallbackUsed=false;
context.warmStartUsed=warm;context.coldRestarted=seedIndex>1;context.localWallTransport=true;
context.gaugeTransportAssembly=r.gaugeTransportAssembly;context.netWallFullRelativeDifference=wallDifference;
context.netRootPolished=polished;

    function [g,b]=evaluate(z)
        evaluations=evaluations+1;
        if ~isequal(verticalCL,z(1)),verticalCache=[];end
        [b,verticalCache]=ipm.field.weno5WallTransport( ...
            rho,physicalU1,flow.u2,z(1),z(2),ops,1:ops.nx,verticalCache);
        verticalCL=z(1);v=observe(b);g=[v(1)-v(2);v(3)+phase];
        assert(all(isfinite(g)),'ipm:NetGaugeFinite','The net wall gauge must remain finite.');
    end
    function J=jacobian(z)
        J=zeros(2,2);
        for j=1:2
            h=1e-4*max(1,abs(z(j)));e=zeros(2,1);e(j)=h;
            [plus,~]=evaluate(z+e);[minus,~]=evaluate(z-e);J(:,j)=(plus-minus)/(2*h);
        end
    end
end
