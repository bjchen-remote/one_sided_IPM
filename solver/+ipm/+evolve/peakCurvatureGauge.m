function [cL,cR,nodes,firstWeights,secondWeights,phaseContext,peakP,peakK] = ...
        peakCurvatureGauge(rho,flow,ops,details,cLGuess)
%IPM.EVOLVE.PEAKCURVATUREGAUGE Pin wall-gradient peak curvature and position.
%   On the active three-node quadratic q=R_X|_{Y=0}, let P=q(a)>0 and
%   a be its vertex. The length condition uses a fixed physical-window
%   polynomial curvature C=q_fit''(a)<0, whose moving-vertex rate is
%       C_tau=W_C*q_tau+q_fit'''(a)*a_tau.
%   The two independent discrete conditions are
%       C_tau/C-P_tau/P=0,  a_tau+kappa*a=0.
%   The first is independent of c_omega; the existing quadratic-peak
%   amplitude gauge subsequently sets c_omega=-P_tau/P.  Thus both P and
%   K have zero rate under the complete RHS.  The second is the gradient
%   flow of the position penalty kappa*a^2/2.  All trial rates use the
%   actual nonlinear wall transport, including WENO flux splitting.
%   User-facing choices (length/amplitude/position gauges, kappa and
%   curvature window radius) are collected in main.m. The
%   tolerance and Jacobian probes below are numerical root-solver controls,
%   not parameters in the IPM equation or a prescribed value of C_s.

if ~isnumeric(cLGuess) || ~isreal(cLGuess) || ...
        ~isscalar(cLGuess) || ~isfinite(cLGuess)
    error('ipm:DegenerateLengthGauge', ...
        'The peak-curvature length-rate initial guess must be finite.');
end
r = ops.rescaling;
if ~isfield(r,'peakPositionRestoringRate') || ...
        ~isnumeric(r.peakPositionRestoringRate) || ...
        ~isreal(r.peakPositionRestoringRate) || ...
        ~isscalar(r.peakPositionRestoringRate) || ...
        ~isfinite(r.peakPositionRestoringRate) || ...
        r.peakPositionRestoringRate < 0
    error('ipm:BadPeakPositionRestoringRate', ...
        'The peak-position restoring rate must be nonnegative and finite.');
end
kappa = r.peakPositionRestoringRate;
peak = ipm.evolve.quadraticPeakFunctional( ...
    flow.source(1,:),ops.x,details.wallOmegaPeakIndex);
nodes = peak.indices;
a = peak.x;
peakP = peak.value;
peakK = -2*peak.curvature;
if peakK <= 0 || ~isfinite(peakK)
    error('ipm:DegenerateLengthGauge', ...
        'The tracked wall-gradient peak has no resolved negative curvature.');
end
[curvatureNodes,smoothCurvature,curvatureThird, ...
    curvatureWeights,curvatureCondition] = ...
    ipm.evolve.smoothPeakCurvatureFunctional( ...
    flow.source(1,:),ops.x,peak,r.peakCurvatureRadius);

firstWeights = zeros(1,3);
secondWeights = zeros(1,3);
for j = 1:3
    other = setdiff(1:3,j);
    xj = ops.x(nodes(j));
    xa = ops.x(nodes(other(1)));
    xb = ops.x(nodes(other(2)));
    denominator = (xj-xa)*(xj-xb);
    firstWeights(j) = (2*a-xa-xb)/denominator;
    secondWeights(j) = 2/denominator;
end
if abs(firstWeights*flow.source(1,nodes)') > ...
        1e-8*max(sum(abs(firstWeights).* ...
        abs(flow.source(1,nodes))),realmin)
    error('ipm:DegenerateLengthGauge', ...
        'The tracked quadratic wall-gradient peak is not stationary.');
end

useLocal = strcmp(ops.transportScheme,'weno5_fd') && ...
    strcmp(ops.transportBoundaryMode,'open') && ...
    ~ipm.mesh.isQuadrant(ops);
columns = find(any(ops.Dx([nodes,curvatureNodes(:)'],:) ~= 0,1));
phaseContext = struct('columns',columns,'vertical',[], ...
    'peakX',a,'peakP',peakP,'peakK',peakK, ...
    'smoothCurvature',smoothCurvature, ...
    'curvatureNodes',curvatureNodes(:)', ...
    'curvatureWeights',curvatureWeights, ...
    'curvatureThird',curvatureThird, ...
    'curvatureCondition',curvatureCondition, ...
    'positionRestoringRate',kappa,'iterations',0, ...
    'evaluations',0,'localWallTransport',useLocal);

peakU1 = interp1(ops.x,flow.u1(1,:),a,'linear');
cL = cLGuess;
cR = -peakU1-cL*a-kappa*a;
cRInitial = cR;
[residual,vertical,wallRhs] = evaluate(cL,cR,[]);
evaluationCount = 1;
% The later complete-RHS audit checks absolute peak and curvature rates.
% Scaling this gate with an outlying trial rate can admit a spurious root.
rateTolerance = @(~,~)[1e-9;1e-9];
trustRadiusL = 4*max(1,abs(cLGuess));
trustRadiusR = 4*max(1,abs(cRInitial));
success = false;
fallbackUsed = false;
for iteration = 1:8
    tolerance = rateTolerance(cL,cR);
    score = max(abs(residual)./tolerance);
    if score <= 1
        success = true;
        break;
    end

    % Jacobian probe for the TWO discrete gauge equations in (c_l,c_r).
    % A continuous linear identity is insufficient because local-alpha
    % WENO flux splitting changes with these rates. These finite-difference
    % steps were tested against the signed-bracket reference; that fallback
    % and the strict 1e-9 full-residual acceptance gate remain intact.
    hL = 1e-5*max(1,abs(cL));
    hR = 1e-6*max(1,abs(cR));
    [leftMinus,~,~] = evaluate(cL-hL,cR,[]);
    [leftPlus,~,~] = evaluate(cL+hL,cR,[]);
    [rightMinus,~,~] = evaluate(cL,cR-hR,vertical);
    [rightPlus,~,~] = evaluate(cL,cR+hR,vertical);
    evaluationCount = evaluationCount+4;
    jacobian = [(leftPlus-leftMinus)/(2*hL), ...
        (rightPlus-rightMinus)/(2*hR)];
    if any(~isfinite(jacobian),'all') || rcond(jacobian) < 1e-8
        error('ipm:DegeneratePeakCurvatureGauge', ...
            ['The coupled peak-curvature/position response is singular ' ...
            'or poorly conditioned.']);
    end
    step = -jacobian\residual;
    stepLimit = 4*max([1,abs(cL),abs(cR)]);
    if max(abs(step)) > stepLimit
        step = step*(stepLimit/max(abs(step)));
    end

    accepted = false;
    bestTrialScore = Inf;
    for backtrack = 0:8
        damping = 2^(-backtrack);
        trialL = cL+damping*step(1);
        trialR = cR+damping*step(2);
        if abs(trialL-cLGuess) > trustRadiusL || ...
                abs(trialR-cRInitial) > trustRadiusR
            continue;
        end
        [trialResidual,trialVertical,trialWallRhs] = ...
            evaluate(trialL,trialR,[]);
        evaluationCount = evaluationCount+1;
        trialScore = max(abs(trialResidual)./ ...
            rateTolerance(trialL,trialR));
        bestTrialScore = min(bestTrialScore,trialScore);
        if isfinite(trialScore) && ...
                (trialScore < score*(1-1e-3*damping) || trialScore <= 1)
            cL = trialL;
            cR = trialR;
            residual = trialResidual;
            vertical = trialVertical;
            wallRhs = trialWallRhs;
            accepted = true;
            break;
        end
    end
    if ~accepted
        % The local-alpha WENO response is nonlinear in both rates. A
        % finite-difference Newton direction can stall even while the two
        % original discrete equations have a nearby signed root. Keep the
        % strict residual gate and seek that root on the position-null curve.
        % Anchor the inner position search at the physical peak-transport
        % estimate: the failed Newton iterate can cross to another root
        % branch before this fallback is reached.
        [bracketed,trialL,trialR,trialResidual, ...
            trialVertical,trialWallRhs] = ...
            bracket_original_gauge(cL,cRInitial);
        if bracketed
            cL = trialL;
            cR = trialR;
            residual = trialResidual;
            vertical = trialVertical;
            wallRhs = trialWallRhs;
            fallbackUsed = true;
            success = true;
            break;
        end
        error('ipm:UnresolvedPeakCurvatureGauge', ...
            ['The coupled peak-curvature/position Newton step did ' ...
            'not reduce the discrete gauge residual (iteration %d, ' ...
            'score %.3e, best trial %.3e, G %.3e, F %.3e, ' ...
            'c_l %.12g, c_r %.12g, Jacobian ' ...
            '[%.3e %.3e; %.3e %.3e], step [%.3e %.3e]).'], ...
            iteration,score,bestTrialScore,residual(1),residual(2), ...
            cL,cR,jacobian(1,1),jacobian(1,2), ...
            jacobian(2,1),jacobian(2,2),step(1),step(2));
    end
end
if ~success
    tolerance = rateTolerance(cL,cR);
    success = max(abs(residual)./tolerance) <= 1;
end
if ~success
    [bracketed,trialL,trialR,trialResidual, ...
        trialVertical,trialWallRhs] = ...
        bracket_original_gauge(cL,cRInitial);
    if bracketed
        cL = trialL;
        cR = trialR;
        residual = trialResidual;
        vertical = trialVertical;
        wallRhs = trialWallRhs;
        fallbackUsed = true;
        success = true;
    end
end
if ~success
    error('ipm:UnresolvedPeakCurvatureGauge', ...
        ['The coupled discrete gauge did not converge: normalized ' ...
        'curvature-rate residual %.3e, position-rate residual %.3e.'], ...
        residual(1),residual(2));
end

phaseContext.vertical = vertical;
phaseContext.localWallRhs = wallRhs;
phaseContext.curvatureRatioRate = residual(1);
phaseContext.positionRateResidual = residual(2);
phaseContext.iterations = iteration;
phaseContext.evaluations = evaluationCount;
phaseContext.bracketFallbackUsed = fallbackUsed;

    function [value,verticalCache,selectedWallRhs] = ...
            evaluate(trialL,trialR,verticalCache)
        if useLocal
            [selectedWallRhs,verticalCache] = ...
                ipm.field.weno5WallTransport( ...
                rho,flow.u1,flow.u2,trialL,trialR,ops, ...
                columns,verticalCache);
            qRate = selectedWallRhs*ops.Dx(nodes,columns)';
            curvatureQRate = selectedWallRhs* ...
                ops.Dx(curvatureNodes,columns)';
        else
            [~,baseRhs] = ipm.evolve.assembleRhs( ...
                rho,flow.u1,flow.u2,trialL,trialL,0,trialR, ...
                ops,ops.transportBoundaryMode);
            selectedWallRhs = baseRhs(1,:);
            qRate = selectedWallRhs*ops.Dx(nodes,:)';
            curvatureQRate = selectedWallRhs* ...
                ops.Dx(curvatureNodes,:)';
            verticalCache = [];
        end
        pRate = peak.weights*qRate';
        phaseRate = firstWeights*qRate';
        aRate = phaseRate/peakK;
        curvatureRate = curvatureWeights*curvatureQRate' + ...
            curvatureThird*aRate;
        value = [curvatureRate/smoothCurvature-pRate/peakP; ...
            aRate+kappa*a];
    end

    function [found,rootL,rootR,rootResidual,rootVertical,rootWallRhs] = ...
            bracket_original_gauge(centerL,centerR)
        found = false;
        rootL = NaN;
        rootR = NaN;
        rootResidual = [NaN;NaN];
        rootVertical = [];
        rootWallRhs = [];
        % Search the strain-based initial guess before the failed Newton
        % iterate. A Newton step can leave the original root branch by an
        % order-one distance, so a narrow search about that iterate alone
        % misses an otherwise well-resolved signed root.
        scaleL = max(1,abs(cLGuess));
        trialGroups = { ...
            cLGuess+[-1,-0.5,0,0.5,1]*scaleL, ...
            cLGuess+[-4,-2,-1,-0.5,-0.25,-0.1,-0.05, ...
                -0.02,0,0.02,0.05,0.1,0.25,0.5,1,2,4]*scaleL, ...
            centerL+[-0.02,-0.01,-0.005,0,0.005,0.01,0.02]};
        options = optimset('TolX',1e-12,'Display','off');
        for group = 1:numel(trialGroups)
            trialL = sort(unique(trialGroups{group}));
            curvatureOnCurve = NaN(size(trialL));
            for sample = 1:numel(trialL)
                try
                    [~,pair] = position_root(trialL(sample),centerR);
                    curvatureOnCurve(sample) = pair(1);
                catch
                    % A failed inner bracket is not an accepted gauge state.
                end
            end
            for sample = 1:numel(trialL)-1
                if ~all(isfinite(curvatureOnCurve(sample:sample+1))) || ...
                        curvatureOnCurve(sample)* ...
                        curvatureOnCurve(sample+1)>0
                    continue;
                end
                try
                    candidateL = fzero(@(value) ...
                        curvature_on_position_curve(value,centerR), ...
                        trialL(sample:sample+1),options);
                    [candidateR,~] = position_root(candidateL,centerR);
                    [candidateResidual,candidateVertical,candidateWallRhs] = ...
                        evaluate(candidateL,candidateR,[]);
                    evaluationCount = evaluationCount+1;
                    candidateScore = max(abs(candidateResidual)./ ...
                        rateTolerance(candidateL,candidateR));
                    if isfinite(candidateScore) && candidateScore <= 1 && ...
                            abs(candidateL-cLGuess) <= trustRadiusL && ...
                            abs(candidateR-cRInitial) <= trustRadiusR
                        found = true;
                        rootL = candidateL;
                        rootR = candidateR;
                        rootResidual = candidateResidual;
                        rootVertical = candidateVertical;
                        rootWallRhs = candidateWallRhs;
                        return;
                    end
                catch
                    % Try the next signed outer bracket, then retain the error.
                end
            end
        end
    end

    function [rootR,pair] = position_root(trialL,centerR)
        halfWidth = 0.01;
        bracket = [];
        for attempt = 1:14
            candidate = centerR+[-halfWidth,halfWidth];
            leftValue = position_residual(trialL,candidate(1));
            rightValue = position_residual(trialL,candidate(2));
            if isfinite(leftValue) && isfinite(rightValue) && ...
                    leftValue*rightValue <= 0
                bracket = candidate;
                break;
            end
            halfWidth = 2*halfWidth;
        end
        if isempty(bracket)
            error('ipm:UnbracketedPeakPosition', ...
                'No signed position root near the failed Newton state.');
        end
        options = optimset('TolX',1e-12,'Display','off');
        rootR = fzero(@(value) position_residual(trialL,value), ...
            bracket,options);
        [pair,~,~] = evaluate(trialL,rootR,[]);
        evaluationCount = evaluationCount+1;
    end

    function value = position_residual(trialL,trialR)
        pair = evaluate(trialL,trialR,[]);
        evaluationCount = evaluationCount+1;
        value = pair(2);
    end

    function value = curvature_on_position_curve(trialL,centerR)
        [~,pair] = position_root(trialL,centerR);
        value = pair(1);
    end
end
