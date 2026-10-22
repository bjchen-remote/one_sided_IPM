function [cL,cR,nodes,firstWeights,secondWeights,phaseContext,peakP,peakK] = ...
        wallIntegralGauge(rho,flow,ops,details,cLGuess)
%IPM.EVOLVE.WALLINTEGRALGAUGE Preserve a local wall-gradient integral.
%   I=integral_{-L}^{L} R_X(X,0) dX uses fixed rescaled coordinates and the
%   exact piecewise-linear clipped quadrature in wallGradientIntegral. Let
%   P=q(a)>0 be the declared peak observable of q=R_X|_{Y=0}. The original
%   selector uses a three-node quadratic; candidate tracked_spline uses one
%   continued C2 natural-spline peak and its exact linear functionals.
%   The coupled conditions on the actual nonlinear wall transport B are
%       I_B/I-P_B/P=0,  a_tau+kappa*a=0.
%   The later amplitude gauge sets c_omega=-P_B/P, so the complete RHS
%   has I_tau=0 and P_tau=0. Curvature is used only to locate the peak and
%   evaluate its phase; no curvature observable is preserved.
%   The nonlinear WENO response is solved with the same strict residual,
%   trust-region and signed-bracket safeguards as the peak phase gauge.

if ~isnumeric(cLGuess) || ~isreal(cLGuess) || ...
        ~isscalar(cLGuess) || ~isfinite(cLGuess)
    error('ipm:DegenerateLengthGauge', ...
        'The wall-integral length-rate initial guess must be finite.');
end
r = ops.rescaling;
% Optional candidate diagnostics never alter the gauge equations or gates.
captureFile=getenv('IPM_GAUGE_CAPTURE_FILE');
evaluationTrace=zeros(0,4);
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
trackedSpline=strcmp(r.cOmegaGauge,'wall_omega_tracked_spline_peak');
if trackedSpline
    peak=details.trackedPeakFunctional;
else
    peak = ipm.evolve.quadraticPeakFunctional( ...
        flow.source(1,:),ops.x,details.wallOmegaPeakIndex);
end
nodes = peak.indices;
a = peak.x;
peakP = peak.value;
peakK = -2*peak.curvature;
if peakK <= 0 || ~isfinite(peakK)
    error('ipm:DegenerateLengthGauge', ...
        'The tracked wall-gradient peak has no resolved negative curvature.');
end
integral = ipm.evolve.wallGradientIntegral(flow.source(1,:),ops,r);
integralNodes = find(integral.support);

if trackedSpline
    firstWeights=peak.firstWeights;
    secondWeights=peak.secondWeights;
else
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
end
if abs(firstWeights*flow.source(1,nodes)') > ...
        1e-8*max(sum(abs(firstWeights).* ...
        abs(flow.source(1,nodes))),realmin)
    error('ipm:DegenerateLengthGauge', ...
        'The declared wall-gradient peak functional is not stationary.');
end

useLocal = strcmp(ops.transportScheme,'weno5_fd') && ...
    strcmp(ops.transportBoundaryMode,'open') && ...
    ~ipm.mesh.isQuadrant(ops);
columns = find(any(ops.Dx([nodes,integralNodes(:)'],:) ~= 0,1));
phaseContext = struct('columns',columns,'vertical',[], ...
    'peakX',a,'peakP',peakP,'peakK',peakK, ...
    'integralNodes',integralNodes(:)', ...
    'integralWeights',integral.weights, ...
    'integralValue',integral.value, ...
    'integralRadius',integral.radius, ...
    'integralSupportPoints',integral.supportPoints, ...
    'integralCondition',integral.condition, ...
    'positionRestoringRate',kappa,'peakFunctional',peak,'iterations',0, ...
    'evaluations',0,'localWallTransport',useLocal);

if any(strcmp(r.scalingContract, ...
        {'wall_integral_generators_v2','wall_integral_centered_generators_v3'}))
    [cL,cR,phaseContext] = ipm.evolve.generatorWallIntegralGauge( ...
        rho,flow,ops,peak,firstWeights,integral,phaseContext,cLGuess);
    return
end

peakU1 = interp1(ops.x,flow.u1(1,:),a,'linear');
cL = cLGuess;
cR = -peakU1-cL*a-kappa*a;
cRInitial = cR;
[residual,vertical,wallRhs] = evaluate(cL,cR,[]);
evaluationCount = 1;
% The later complete-RHS audit checks absolute peak and integral rates.
% Scaling this gate with an outlying trial rate can admit a spurious root.
rateTolerance = @(~,~)[1e-9;1e-9];
trustRadiusL = 4*max(1,abs(cLGuess));
trustRadiusR = 4*max(1,abs(cRInitial));
warmStartUsed = false;
coldRestarted = false;
if isfield(r,'gaugeRateGuess') && isnumeric(r.gaugeRateGuess) && ...
        isreal(r.gaugeRateGuess) && numel(r.gaugeRateGuess) == 2 && ...
        all(isfinite(r.gaugeRateGuess))
    warm = r.gaugeRateGuess(:)';
    if abs(warm(1)-cLGuess) <= trustRadiusL && ...
            abs(warm(2)-cRInitial) <= trustRadiusR
        [warmResidual,warmVertical,warmWallRhs] = evaluate(warm(1),warm(2),[]);
        evaluationCount = evaluationCount+1;
        if max(abs(warmResidual)) < max(abs(residual))
            cL = warm(1);
            cR = warm(2);
            residual = warmResidual;
            vertical = warmVertical;
            wallRhs = warmWallRhs;
            warmStartUsed = true;
        end
    end
end
success = false;
fallbackUsed = false;
% A warm predictor has its own budget; it must not shorten the cold backup.
for iteration = 1:(16+8*double(warmStartUsed))
    tolerance = rateTolerance(cL,cR);
    score = max(abs(residual)./tolerance);
    if score <= 1
        success = true;
        break;
    end
    if warmStartUsed && ~coldRestarted && iteration == 9
        cL = cLGuess;
        cR = cRInitial;
        [residual,vertical,wallRhs] = evaluate(cL,cR,[]);
        evaluationCount = evaluationCount+1;
        coldRestarted = true;
        score = max(abs(residual)./rateTolerance(cL,cR));
        if score <= 1
            success = true;
            break
        end
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
        if warmStartUsed && ~coldRestarted
            cL = cLGuess;
            cR = cRInitial;
            [residual,vertical,wallRhs] = evaluate(cL,cR,[]);
            evaluationCount = evaluationCount+1;
            coldRestarted = true;
            continue
        end
        % A bad Jacobian is not proof that the original equations have no
        % root. Try the signed search before reporting a degenerate gauge.
        [bracketed,trialL,trialR,trialResidual, ...
            trialVertical,trialWallRhs] = bracket_original_gauge(cL,cRInitial);
        if bracketed
            cL = trialL;
            cR = trialR;
            residual = trialResidual;
            vertical = trialVertical;
            wallRhs = trialWallRhs;
            fallbackUsed = true;
            success = true;
            break
        end
        error('ipm:DegenerateWallIntegralGauge', ...
            ['The coupled wall-integral/position response is singular ' ...
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
        if warmStartUsed && ~coldRestarted
            cL = cLGuess;
            cR = cRInitial;
            [residual,vertical,wallRhs] = evaluate(cL,cR,[]);
            evaluationCount = evaluationCount+1;
            coldRestarted = true;
            continue
        end
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
        capture_failure('newton_step_not_reduced');
        error('ipm:UnresolvedWallIntegralGauge', ...
            ['The coupled wall-integral/position Newton step did ' ...
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
    capture_failure('nonlinear_gauge_not_converged');
    error('ipm:UnresolvedWallIntegralGauge', ...
        ['The coupled discrete gauge did not converge: normalized ' ...
        'integral-rate residual %.3e, position-rate residual %.3e.'], ...
        residual(1),residual(2));
end

phaseContext.vertical = vertical;
phaseContext.localWallRhs = wallRhs;
phaseContext.integralRatioRate = residual(1);
phaseContext.positionRateResidual = residual(2);
phaseContext.iterations = iteration;
phaseContext.evaluations = evaluationCount;
phaseContext.bracketFallbackUsed = fallbackUsed;
phaseContext.warmStartUsed = warmStartUsed;
phaseContext.coldRestarted = coldRestarted;

    function [value,verticalCache,selectedWallRhs] = ...
            evaluate(trialL,trialR,verticalCache)
        if useLocal
            [selectedWallRhs,verticalCache] = ...
                ipm.field.weno5WallTransport( ...
                rho,flow.u1,flow.u2,trialL,trialR,ops, ...
                columns,verticalCache);
            qRate = selectedWallRhs*ops.Dx(nodes,columns)';
            integralSamples = selectedWallRhs* ...
                ops.Dx(integralNodes,columns)';
            integralRate = integral.weights(integralNodes)*integralSamples';
        else
            [~,baseRhs] = ipm.evolve.assembleRhs( ...
                rho,flow.u1,flow.u2,trialL,trialL,0,trialR, ...
                ops,ops.transportBoundaryMode);
            selectedWallRhs = baseRhs(1,:);
            qRate = selectedWallRhs*ops.Dx(nodes,:)';
            integralSamples = selectedWallRhs*ops.Dx(integralNodes,:)';
            integralRate = integral.weights(integralNodes)*integralSamples';
            verticalCache = [];
        end
        pRate = peak.weights*qRate';
        phaseRate = firstWeights*qRate';
        aRate = phaseRate/peakK;
        value = [integralRate/integral.value-pRate/peakP; ...
            aRate+kappa*a];
        if ~isempty(captureFile)
            evaluationTrace(end+1,:)=[trialL,trialR,value(:)']; %#ok<AGROW>
        end
    end

    function capture_failure(label)
        if isempty(captureFile),return;end
        try
            debugOps=ops;
            if isfield(debugOps,'poisson'),debugOps=rmfield(debugOps,'poisson');end
            captured=struct('kind','candidate_tracked_spline_gauge_failure', ...
                'peakReconstruction','tracked_natural_cubic_spline', ...
                'failureLabel',label,'rho',rho,'flow',flow,'details',details, ...
                'ops',debugOps,'peak',peak,'integral',integral, ...
                'phaseContext',phaseContext,'residual',residual, ...
                'rates',[cL,cR],'guess',[cLGuess,cRInitial], ...
                'trustRadius',[trustRadiusL,trustRadiusR], ...
                'iteration',iteration,'evaluations',evaluationTrace, ...
                'warmStartUsed',warmStartUsed,'coldRestarted',coldRestarted, ...
                'rateTolerance',[1e-9;1e-9],'useLocal',useLocal);
            if exist('jacobian','var')
                captured.jacobian=jacobian;
                captured.jacobianDeterminant=det(jacobian);
                captured.jacobianRcond=rcond(jacobian);
            end
            save(captureFile,'-struct','captured','-v7.3');
        catch captureError
            warning('ipmtests:GaugeFailureCapture','Capture failed: %s',captureError.message);
        end
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
            integralOnCurve = NaN(size(trialL));
            for sample = 1:numel(trialL)
                try
                    [~,pair] = position_root(trialL(sample),centerR);
                    integralOnCurve(sample) = pair(1);
                catch
                    % A failed inner bracket is not an accepted gauge state.
                end
            end
            for sample = 1:numel(trialL)-1
                if ~all(isfinite(integralOnCurve(sample:sample+1))) || ...
                        integralOnCurve(sample)* ...
                        integralOnCurve(sample+1)>0
                    continue;
                end
                try
                    candidateL = fzero(@(value) ...
                        integral_on_position_curve(value,centerR), ...
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

    function value = integral_on_position_curve(trialL,centerR)
        [~,pair] = position_root(trialL,centerR);
        value = pair(1);
    end
end
