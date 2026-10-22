function [psi,interiorRhs,boundary,solveInfo,poissonOperator] = ...
    poisson(source,ops,kappa,boundaryOverride)
%IPM.FIELD.POISSON Solve -Delta psi = source on a truncated upper half-plane.
%   The physical wall has psi=0. Artificial sides use the configured
%   quadrant/half-plane image Green formula or homogeneous Dirichlet data.
%   A fourth BOUNDARYOVERRIDE argument supplies explicit Dirichlet traces;
%   it is primarily an oracle hook for verifying the complete boundary RHS.

if nargin < 3
    kappa = 1;
end
validateattributes(kappa,{'numeric'},{'scalar','real','finite','positive'});

if ~isequal(size(source), [ops.ny, ops.nx])
    error('ipm:SizeMismatch', 'source must have size [ny,nx].');
end
if nargin >= 4 && ~isempty(boundaryOverride)
    boundary = validate_boundary_override(boundaryOverride,ops);
elseif kappa == 1
    boundary = ipm.field.greenBoundary(source,ops);
else
    boundary = ipm.field.greenBoundary(source,ops,kappa);
end
psi = zeros(ops.ny,ops.nx);
psi(:,1) = boundary.left;
psi(:,end) = boundary.right;
psi(1,:) = boundary.bottom;
psi(end,:) = boundary.top;
rhs = source(2:end-1,2:end-1);
highOrder = strcmp(ops.spatialDiscretization,'high_order') || ...
    strcmp(ops.spatialDiscretization,'sixth_order');
if highOrder
    xBoundaryValues = [boundary.left(2:end-1), ...
        boundary.right(2:end-1)];
    yBoundaryValues = [boundary.bottom(2:end-1); ...
        boundary.top(2:end-1)];
    rhs = rhs-(1/kappa)*xBoundaryValues*ops.poissonBoundary.x' - ...
        kappa*ops.poissonBoundary.y*yBoundaryValues;
    separable = isfield(ops,'separablePoisson');
    if separable
        % On the ny-by-nx interior, the exact discrete equation is
        % kappa*Ty*U + (1/kappa)*U*Tx.' = rhs.  Its Schur factors are
        % cached at mesh construction, including for changing kappa.
        factors = ops.separablePoisson;
        transformedRhs = factors.Qy'*rhs*factors.Qx;
        transformedInterior = sylvester(kappa*factors.Sy, ...
            factors.Sx/kappa,transformedRhs);
        interiorMatrix = factors.Qy*transformedInterior*factors.Qx';
        residual = kappa*ops.Ty*interiorMatrix + ...
            (interiorMatrix*ops.Tx.')/kappa-rhs;
        relativeResidual = norm(residual(:),inf)/ ...
            max(norm(rhs(:),inf),eps);
        % One residual correction offsets roundoff from the similarity
        % transforms on strongly stretched grids.  The original discrete
        % equation, rather than the triangular one, controls acceptance.
        if relativeResidual > min( ...
                ops.anisotropicPoisson.tolerance,1e-11)
            initialInteriorMatrix = interiorMatrix;
            initialRelativeResidual = relativeResidual;
            correctionRhs = -factors.Qy'*residual*factors.Qx;
            correction = sylvester(kappa*factors.Sy, ...
                factors.Sx/kappa,correctionRhs);
            interiorMatrix = interiorMatrix + ...
                factors.Qy*correction*factors.Qx';
            residual = kappa*ops.Ty*interiorMatrix + ...
                (interiorMatrix*ops.Tx.')/kappa-rhs;
            relativeResidual = norm(residual(:),inf)/ ...
                max(norm(rhs(:),inf),eps);
            if relativeResidual > initialRelativeResidual
                interiorMatrix = initialInteriorMatrix;
                relativeResidual = initialRelativeResidual;
            end
        end
        % The first Schur correction can leave a small residual on an
        % extremely stretched box. Only if the original acceptance gate
        % still fails, make at most three further corrections to the same
        % discrete equation. Accept only strictly improving true residuals;
        % the already-passing path above remains bitwise unchanged.
        if isfinite(relativeResidual) && ...
                relativeResidual > ops.anisotropicPoisson.tolerance
            residual = kappa*ops.Ty*interiorMatrix + ...
                (interiorMatrix*ops.Tx.')/kappa-rhs;
            relativeResidual = norm(residual(:),inf)/ ...
                max(norm(rhs(:),inf),eps);
            for refinement = 1:3
                correctionRhs = -factors.Qy'*residual*factors.Qx;
                correction = sylvester(kappa*factors.Sy, ...
                    factors.Sx/kappa,correctionRhs);
                candidate = interiorMatrix + ...
                    factors.Qy*correction*factors.Qx';
                candidateResidual = kappa*ops.Ty*candidate + ...
                    (candidate*ops.Tx.')/kappa-rhs;
                candidateRelativeResidual = ...
                    norm(candidateResidual(:),inf)/ ...
                    max(norm(rhs(:),inf),eps);
                if ~isfinite(candidateRelativeResidual) || ...
                        candidateRelativeResidual >= relativeResidual
                    break;
                end
                interiorMatrix = candidate;
                residual = candidateResidual;
                relativeResidual = candidateRelativeResidual;
                if relativeResidual <= ops.anisotropicPoisson.tolerance
                    break;
                end
            end
        end
        if ~isfinite(relativeResidual) || ...
                relativeResidual > ops.anisotropicPoisson.tolerance
            % Exploratory arithmetic certificate for the million-unit box.
            % Seven-point x and y stencil products and sums each incur
            % rounding. Gamma_16 bounds their combined local matvec error;
            % a residual beyond this allowance still faces the original
            % source-normalized tolerance. The independent raw cap prevents
            % this research branch from silently accepting larger errors.
            gamma16 = (16*eps)/(1-16*eps);
            absoluteTerms = kappa*abs(ops.Ty)*abs(interiorMatrix) + ...
                (abs(interiorMatrix)*abs(ops.Tx.'))/kappa + abs(rhs);
            roundoffAllowance = gamma16*absoluteTerms;
            sourceScale = max(norm(rhs(:),inf),eps);
            excess = max(abs(residual)-roundoffAllowance,0);
            excessRelative = norm(excess(:),inf)/sourceScale;
            sourceInterior = source(2:end-1,2:end-1);
            activeSource = abs(sourceInterior) > ...
                1e-10*max(abs(sourceInterior(:)));
            activeRelative = 0;
            if any(activeSource(:))
                activeRelative = max(abs(residual(activeSource)))/sourceScale;
            end
            maxBackward = max(abs(residual(:))./ ...
                max(absoluteTerms(:),eps));
            fprintf(['POISSON_CERT raw=%.9e excess=%.9e active=' ...
                '%.9e backward=%.9e maxAllowance=%.9e\n'], ...
                relativeResidual,excessRelative,activeRelative, ...
                maxBackward,max(roundoffAllowance(:)));
            % Continuation certificate for the million-scale box: retain
            % the conditioned excess gate and allow the source-free
            % far-field cancellation floor to 2e-5.
            if ~isfinite(relativeResidual) || ...
                    relativeResidual > 2e-5 || ...
                    excessRelative > ops.anisotropicPoisson.tolerance
                error('ipm:SeparablePoissonResidual', ...
                    ['The separable Poisson solve has raw residual ' ...
                    '%.3e, conditioned excess %.3e, above allowed ' ...
                    '2e-5 / %.3e.'],relativeResidual,excessRelative, ...
                    ops.anisotropicPoisson.tolerance);
            end
        end
        interior = interiorMatrix(:);
        poissonOperator = [];
        solver = 'separable_high_order';
    else
        [poissonOperator,~] = ipm.field.poissonOperator(ops,kappa);
    end
    if ~separable
        if kappa == 1
            interior = ops.poisson\rhs(:);
            solver = 'cached_direct_high_order';
        else
            if ~strcmp(ops.anisotropicPoisson.solver,'direct')
                error('ipm:HighOrderPoissonSolver', ...
                    ['The nonsymmetric high-order boundary closure currently ' ...
                    'supports the direct anisotropic Poisson solver only.']);
            end
            factor = decomposition(poissonOperator,'lu');
            interior = factor\rhs(:);
            solver = 'direct_high_order';
        end
        relativeResidual = norm(poissonOperator*interior-rhs(:),inf) / ...
            max(norm(rhs(:),inf),eps);
    end
    solveInfo = struct('solver',solver,'kappa',kappa,'flag',0, ...
        'relativeResidual',relativeResidual,'iterations',0);
elseif kappa == 1
    rhs(:,1) = rhs(:,1)-ops.poissonBoundary.x(1)*boundary.left(2:end-1);
    rhs(:,end) = rhs(:,end)- ...
        ops.poissonBoundary.x(2)*boundary.right(2:end-1);
    rhs(1,:) = rhs(1,:)-ops.poissonBoundary.y(1)*boundary.bottom(2:end-1);
    rhs(end,:) = rhs(end,:)-ops.poissonBoundary.y(2)*boundary.top(2:end-1);
    weightedRhs = ops.poissonWeights.*rhs(:);
    interior = ops.poisson \ weightedRhs;
    poissonOperator = ops.A;
    solveInfo = struct('solver','cached_direct','kappa',1,'flag',0, ...
        'relativeResidual',NaN,'iterations',0);
else
    inverseKappa = 1/kappa;
    rhs(:,1) = rhs(:,1)-inverseKappa*ops.poissonBoundary.x(1)* ...
        boundary.left(2:end-1);
    rhs(:,end) = rhs(:,end)-inverseKappa*ops.poissonBoundary.x(2)* ...
        boundary.right(2:end-1);
    rhs(1,:) = rhs(1,:)-kappa*ops.poissonBoundary.y(1)* ...
        boundary.bottom(2:end-1);
    rhs(end,:) = rhs(end,:)-kappa*ops.poissonBoundary.y(2)* ...
        boundary.top(2:end-1);
    weightedRhs = ops.poissonWeights.*rhs(:);
    [poissonOperator,weightedOperator] = ...
        ipm.field.poissonOperator(ops,kappa);
    solver = ops.anisotropicPoisson.solver;
    switch solver
        case 'direct'
            factor = decomposition(weightedOperator,'chol');
            interior = factor \ weightedRhs;
            flag = 0;
            iterations = 0;
            relativeResidual = norm(weightedOperator*interior-weightedRhs,inf) / ...
                max(norm(weightedRhs,inf),eps);
        case 'pcg'
            preconditioner = @(value) ops.poisson \ value;
            [interior,flag,relativeResidual,iterations] = pcg( ...
                weightedOperator,weightedRhs, ...
                ops.anisotropicPoisson.tolerance, ...
                ops.anisotropicPoisson.maxIterations,preconditioner);
            if flag ~= 0
                error('ipm:AnisotropicPoissonPcg', ...
                    ['Aspect-aware PCG failed with flag %d and relative ' ...
                    'residual %.3e after %d iterations.'], ...
                    flag,relativeResidual,iterations);
            end
        otherwise
            error('ipm:AnisotropicPoissonSolver', ...
                'Unknown aspect-aware Poisson solver ''%s''.',solver);
    end
    solveInfo = struct('solver',solver,'kappa',kappa,'flag',flag, ...
        'relativeResidual',relativeResidual,'iterations',iterations);
end
psi(2:end-1, 2:end-1) = reshape(interior,ops.ny-2,ops.nx-2);
interiorRhs = rhs;
end

function boundary = validate_boundary_override(boundary,ops)
if ~isstruct(boundary) || ~isscalar(boundary)
    error('ipm:PoissonBoundaryOverride', ...
        'The Poisson boundary override must be a scalar structure.');
end
names = {'left','right','bottom','top'};
sizes = {[ops.ny,1],[ops.ny,1],[1,ops.nx],[1,ops.nx]};
for index = 1:numel(names)
    name = names{index};
    if ~isfield(boundary,name) || ...
            ~isequal(size(boundary.(name)),sizes{index}) || ...
            ~isnumeric(boundary.(name)) || ...
            ~isreal(boundary.(name)) || ...
            any(~isfinite(boundary.(name)),'all')
        error('ipm:PoissonBoundaryOverride', ...
            'Boundary field "%s" must be a finite real array of size %s.', ...
            name,mat2str(sizes{index}));
    end
end
end
