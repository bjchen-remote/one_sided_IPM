function [A,weightedA] = poissonOperator(ops,kappa)
%IPM.FIELD.POISSONOPERATOR Aspect-aware determinant-one elliptic operator.

if nargin < 2
    kappa = 1;
end
validateattributes(kappa,{'numeric'},{'scalar','real','finite','positive'});

if kappa == 1
    A = ops.A;
    if isempty(A)
        % Explicit diagnostics may request the matrix even when the
        % production solve deliberately keeps only separable factors.
        nxi = ops.nx-2;
        nyi = ops.ny-2;
        A = kron(ops.Tx,speye(nyi))+kron(speye(nxi),ops.Ty);
    end
else
    nxi = ops.nx-2;
    nyi = ops.ny-2;
    Ax = kron(ops.Tx,speye(nyi));
    Ay = kron(speye(nxi),ops.Ty);
    A = (1/kappa)*Ax+kappa*Ay;
end
if strcmp(ops.spatialDiscretization,'high_order') || ...
        strcmp(ops.spatialDiscretization,'sixth_order')
    % The polynomial boundary closures are deliberately full-order rather
    % than diagonal-norm SBP closures, so the interior matrix is generally
    % nonsymmetric.  Direct high-order solves therefore use this operator
    % itself rather than a fictitious symmetrization.
    weightedA = A;
    return;
end
numberOfUnknowns = numel(ops.poissonWeights);
weight = spdiags(ops.poissonWeights,0,numberOfUnknowns,numberOfUnknowns);
weightedA = weight*A;
weightedA = (weightedA+weightedA')/2;
end
