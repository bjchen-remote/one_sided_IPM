function [rhs,flow] = flow(rho,ops,scale)
%IPM.EVOLVE.FLOW Evaluate the RHS with the geometry's required state.

if strcmp(ops.dynamicScaleGeometry,'anisotropic')
    [rhs,flow] = ipm.evolve.rhs(rho,ops,scale);
else
    [rhs,flow] = ipm.evolve.rhs(rho,ops);
end
end
