function [rhs,flow] = flow(rho,ops,scale)
%IPM.EVOLVE.FLOW Evaluate the RHS with the geometry's required state.

if isfield(ops,'spatialOperators')
    assert(isfield(ops.spatialOperators,'evaluate')&& ...
        isa(ops.spatialOperators.evaluate,'function_handle'), ...
        'ipm:SpatialOperators','The spatial backend must supply complete evaluate(rho,ops,scale).');
    if nargin<3,scale=ipm.evolve.initialScale(ops);end
    [rhs,flow]=ops.spatialOperators.evaluate(rho,ops,scale);
elseif strcmp(ops.dynamicScaleGeometry,'anisotropic')
    [rhs,flow] = ipm.evolve.rhs(rho,ops,scale);
elseif nargin>=3
    [rhs,flow] = ipm.evolve.rhs(rho,ops,scale);
else
    [rhs,flow] = ipm.evolve.rhs(rho,ops);
end
end
