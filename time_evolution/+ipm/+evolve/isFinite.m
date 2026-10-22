function finite = isFinite(state)
%IPM.EVOLVE.ISFINITE Check the numerical fields required by the next step.
%   Gauge diagnostics may intentionally contain NaN when inactive, so this
%   predicate checks only fields that participate in evolution or physical
%   reconstruction.

requiredFlow = {'u1','u2','psi','source','transportU1','transportU2', ...
    'poissonResidual','timeSpeed','c_l','c_x','c_y','c_omega','c_r', ...
    'canonicalCL','canonicalCX','canonicalCY','canonicalCOmega', ...
    'canonicalCR','conservativeSource','canonicalConservativeSource'};

finite = isfield(state,'rho') && finite_numeric(state.rho) && ...
    isfield(state,'flow') && isstruct(state.flow) && ...
    isfield(state,'scale') && isstruct(state.scale);
if ~finite
    return;
end
for index = 1:numel(requiredFlow)
    name = requiredFlow{index};
    if ~isfield(state.flow,name) || ~finite_numeric(state.flow.(name))
        finite = false;
        return;
    end
end

scaleNames = fieldnames(state.scale);
for index = 1:numel(scaleNames)
    if ~finite_numeric(state.scale.(scaleNames{index}))
        finite = false;
        return;
    end
end

stateScalars = {'normalizedTime','step','mass0','rhoRange0'};
for index = 1:numel(stateScalars)
    name = stateScalars{index};
    if isfield(state,name) && ~finite_numeric(state.(name))
        finite = false;
        return;
    end
end
end

function finite = finite_numeric(value)
finite = isnumeric(value) && isreal(value) && ...
    ~isempty(value) && all(isfinite(value),'all');
end
