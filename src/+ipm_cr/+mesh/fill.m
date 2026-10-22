function [rho,info]=fill(input,g)
%SYNCHRONIZE Exchange only the overset ghost collars from old ACTIVE data.
% Call before EVERY explicit stage RHS, including its Poisson and gauges.
% Active overlap samples are independent evolved degrees of freedom and
% are never projected onto the other patch by this operation.
assert(isstruct(input)&&isscalar(input)&&all(isfield(input,{'core','polar'}))&& ...
    isnumeric(input.core)&&isreal(input.core)&&isequal(size(input.core),size(g.core.X))&& ...
    isnumeric(input.polar)&&isreal(input.polar)&&isequal(size(input.polar),size(g.polar.X))&& ...
    all(isfinite(input.core),'all')&&all(isfinite(input.polar),'all'), ...
    'ipm_cr:ExplicitSynchronization','Use finite native Cartesian and annular density arrays.');
e=g.exchange;
if isfield(e,'coreDifferencePlan')
    newCoreGhost=apply(e.coreDifferencePlan,input.polar(e.polarDonorIndices));
    newPolarGhost=apply(e.polarDifferencePlan,input.core(e.coreDonorIndices));
else
    % Historical read-only checkpoints retain their original geometry.
    newCoreGhost=e.coreFromPolar*input.polar(e.polarDonorIndices);
    newPolarGhost=e.polarFromCore*input.core(e.coreDonorIndices);
end
assert(all(isfinite(newCoreGhost))&&all(isfinite(newPolarGhost)), ...
    'ipm_cr:ExplicitSynchronization','An active-donor ghost interpolation produced a nonfinite value.');
rho=input;
rho.core(e.coreGhostIndices)=newCoreGhost;
rho.polar(e.polarGhostIndices)=newPolarGhost;
if nargout>1
    info=struct('kind','simultaneous_active_donor_overset_ghost_exchange_v2', ...
        'coreGhostCount',numel(newCoreGhost),'polarGhostCount',numel(newPolarGhost), ...
        'maximumCoreGhostUpdate',max(abs(newCoreGhost-input.core(e.coreGhostIndices))), ...
        'maximumPolarGhostUpdate',max(abs(newPolarGhost-input.polar(e.polarGhostIndices))), ...
        'activeValuesChanged',false,'donorsAreActiveOnly',true, ...
        'coreMapInfinityNorm',e.coreMapInfinityNorm,'polarMapInfinityNorm',e.polarMapInfinityNorm, ...
        'interpolationWidth',e.interpolationWidth);
    info.localDifferenceForm=isfield(e,'coreDifferencePlan');
end
end

function v=apply(p,donor)
anchor=donor(p.anchor);
v=anchor+p.sum*(donor(p.column)-anchor(p.row));
end
