function solution=withInnerTrace(p,base,inner)
%WITHINNERTRACE Change only the inner RHS using cached homogeneous modes.
% BASE already solves the SAME modal systems with the desired source and
% outer RHS. Updating its first RHS row changes neither the source nor the
% outer condition. This exact linear-system decomposition avoids repeating
% source transforms and radial factorizations for a second direct solve.
% The caller must retain the original complete-system accuracy audit.
if nargin<3||isempty(inner),inner=zeros(1,p.nt);end
if isscalar(inner),inner=repmat(inner,1,p.nt);end
assert(isnumeric(inner)&&isreal(inner)&&numel(inner)==p.nt&&all(isfinite(inner(:))), ...
    'ipm_cr:polar:TraceSize','inner must contain one finite real value per angular node.');
inner=reshape(inner,1,p.nt);
if any(abs(inner([1,end]))>1e-11*max(1,max(abs(inner))))
    error('ipm_cr:polar:WallTrace', ...
        'inner trace must vanish at the physical walls theta=0,pi.');
end
inner([1,end])=0;
assert(isstruct(base)&&all(isfield(base,{'psi','modal','rhs'}))&& ...
    isequal(size(base.psi),[p.nr,p.nt])&& ...
    isequal(size(base.modal),[p.nr,numel(p.mu)])&& ...
    isequal(size(base.rhs),size(base.modal)), ...
    'ipm_cr:polar:BaseSolution','Use a complete base solution on this modal geometry.');
innerModes=inner(p.interior)*p.Vinv.';
change=innerModes-base.rhs(1,:);
solution=base;
solution.rhs(1,:)=innerModes;
if ~any(change)&&isequal(base.modal(1,:),innerModes),return;end
if isfield(solution,'refinement')
    solution.refinement.scope='base_source_outer_solve_before_this_inner_trace_update';
end
solution.modal=base.modal+p.harmonicModes.*change;
% The first original modal row is Dirichlet, even if the requested RHS did
% not change but BASE carries roundoff in that row.
solution.modal(1,:)=innerModes;
solution.psi=zeros(p.nr,p.nt);
solution.psi(:,p.interior)=solution.modal*p.V.';
end
