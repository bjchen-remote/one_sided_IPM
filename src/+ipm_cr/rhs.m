function evaluation=rhs(s,ctx,wantDiagnostics)
%RHS Complete same-state RHS on the context's immutable mathematical setup.
if nargin<3,wantDiagnostics=true;end
assert(isfield(s,'ctxId')&&strcmp(s.ctxId,ctx.id),'ipm_cr:StaleState','State and context differ.');
guess=[s.rhs.c_l,s.rhs.c_r];
if isfield(s,'eta')
    evaluation=ipm_cr.evolve.residual(s.eta,s.scale,ctx.mesh,ctx.wall,ctx.opts,guess,wantDiagnostics,ctx);
else
    evaluation=ipm_cr.evolve.rhs(s.W,ctx.mesh,ctx.wall,ctx.opts,guess,[],wantDiagnostics,ctx);
end
end
