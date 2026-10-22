function s=refresh(s,ctx)
%REFRESH Rebind a recovered state only to the exact same mathematical setup.
assert(isfield(s,'ctxKey')&&isequaln(s.ctxKey,ctx.key),'ipm_cr:StateSetup', ...
    'Mesh, wall reference, background or RHS method changed; the old state cannot be rebound.');
s.ctxId=ctx.id;
if isfield(s,'eta')
    [s.rhs,s.W]=ipm_cr.evolve.residual(s.eta,s.scale,ctx.mesh,ctx.wall,ctx.opts, ...
        [s.rhs.c_l,s.rhs.c_r],true,ctx);
else
    [s.rhs,s.W]=ipm_cr.evolve.rhs(s.W,ctx.mesh,ctx.wall,ctx.opts, ...
        [s.rhs.c_l,s.rhs.c_r],[],true,ctx);
end
end
