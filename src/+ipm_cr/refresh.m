function s=refresh(s,ctx)
%REFRESH Rebind a recovered state only to the exact same mathematical setup.
assert(isfield(s,'ctxKey')&&isequaln(s.ctxKey,ctx.key),'ipm_cr:StateSetup', ...
    'Mesh, wall reference, background or RHS method changed; the old state cannot be rebound.');
s.ctxId=ctx.id;
[s.rhs,s.W]=ipm_cr.rhs(s,ctx);
end
