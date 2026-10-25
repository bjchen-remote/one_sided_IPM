function values=sample(s,ctx,x,y)
%SAMPLE Read W, q, Q, F and velocities at matching coordinate arrays.
% Analytic background terms are evaluated at the queries; residuals are
% interpolated. No advance, patch projection or comparison mask is applied.
assert(isfield(s,'ctxId')&&strcmp(s.ctxId,ctx.id),'ipm_cr:StaleState', ...
    'State belongs to another context.');
values=ipm_cr.diag.fields(s,ctx.mesh,ctx.opts,x,y);
end
