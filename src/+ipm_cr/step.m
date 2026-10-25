function [s,info]=step(s,ctx,tauEnd,limits)
%STEP One accepted explicit step, with complete-stage CFL and real RK events.
opts=ctx.opts;
if nargin>3
    assert(isstruct(limits)&&isscalar(limits),'ipm_cr:StepOptions','Use scalar time controls.');
    names=fieldnames(limits);
    assert(all(ismember(names,{'cfl','maxDt','minDt','timeIntegrator'})), ...
        'ipm_cr:StepOptions','Only time controls can change on an existing context.');
    for k=1:numel(names),opts.(names{k})=limits.(names{k});end
end
assert(isfield(s,'ctxId')&&strcmp(s.ctxId,ctx.id),'ipm_cr:StaleState', ...
    'State belongs to another context. Use refresh only for an identical mathematical setup.');
[s,info]=ipm_cr.time.step(s,ctx.mesh,ctx.wall,opts,tauEnd,ctx);
end
