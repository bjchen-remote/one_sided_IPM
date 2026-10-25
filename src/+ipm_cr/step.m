function [s,info]=step(s,ctx,tauEnd,controls)
%STEP One accepted explicit step, with complete-stage CFL and real RK events.
% A prepared context is required; controls contain time policy only.
assert(nargin>=3&&isa(ctx,'ipm_cr.context')&& ...
    isfield(s,'ctxId')&&strcmp(s.ctxId,ctx.id),'ipm_cr:StaleState', ...
    'State belongs to another context. Use refresh only for an identical mathematical setup.');
if nargin<4,controls=struct();end
[s,info]=ipm_cr.time.step(s,ctx,tauEnd,controls);
end
