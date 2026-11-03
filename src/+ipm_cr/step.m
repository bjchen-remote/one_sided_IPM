function [s,info]=step(s,ctx,tauEnd,controls)
%STEP One direct explicit step on the supplied prepared context.
if nargin<4,controls=struct();end
[s,info]=ipm_cr.time.step(s,ctx,tauEnd,controls);
end
