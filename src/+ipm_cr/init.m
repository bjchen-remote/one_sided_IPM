function [s,ctx]=init(opts)
%INIT Original datum and its own prepared simulation context.
if nargin<1,opts=struct();end
[s,~,~,~,ctx]=ipm_cr.evolve.init(opts);
s=ipm_cr.internal.state(s,'new');
end
