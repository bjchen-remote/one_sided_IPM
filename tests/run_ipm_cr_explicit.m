function [summary,state,g,wallOps]=run_ipm_cr_explicit(options)
%RUN_IPM_CR_EXPLICIT Compatibility wrapper for historical research callers.
% New applications use [summary,state,ctx]=ipm_cr.run(options).
if nargin<1,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
[summary,public,ctx]=ipm_cr.run(options);state=[];g=[];wallOps=[];
if ~isempty(public),state=ipm_cr.internal.state(public,'old');end
if ~isempty(ctx),g=ctx.mesh;wallOps=ctx.wall;end
end
