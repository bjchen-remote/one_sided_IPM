function [summary,state,g,wallOps]=run_ipm_cr_explicit(options)
%RUN_IPM_CR_EXPLICIT Compatibility wrapper for historical research callers.
% New applications use [summary,state,ctx]=ipm_cr.run(options).
if nargin<1,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
[summary,state,ctx]=ipm_cr.run(options);g=[];wallOps=[];
if ~isempty(ctx),g=ctx.mesh;wallOps=ctx.wall;end
end
