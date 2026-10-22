function [flow,rho]=residual(residual,scale,g,wallOps,options,guess,wantDiagnostics,ctx)
%EVALUATERESIDUAL Full IPM RHS in a prescribed analytic background.
% The evolving array is eta=W-(G-G(0)); rho returned here remains full W.
% Each stage supplies its own scale. Only eta ghosts are interpolated.
% Background timeH/timeX already include dt/dtau=exp(logC_omega-logC_l).
% evaluate subtracts timeH from eta_t and restores it in the full-W rate;
% source derivatives and Poisson always use the same complete background.
if nargin<6,guess=[0,0];end
if nargin<7,wantDiagnostics=true;end
if nargin<8||isempty(ctx),ctx=ipm_cr.context(g,wallOps,options);end
model=options.background;
context=struct('residual',residual,'scale',scale, ...
    'core',ipm_cr.model.background(model,scale,g.core.X,g.core.Y,wantDiagnostics), ...
    'polar',ipm_cr.model.background(model,scale,g.polar.X,g.polar.Y,wantDiagnostics), ...
    'source',@(X,Y)ipm_cr.model.source(model,scale,X,Y));
[flow,rho]=ipm_cr.evolve.rhs([],g,wallOps,options,guess,context,wantDiagnostics,ctx);
end
