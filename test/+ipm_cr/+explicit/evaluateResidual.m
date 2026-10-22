function [flow,rho]=evaluateResidual(residual,scale,g,wallOps,options,guess,wantDiagnostics)
%EVALUATERESIDUAL Full IPM RHS in a prescribed analytic background.
% The evolving array is eta=W-(G-G(0)); rho returned here remains full W.
% Each stage supplies its own scale. Only eta ghosts are interpolated.
% Background timeH/timeX already include dt/dtau=exp(logC_omega-logC_l).
% evaluate subtracts timeH from eta_t and restores it in the full-W rate;
% source derivatives and Poisson always use the same complete background.
if nargin<6,guess=[0,0];end
if nargin<7,wantDiagnostics=true;end
model=options.background;
context=struct('residual',residual,'scale',scale, ...
    'core',ipm_cr.explicit.background(model,scale,g.core.X,g.core.Y,wantDiagnostics), ...
    'polar',ipm_cr.explicit.background(model,scale,g.polar.X,g.polar.Y,wantDiagnostics), ...
    'source',@(X,Y)ipm_cr.explicit.backgroundSource(model,scale,X,Y));
[flow,rho]=ipm_cr.explicit.evaluate([],g,wallOps,options,guess,context,wantDiagnostics);
end
