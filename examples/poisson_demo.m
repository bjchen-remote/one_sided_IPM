function [solver,solution,info]=poisson_demo(options)
%POISSON_DEMO Reproducible mixed Gaussian elliptic demonstration (no evolution).
%   addpath('src','examples'); [solver,solution,info]=poisson_demo;
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,options=struct();end
options.outerBoundary='dirichlet';solver=ipm_cr.psolver(options);
[solution,info]=solver.solve(@source,@boundary);
fprintf('ipm_cr: build %.3fs, solve %.3fs, overlap psi %.3g, u %.3g, strong %.3g\n', ...
    info.buildSeconds,info.solveSeconds,info.overlapValueMismatch, ...
    info.overlapVelocityMismatch,info.strongResidual);
end

function q=source(x,y)
beta=1.7;P=y+.2*x.*y+.3*y.^2+.1*x.^2.*y;
Px=.2*y+.2*x.*y;Py=1+.2*x+.6*y+.1*x.^2;
q=-exp(-x.^2-beta*y.^2).*(.2*y+.6-4*(x.*Px+beta*y.*Py)+ ...
    (4*x.^2+4*beta^2*y.^2-2-2*beta).*P);
end

function psi=boundary(x,y)
psi=(y+.2*x.*y+.3*y.^2+.1*x.^2.*y).*exp(-x.^2-1.7*y.^2);
end
