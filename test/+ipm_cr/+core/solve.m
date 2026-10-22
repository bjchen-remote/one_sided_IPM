function s=solve(c,source,trace,refineOriginalPde)
%SOLVE -Delta psi=source, zero physical-wall and supplied artificial trace.
% Optional refinement is requested only by the owning Solver's original
% direct-fallback path. Ordinary cached solves retain their existing cost.
% Corrections solve the ORIGINAL PDE residual with homogeneous boundaries
% and the same cached Sylvester factors. No acceptance gate is changed.
if nargin<4,refineOriginalPde=false;end
psi=zeros(c.ny,c.nx);psi(c.traceIndex)=trace;
boundaryLaplacian=-(c.Dyy*psi+psi*c.Dxx');
rhs=source(c.iy,c.ix)-boundaryLaplacian(c.iy,c.ix);
psi(c.iy,c.ix)=c.Vy*((c.Vyi*rhs*c.Vxi.')./c.denominator)*c.Vx.';
s=struct('psi',psi);
if ~refineOriginalPde,return;end
% This roundoff guard only controls whether to attempt a correction. The
% public linearTolerance and all owning-Solver acceptance tests are intact.
guard=1e-12;maximumAttempts=3;
[backward,residual]=original_backward(c,psi,source);
initialBackward=backward;attempts=0;accepted=0;candidateErrors=zeros(1,0);
for k=1:maximumAttempts
    if ~isfinite(backward)||backward<=guard,break;end
    attempts=attempts+1;
    correctionRhs=-residual(c.iy,c.ix);
    correction=c.Vy*((c.Vyi*correctionRhs*c.Vxi.')./c.denominator)*c.Vx.';
    candidate=psi;
    candidate(c.iy,c.ix)=psi(c.iy,c.ix)+correction;
    [candidateBackward,candidateResidual]=original_backward(c,candidate,source);
    candidateErrors(end+1)=candidateBackward; %#ok<AGROW>
    if ~isfinite(candidateBackward)||candidateBackward>=backward,break;end
    psi=candidate;backward=candidateBackward;residual=candidateResidual;
    accepted=accepted+1;
end
s.psi=psi;
s.refinement=struct('requested',true,'originalPdeComponentwiseMetric',true, ...
    'triggerRoundoffGuard',guard,'maximumAttempts',maximumAttempts, ...
    'attempts',attempts,'acceptedCorrections',accepted, ...
    'initialBackwardError',initialBackward,'finalBackwardError',backward, ...
    'candidateBackwardErrors',candidateErrors,'homogeneousCorrectionBoundary',true, ...
    'cachedSylvesterFactorsReused',true,'acceptanceToleranceChanged',false);
end

function [backward,residual]=original_backward(c,psi,source)
% Keep the original owning-Solver audit's arithmetic and denominator.
residual=-c.Dyy*psi-psi*c.Dxx'-source;
denominator=abs(c.Dyy)*abs(psi)+abs(psi)*abs(c.Dxx')+abs(source);
backward=max(abs(residual(c.iy,c.ix))./ ...
    max(denominator(c.iy,c.ix),realmin),[],'all');
end
