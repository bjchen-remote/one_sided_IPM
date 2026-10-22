function solution = solve(p, source, inner, outer)
%SOLVE Solve -Delta psi=source using the cached polar factorizations.
%   INNER is the trace at rInner. OUTER is a Dirichlet trace, or the right
%   hand side of psi_t+sqrt(-d_theta_theta)psi for outerBoundary='decay'.
%   Both physical walls have zero streamfunction, independently of SOURCE.

if nargin < 3 || isempty(inner), inner = zeros(1,p.nt); end
if nargin < 4 || isempty(outer), outer = zeros(1,p.nt); end
if isa(source,'function_handle'), source = source(p.X,p.Y); end
if isscalar(source), source = repmat(source,p.nr,p.nt); end
if ~isequal(size(source), [p.nr,p.nt])
    error('ipm_cr:polar:SourceSize', ...
        'source must have size [%d,%d].',p.nr,p.nt);
end
inner = traceRow(inner,p.nt,'inner');
outer = traceRow(outer,p.nt,'outer');
innerModes = inner(p.interior)*p.Vinv.';
outerModes = outer(p.interior)*p.Vinv.';
sourceModes = source(:,p.interior)*p.Vinv.';
rhs = zeros(p.nr,numel(p.mu));
rhs(p.equationRows,:) = bsxfun(@times, ...
    sourceModes(p.equationRows,:),p.r(p.equationRows).^2);
rhs(1,:) = innerModes;
rhs(end,:) = outerModes;
modal = zeros(size(rhs));
roundoffGuard = 128*eps;
initialBackward = zeros(1,numel(p.mu));
finalBackward = initialBackward;
refinementAttempts = initialBackward;
refinementUpdates = initialBackward;
for k = 1:numel(p.mu)
    modal(:,k) = p.factors{k}\(p.rowScale(:,k).*rhs(:,k));
    % The first ORIGINAL modal row is exactly the inner Dirichlet value.
    % Enforce it before auditing, including on every correction candidate.
    modal(1,k) = innerModes(k);
    % Use the owning Solver's original componentwise backward error. An
    % unscaled residual norm overweights the large PDE/flux rows and can
    % spend several solves improving an already roundoff-accurate solution.
    A = p.operators{k};
    absA = p.absOperators{k};
    residual = rhs(:,k)-A*modal(:,k);
    backward = max(abs(residual)./max(absA*abs(modal(:,k))+abs(rhs(:,k)),realmin));
    initialBackward(k) = backward;
    for attempt = 1:3
        if backward <= roundoffGuard, break; end
        refinementAttempts(k) = refinementAttempts(k)+1;
        candidate = modal(:,k) ...
            + p.factors{k}\(p.rowScale(:,k).*residual);
        candidate(1) = innerModes(k);
        candidateResidual = rhs(:,k)-A*candidate;
        candidateBackward = max(abs(candidateResidual)./ ...
            max(absA*abs(candidate)+abs(rhs(:,k)),realmin));
        if ~isfinite(candidateBackward) || candidateBackward >= backward, break; end
        modal(:,k) = candidate;
        residual = candidateResidual;
        backward = candidateBackward;
        refinementUpdates(k) = refinementUpdates(k)+1;
    end
    finalBackward(k) = backward;
end
psi = zeros(p.nr,p.nt);
psi(:,p.interior) = modal*p.V.';
solution = struct('psi',psi,'modal',modal,'rhs',rhs);
solution.refinement = struct('scope','this_direct_modal_rhs_before_any_later_trace_update', ...
    'roundoffGuard',roundoffGuard, ...
    'initialBackwardError',initialBackward,'finalBackwardError',finalBackward, ...
    'attempts',refinementAttempts,'acceptedCorrections',refinementUpdates, ...
    'exactModalInnerBoundary',true,'acceptanceToleranceChanged',false);
end

function value = traceRow(value, nt, name)
if isscalar(value), value = repmat(value,1,nt); end
if numel(value) ~= nt
    error('ipm_cr:polar:TraceSize','%s must have %d values.',name,nt);
end
value = reshape(value,1,nt);
if any(abs(value([1,end])) > 1e-11*max(1,max(abs(value))))
    error('ipm_cr:polar:WallTrace', ...
        '%s trace must vanish at the physical walls theta=0,pi.',name);
end
value([1,end]) = 0;
end
