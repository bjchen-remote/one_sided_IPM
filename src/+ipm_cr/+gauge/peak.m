function peak = peak(values,ops)
%IPM.EVOLVE.GLOBALSPLINEMAXIMUM Global maximum of a complete-domain C2 spline.
% Natural boundary conditions are imposed at the actual first/last x nodes.
% Every node, endpoint, and real stationary point is considered. There is
% no tracked branch, curvature/positivity floor, or phase/time restriction.
% Flat maxima and exchanges of the global maximizing branch are allowed.
% The height of the maximum is a continuous observable of the nodal data.
% At a tie its derivative need not be unique: returned weights evaluate the
% selected maximizer, not a unique derivative of the maximum functional.
% Numerical ties (within peak-scaled roundoff) prefer the location nearest
% pinX, then the smaller coordinate; the numerical selection gap is exposed.

assert(isstruct(ops) && isscalar(ops) && isfield(ops,'x'), ...
    'ipm:GlobalSplineInput','Supply a scalar operator structure with x.');
x = ops.x(:)'; q = values(:)';
assert(isnumeric(values) && isreal(values) && isnumeric(ops.x) && ...
    isreal(ops.x) && numel(q) == numel(x) && numel(x) >= 2 && ...
    all(isfinite(q)) && all(isfinite(x)) && all(diff(x) > 0), ...
    'ipm:GlobalSplineInput','The spline needs finite real values on increasing x.');
pinX = 0;
if isfield(ops,'rescaling') && isfield(ops.rescaling,'pinX')
    pinX = ops.rescaling.pinX;
end
assert(isnumeric(pinX) && isreal(pinX) && isscalar(pinX) && isfinite(pinX), ...
    'ipm:GlobalSplineInput','The optional tie-breaking pinX must be finite.');
g = spline_geometry(x);
M = zeros(numel(x),1);
if numel(x) > 2
    M(2:end-1) = g.L'\(g.L\(g.B*q'));
end
% Rows: physical coordinate, interval, normalized coordinate, value, flat.
n = numel(x);
candidates = [x(:),[1:n-1,n-1]',[zeros(1,n-1),1]',q(:),zeros(n,1)];
constant = all(q == q(1));
for i = 1:n-1
    H = g.h(i); m0 = M(i); m1 = M(i+1);
    c0 = q(i+1)-q(i)+H^2*(-2*m0-m1)/6;
    c1 = H^2*m0; c2 = H^2*(m1-m0)/2;
    if constant || all([c2,c1,c0] == 0)
        a = min(x(i+1),max(x(i),pinX));
        t = (a-x(i))/H;
        candidates(end+1,:) = [a,i,t,q(i),1]; %#ok<AGROW>
        continue
    end
    critical = real_quadratic_roots(c2,c1,c0);
    for j = 1:numel(critical)
        t = critical(j);
        % Physical-coordinate rounding can exceed eps(t) when a root lies
        % on a source node. Retain the previous endpoint tolerance repair.
        tolerance = max(128*eps(max(1,abs(t))), ...
            128*eps(max([1,abs(x(i)),abs(x(i+1))]))/H);
        if t < -tolerance || t > 1+tolerance, continue; end
        t = min(max(t,0),1);
        a = x(i)+H*t;
        if t == 0
            value = q(i);
        elseif t == 1
            value = q(i+1);
        else
            value = (1-t)*q(i)+t*q(i+1)+H^2/6* ...
                (((1-t)^3-(1-t))*m0+(t^3-t)*m1);
        end
        candidates(end+1,:) = [a,i,t,value,0]; %#ok<AGROW>
    end
end
maximumCandidateValue = max(candidates(:,4));
tieTolerance = 64*eps(max(abs(maximumCandidateValue),realmin));
tied = candidates(candidates(:,4) >= maximumCandidateValue-tieTolerance,:);
% This order makes exact/numerical tie selection independent of root order.
[~,order] = sortrows([abs(tied(:,1)-pinX),tied(:,1)],[1,2]);
selected = tied(order(1),:);
a = selected(1); i = selected(2); t = selected(3);
[w,wx,wxx] = functional_rows(g,i,t);
indices = find(w ~= 0 | wx ~= 0 | wxx ~= 0);
value = w(indices)*q(indices)';
first = wx(indices)*q(indices)';
second = wxx(indices)*q(indices)';
if constant
    % These identities are exact for a constant natural spline. Avoid
    % presenting cancellation in a derivative row as physical curvature.
    value = q(1); first = 0; second = 0;
end
coordinateTolerance = 256*eps(max([1,abs(a)]));
tieNonunique = any(abs(tied(:,1)-a) > coordinateTolerance) || ...
    any(tied(:,5) ~= 0);
[~,node] = min(abs(x-a));
peak = struct('definition','global_natural_cubic_spline_maximum', ...
    'index',node,'indices',indices,'weights',w(indices), ...
    'firstWeights',wx(indices),'secondWeights',wxx(indices), ...
    'x',a,'value',value,'curvature',second/2, ...
    'stationarityResidual',first,'tieNonunique',logical(tieNonunique), ...
    'valueDerivativeUnique',~tieNonunique, ...
    'derivativeInterpretation','selected_maximizer_evaluation', ...
    'tieTolerance',tieTolerance,'maximumCandidateValue',maximumCandidateValue, ...
    'selectionGap',maximumCandidateValue-value, ...
    'window',x([1,end]),'cellIndex',i,'cellFraction',t, ...
    'gridValue',q(node),'weightCondition',1/sum(abs(w(indices))), ...
    'splineKnots',x,'splineSecondDerivatives',M);
end

function g = spline_geometry(x)
persistent cached
if ~isempty(cached) && isequal(cached.x,x), g = cached; return; end
n = numel(x); h = diff(x);
A = spalloc(n-2,n-2,max(0,3*(n-2)));
B = spalloc(n-2,n,max(0,3*(n-2)));
for j = 2:n-1
    row = j-1; hl = h(j-1); hr = h(j);
    A(row,row) = 2*(hl+hr);
    if row > 1, A(row,row-1) = hl; end
    if row < n-2, A(row,row+1) = hr; end
    B(row,j-1:j+1) = [6/hl,-6*(1/hl+1/hr),6/hr];
end
if n > 2, L = chol(A,'lower'); else, L = sparse(0,0); end
g = struct('x',x,'h',h,'B',B,'L',L);
cached = g;
end

function roots = real_quadratic_roots(a,b,c)
scale = max(abs([a,b,c])); roots = [];
if scale == 0, return; end
a = a/scale; b = b/scale; c = c/scale;
if abs(a) <= 64*eps(1)
    if abs(b) > 64*eps(1), roots = -c/b; end
    return
end
discriminant = b*b-4*a*c;
if discriminant < 0
    if discriminant >= -64*eps(1)*(b*b+4*abs(a*c))
        discriminant = 0;
    else
        return
    end
end
signed = 1; if b < 0, signed = -1; end
v = -.5*(b+signed*sqrt(discriminant));
if v == 0, roots = -b/(2*a); else, roots = [v/a,c/v]; end
end

function [w,wx,wxx] = functional_rows(g,index,t)
n = numel(g.x); H = g.h(index);
r = zeros(n,1); r(index:index+1) = [1-t;t];
c = zeros(n,1); c(index:index+1) = H^2/6*[(1-t)^3-(1-t);t^3-t];
w = transpose_row(g,r,c);
r(:) = 0; r(index:index+1) = [-1/H;1/H];
c(:) = 0; c(index:index+1) = H/6*[-3*(1-t)^2+1;3*t^2-1];
wx = transpose_row(g,r,c);
r(:) = 0; c(:) = 0; c(index:index+1) = [1-t;t];
wxx = transpose_row(g,r,c);
end

function row = transpose_row(g,r,c)
if numel(g.x) > 2
    row = (r+g.B'*(g.L'\(g.L\c(2:end-1))))';
else
    row = r';
end
end
