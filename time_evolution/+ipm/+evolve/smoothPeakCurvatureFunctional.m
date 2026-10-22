function [nodes,K,J,curvatureWeights,conditioning] = ...
        smoothPeakCurvatureFunctional(q,x,peak,radius)
%IPM.EVOLVE.SMOOTHPEAKCURVATUREFUNCTIONAL Curvature in a fixed wall window.
%   [NODES,K,J,WEIGHTS,CONDITIONING] = ... (Q,X,PEAK,RADIUS) fits a degree-4
%   polynomial to the wall-gradient Q at all nodes with |X|<=RADIUS. The
%   dimensionless fit coordinate is xi=X/RADIUS and the fixed Gaussian
%   weights are exp(-0.5*(xi/0.5)^2). PEAK is the current output of
%   quadraticPeakFunctional and supplies the evaluation location a.
%   K=q_fit''(a)<0 and J=q_fit'''(a). WEIGHTS is a row vector satisfying
%   K=WEIGHTS*Q(NODES)' to floating-point accuracy. While the fixed window
%   is unchanged, the moving-peak derivative is
%       K_tau=WEIGHTS*Q_tau(NODES)' + J*a_tau.

validateattributes(q,{'single','double'}, ...
    {'vector','real','finite','nonempty'},mfilename,'q',1);
validateattributes(x,{'single','double'}, ...
    {'vector','real','finite','nonempty'},mfilename,'x',2);
validateattributes(radius,{'single','double'}, ...
    {'scalar','real','finite','positive'},mfilename,'radius',4);
q = q(:);
x = x(:);
if numel(q) ~= numel(x) || numel(x) < 9 || any(diff(x) <= 0)
    error('ipm:SmoothPeakCurvatureInput', ...
        'Q and X must have equal lengths, at least nine nodes, and ordered X.');
end
if ~isstruct(peak) || ~isscalar(peak) || ...
        ~isfield(peak,'x') || ~isnumeric(peak.x) || ...
        ~isreal(peak.x) || ~isscalar(peak.x) || ...
        ~isfinite(peak.x)
    error('ipm:SmoothPeakCurvaturePeak', ...
        'PEAK must contain one finite quadratic vertex PEAK.x.');
end
a = double(peak.x);
radius = double(radius);
if abs(a) >= radius
    error('ipm:SmoothPeakCurvatureWindow', ...
        'The quadratic wall-gradient peak must lie inside |X|<RADIUS.');
end
nodes = find(abs(x) <= radius);
if numel(nodes) < 9
    error('ipm:SmoothPeakCurvatureWindow', ...
        'The fixed curvature window must contain at least nine nodes.');
end

xi = double(x(nodes))/radius;
qWindow = double(q(nodes));
sqrtWeight = exp(-0.25*(xi/0.5).^2);
design = [ones(size(xi)),xi,xi.^2,xi.^3,xi.^4];
weightedDesign = design.*sqrtWeight;
[orthogonal,triangular] = qr(weightedDesign,0);
designRcond = rcond(triangular);
if ~isfinite(designRcond) || designRcond < 1e-10
    error('ipm:SmoothPeakCurvatureCondition', ...
        'The fixed-window degree-4 curvature fit is ill conditioned.');
end
coefficients = triangular\(orthogonal'*(sqrtWeight.*qWindow));
xiA = a/radius;
secondDerivative = [0,0,2,6*xiA,12*xiA^2]/radius^2;
thirdDerivative = [0,0,0,6,24*xiA]/radius^3;
K = secondDerivative*coefficients;
J = thirdDerivative*coefficients;
curvatureWeights = ((secondDerivative/triangular)* ...
    orthogonal').*sqrtWeight';
curvatureCondition = abs(K)/ ...
    max(sum(abs(curvatureWeights).*abs(qWindow')),realmin);
relativeFitResidual = norm( ...
    sqrtWeight.*(design*coefficients-qWindow)) / ...
    max(norm(sqrtWeight.*qWindow),realmin);
if ~isfinite(K) || ~isfinite(J) || K >= 0 || ...
        ~isfinite(curvatureCondition) || curvatureCondition < 1e-8
    error('ipm:SmoothPeakCurvatureDegenerate', ...
        'The fixed-window fit has no resolved concave peak curvature.');
end
conditioning = struct('designRcond',designRcond, ...
    'curvatureCondition',curvatureCondition, ...
    'relativeFitResidual',relativeFitResidual, ...
    'nodeCount',numel(nodes), ...
    'windowMargin',1-abs(a)/radius);
end
