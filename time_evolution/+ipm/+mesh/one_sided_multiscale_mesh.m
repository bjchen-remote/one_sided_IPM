function [x,y,audit] = one_sided_multiscale_mesh(nx,ny,left,right,ymax,frontX)
%ONE_SIDED_MULTISCALE_MESH Smooth full-X mesh for a million-scale box.
%   [X,Y,AUDIT] = ONE_SIDED_MULTISCALE_MESH() returns row X and column Y
%   for [-1e6,1e6] x [0,1e6], with 513 x 193 nodes and an exact X=0 node.
%   Optional arguments are NX,NY,LEFT,RIGHT,YMAX,FRONTX, in that order.
%   This is an example coordinate builder; it does not change the solver.
%
%   The horizontal monitor has an analytic cumulative integral. Its asinh
%   tail yields geometric cell growth over many decades; two Gaussians
%   refine X=0 and the one-sided support front. The two integer cell counts
%   nearest the global monitor mass split make X=0 exact, with only a small
%   difference between the left and right mass steps.

if nargin < 1 || isempty(nx), nx = 513; end
if nargin < 2 || isempty(ny), ny = 193; end
if nargin < 3 || isempty(left), left = 1e6; end
if nargin < 4 || isempty(right), right = 1e6; end
if nargin < 5 || isempty(ymax), ymax = 1e6; end
if nargin < 6 || isempty(frontX), frontX = -7^(1/8); end
assert(isscalar(nx) && isfinite(nx) && nx==fix(nx) && nx>=9 && ...
    isscalar(ny) && isfinite(ny) && ny==fix(ny) && ny>=9 && ...
    isscalar(left) && isfinite(left) && left>0 && ...
    isscalar(right) && isfinite(right) && right>0 && ...
    isscalar(ymax) && isfinite(ymax) && ymax>0 && ...
    isscalar(frontX) && isfinite(frontX) && -left<frontX && frontX<0, ...
    'Invalid one-sided multiscale mesh parameters.');

% Monitor: 1/sqrt(x^2+a^2) plus two Gaussian refinements. The chosen
% scales give about 0.02 at X=0 and 0.03 at FRONTX for the default box.
coreScale = 0.4;
originWeight = 0.5;
originSigma = 0.55;
frontWeight = 1.4;
frontSigma = 0.6;
leftEdge = -left;
massAtZero = mass(0);
totalMass = mass(right);
numberOfCells = nx-1;
leftCells = round(numberOfCells*massAtZero/totalMass);
leftCells = max(1,min(numberOfCells-1,leftCells));
rightCells = numberOfCells-leftCells;

% Monotone bisection of the analytic CDF; no huge sampled integration grid.
leftTargets = (0:leftCells)*(massAtZero/leftCells);
rightTargets = massAtZero+(1:rightCells)* ...
    ((totalMass-massAtZero)/rightCells);
targets = [leftTargets,rightTargets];
lower = leftEdge*ones(1,nx);
upper = right*ones(1,nx);
% A fixed 64 iterations leave O(1e-9) absolute coordinate uncertainty on
% a 1e10 box. Resolve local O(1) nodes to floating-point precision instead.
for iteration = 1:128
    middle = lower+(upper-lower)/2;
    if all(middle == lower | middle == upper)
        break
    end
    below = mass(middle)<targets;
    lower(below) = middle(below);
    upper(~below) = middle(~below);
end
x = lower+(upper-lower)/2;
x([1,leftCells+1,end]) = [leftEdge,0,right];

% The asinh monitor on Y integrates to the exact sinh quantile map.
verticalScale = 0.375;
verticalMass = asinh(ymax/verticalScale);
parameter = (0:ny-1)'/(ny-1);
y = verticalScale*sinh(verticalMass*parameter);
y([1,end]) = [0,ymax];

dx = diff(x);
dy = diff(y);
assert(all(isfinite(x)) && all(isfinite(y)) && all(dx>0) && ...
    all(dy>0) && x(leftCells+1)==0, ...
    'One-sided multiscale mesh is not strictly increasing.');
frontCell = find(x>=frontX,1)-1;
assert(frontCell>=1 && frontCell<=numberOfCells, ...
    'The front is outside the horizontal mesh.');
audit = struct( ...
    'zeroIndex',leftCells+1, ...
    'leftCells',leftCells,'rightCells',rightCells, ...
    'originLeftDx',dx(leftCells), ...
    'originRightDx',dx(leftCells+1), ...
    'frontDx',dx(frontCell), ...
    'wallDy',dy(1), ...
    'maximumAdjacentXRatio',max([dx(2:end)./dx(1:end-1), ...
        dx(1:end-1)./dx(2:end)],[],'all'), ...
    'maximumAdjacentYRatio',max([dy(2:end)./dy(1:end-1), ...
        dy(1:end-1)./dy(2:end)],[],'all'), ...
    'massStepRatio',max((massAtZero/leftCells)/ ...
        ((totalMass-massAtZero)/rightCells), ...
        ((totalMass-massAtZero)/rightCells)/(massAtZero/leftCells)));
assert(audit.maximumAdjacentXRatio<1.3 && ...
    audit.maximumAdjacentYRatio<1.3, ...
    'The requested box or node counts exceed the adjacent-cell ratio limit.');

    function value = mass(location)
        value = asinh(location/coreScale)-asinh(leftEdge/coreScale);
        value = value+originWeight*originSigma*sqrt(pi/2)* ...
            (erf(location/(sqrt(2)*originSigma))- ...
             erf(leftEdge/(sqrt(2)*originSigma)));
        value = value+frontWeight*frontSigma*sqrt(pi/2)* ...
            (erf((location-frontX)/(sqrt(2)*frontSigma))- ...
             erf((leftEdge-frontX)/(sqrt(2)*frontSigma)));
    end
end
