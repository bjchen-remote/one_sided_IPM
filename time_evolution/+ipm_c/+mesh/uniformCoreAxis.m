function [r,info]=uniformCoreAxis(options)
%UNIFORMCOREAXIS Smooth odd linear core with an analytic stretched tail.
% r(s)=a*s+(R-a)*(sinh(alpha*s)/sinh(alpha))^3, 0<=s<=1.
% firstSpacing is a lower bound: increase it when the measured core would
% occupy more than 60 percent of the logical axis. This retains enough tail
% intervals to reach rmax smoothly instead of exhausting them as the core
% grows. Set a so the first positive node equals the selected spacing.
% Choose alpha so the tail at s=coreRadius/((nr-1)*spacing) is only the
% declared fraction of coreRadius. Both positive terms increase and are odd.
nr=option(options,'nr',129);R=option(options,'rmax',1e10);
requestedSpacing=option(options,'firstSpacing',.02);
core=option(options,'coreRadius',2);fraction=option(options,'coreTailFraction',.01);
validateattributes(nr,{'numeric'},{'scalar','real','finite','integer','>=',7});
validateattributes(R,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(requestedSpacing,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(core,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(fraction,{'numeric'},{'scalar','real','finite','>',0,'<',1});
N=nr-1;maximumCoreLogicalLocation=.6;
h=max(requestedSpacing,core/(maximumCoreLogicalLocation*N));
sFirst=1/N;sCore=core/(N*h);power=3;
assert(N*h<R&&sCore<1&&core<R,'polar:UniformCoreRange', ...
    'The uniform core must fit inside the logical axis and physical radius.');
targetTail=core*fraction;
lower=0;upper=1;
while tail_at_core(upper)>targetTail
    upper=2*upper;
    assert(isfinite(upper),'polar:UniformCoreParameter','The tail parameter must be finite.');
end
assert(tail_at_core(lower)>=targetTail,'polar:UniformCoreParameter', ...
    'The requested tail fraction lies outside the analytic map family.');
for iteration=1:100
    middle=lower+(upper-lower)/2;
    if tail_at_core(middle)>targetTail,lower=middle;else,upper=middle;end
end
alpha=lower+(upper-lower)/2;
a=linear_coefficient(alpha);tailCoefficient=R-a;
assert(isfinite(a)&&a>0&&isfinite(tailCoefficient)&&tailCoefficient>0, ...
    'polar:UniformCoreParameter','The core and tail coefficients must both be positive.');
s=(0:N)'/N;r=a*s+tailCoefficient*tail_ratio(s,alpha);
r([1,2,end])=[0,h,R];
assert(all(isfinite(r))&&all(diff(r)>0)&&abs(r(2)-h)<=64*eps(h), ...
    'polar:UniformCoreRepresentability','The analytic map must be increasing and retain its first spacing.');
gaps=diff(r);inside=r(2:end)<=core;
adjacentGapRatio=max([gaps(2:end)./gaps(1:end-1);gaps(1:end-1)./gaps(2:end)]);
assert(any(inside),'polar:UniformCoreRange','The core must contain a positive radial node.');
coreGaps=gaps(inside);observedTail=tail_at_core(alpha);
info=struct('kind','positive_linear_plus_cubic_sinh_odd_axis_v2', ...
    'nr',nr,'rmax',R,'requestedFirstSpacing',requestedSpacing, ...
    'firstSpacing',h,'firstSpacingAdjusted',h>requestedSpacing, ...
    'a',a,'alpha',alpha,'power',power, ...
    'tailCoefficient',tailCoefficient,'coreRadius',core,'coreTailFraction',fraction, ...
    'coreLogicalLocation',sCore,'maximumCoreLogicalLocation',maximumCoreLogicalLocation, ...
    'positiveCoreNodeFraction',nnz(inside)/N, ...
    'tailIntervalCount',N-nnz(inside),'maximumAdjacentGapRatio',adjacentGapRatio, ...
    'tailAtCoreLogicalLocation',observedTail, ...
    'targetTailAtCoreLogicalLocation',targetTail, ...
    'coreTailRelativeResidual',abs(observedTail-targetTail)/targetTail, ...
    'positiveCoreNodeCount',nnz(inside),'coreGaps',coreGaps, ...
    'minimumCoreGap',min(coreGaps),'maximumCoreGap',max(coreGaps), ...
    'coreGapRatio',max(coreGaps)/min(coreGaps), ...
    'signedInverseMapOddAndSmooth',true,'originAndRadiusRetained',true, ...
    'mapDefinition','a*s+(rmax-a)*(sinh(alpha*s)/sinh(alpha))^3');

    function value=linear_coefficient(parameter)
        eta=tail_ratio(sFirst,parameter);
        value=(h-R*eta)/(sFirst-eta);
    end
    function value=tail_at_core(parameter)
        value=(R-linear_coefficient(parameter))*tail_ratio(sCore,parameter);
    end
    function value=tail_ratio(location,parameter)
        if parameter==0
            value=location.^power;
        else
            value=zeros(size(location));positive=location>0;
            value(positive)=exp(power*(log_sinh(parameter*location(positive))-log_sinh(parameter)));
        end
    end
end

function value=log_sinh(argument)
value=argument+log(-expm1(-2*argument))-log(2);
end

function value=option(options,name,fallback)
if isfield(options,name),value=options.(name);else,value=fallback;end
end
