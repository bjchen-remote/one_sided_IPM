function axes=tensorAxis(a,b,input)
%TENSORAXIS Smooth sinh tensor axes with prescribed origin metric spacing.
% x=a*sinh(alpha*xi)/sinh(alpha), xi in [-1,1];
% y=b*sinh(beta*eta)/sinh(beta), eta in [0,1].
% The requested spacing is J(0)*dlogical, not a secant approximation.
% coreMap contains only constructor inputs and can be passed back to BUILD;
% derived map parameters are stored separately. No PDE state is changed.
if nargin<3,input=struct();end
validateattributes(a,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(b,{'numeric'},{'scalar','real','finite','positive'});
map=struct('enabled',true,'kind','sinh','intervals',[240,128], ...
    'minSpacing',[.01,.00375]);
assert(isstruct(input)&&isscalar(input)&&all(ismember(fieldnames(input),fieldnames(map))), ...
    'ipm_cr:TensorMap','Unknown tensor map input.');
names=fieldnames(input);
for k=1:numel(names),map.(names{k})=input.(names{k});end
validateattributes(map.enabled,{'logical'},{'scalar'});
assert(map.enabled,'ipm_cr:TensorMap','Disabled mappings use the unchanged uniform BUILD path.');
map.kind=char(validatestring(map.kind,{'sinh'}));
validateattributes(map.intervals,{'numeric'},{'vector','numel',2,'real','finite','integer','>=',8});
validateattributes(map.minSpacing,{'numeric'},{'vector','numel',2,'real','finite','positive'});
map.intervals=reshape(map.intervals,1,2);map.minSpacing=reshape(map.minSpacing,1,2);
assert(mod(map.intervals(1),2)==0,'ipm_cr:TensorMapOrigin', ...
    'An even number of x intervals is required for the exact origin.');
[x,Jx,dx,dxi,xi,alpha]=axis_map(a,map.intervals(1),map.minSpacing(1),true);
[y,Jy,dy,deta,eta,beta]=axis_map(b,map.intervals(2),map.minSpacing(2),false);
parameters=struct('alpha',alpha,'beta',beta,'uniformAxes',[alpha==0,beta==0], ...
    'originMetricSpacing',[min(dx),min(dy)], ...
    'averagePhysicalSpacing',[2*a/map.intervals(1),b/map.intervals(2)], ...
    'logicalSpacing',[dxi,deta],'logicalBounds',[-1,1,0,1], ...
    'definition','analytic_J_times_logical_spacing');
axes=struct('x',x,'y',y,'Jx',Jx,'Jy',Jy,'dxMetric',dx,'dyMetric',dy, ...
    'dxi',dxi,'deta',deta,'xi',xi,'eta',eta,'coreMap',map,'mapParameters',parameters);
end

function [physical,J,metric,h,logical,s]=axis_map(extent,n,requested,symmetric)
span=1+double(symmetric);average=span*extent/n;h=span/n;
if symmetric,logical=linspace(-1,1,n+1);else,logical=linspace(0,1,n+1)';end
if average<=requested
    s=0;
    % Match the original physical LINSPACE exactly when the map degenerates.
    if symmetric,physical=linspace(-extent,extent,n+1);
    else,physical=linspace(0,extent,n+1)';end
    J=extent*ones(size(logical));
else
    target=log(average)-log(requested);lo=0;hi=1;
    while log_sinhc(hi)<target,hi=2*hi;end
    for k=1:100
        mid=(lo+hi)/2;
        if mid==lo||mid==hi,break;end
        if log_sinhc(mid)<target,lo=mid;else,hi=mid;end
    end
    s=(lo+hi)/2;
    if s<40
        physical=extent*sinh(s*logical)/sinh(s);
        J=(extent*s/sinh(s))*cosh(s*logical);
    else
        % Avoid overflow for valid positive spacing ratios with large alpha.
        q=abs(logical);denominator=-expm1(-2*s);
        positive=exp(s*(q-1));negative=exp(-s*(q+1));
        physical=extent*sign(logical).*(positive-negative)/denominator;
        J=(extent*s/denominator)*(positive+negative);
    end
end
if symmetric
    physical([1,end])=[-extent,extent];physical(n/2+1)=0;
else
    physical(1)=0;physical(end)=extent;
end
metric=J*h;
assert(all(isfinite(physical))&&all(diff(physical)>0)&&all(isfinite(J))&&all(J>0)&& ...
    all(isfinite(metric))&&all(metric>0),'ipm_cr:TensorMapRange', ...
    'The mapped axis or analytic Jacobian is not representable.');
assert(min(metric)<=requested*(1+64*eps),'ipm_cr:TensorMapSpacing', ...
    'The analytic origin spacing exceeds the requested cap.');
end

function value=log_sinhc(s)
% log(sinh(s)/s), including its removable origin singularity.
if s==0,value=0;
elseif s<.1
    z=s*s;
    value=log1p(z/6+z^2/120+z^3/5040+z^4/362880+z^5/39916800+z^6/6227020800);
elseif s<40,value=log(sinh(s)/s);
else,value=s-log(2*s)+log1p(-exp(-2*s));end
end
