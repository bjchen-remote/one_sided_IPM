function [r,drMetric,info]=radial(original,options)
%REDISTRIBUTERADIALAXIS Fixed smooth refinement of a uniform_core radial axis.
% u is the OLD unit-index coordinate, 0<=u<=baseNr-1. New integer j samples
% q(u)=u+c*integral_0^u chi(v)dv, chi=betainc(t,8,8) on the transition.
% This is a one-time geometry construction, not ALE or time adaptation.
% Returned drMetric differentiates r with respect to the NEW unit index.
% Complete FD7 rows inside the exact-copy prefix retain old metric/N.
if nargin<2,options=struct();end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options), ...
    {'enabled','factor','preserveRadius','fullRefinementRadius'})), ...
    'ipm_cr:RadialRedistributionOptions','Unknown radial redistribution option.');
enabled=option(options,'enabled',false);factor=option(options,'factor',2);
preserve=option(options,'preserveRadius',4);full=option(options,'fullRefinementRadius',20);
validateattributes(enabled,{'numeric','logical'},{'scalar','real','finite'});
assert(any(enabled==[0,1]),'ipm_cr:RadialRedistributionOptions','enabled must be a logical flag.');
assert(isstruct(original)&&isscalar(original)&&all(isfield(original, ...
    {'nr','r','radialMetric','radialMap'})), ...
    'ipm_cr:RadialRedistributionBase','Use a complete original ipm_cr.mesh.polar geometry.');
baseR=original.r(:);N=original.nr-1;
assert(numel(baseR)==N+1&&N>=6&&baseR(1)==0&&all(isfinite(baseR))&& ...
    all(diff(baseR)>0)&&numel(original.radialMetric)==N+1, ...
    'ipm_cr:RadialRedistributionBase','The base radial axis must be finite and increasing from zero.');
baseDr=original.radialMetric(:)/N;
assert(all(isfinite(baseDr))&&all(baseDr>0),'ipm_cr:RadialRedistributionMetric');
r=baseR;drMetric=baseDr;
info=struct('enabled',logical(enabled),'baseNr',N+1,'actualNr',N+1, ...
    'monitorVersion','old_index_beta_8_8_integral_v1','geometryFixedInTime',true);
if ~enabled,return;end
validateattributes(factor,{'numeric'},{'scalar','real','finite','>=',1});
validateattributes(preserve,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(full,{'numeric'},{'scalar','real','finite','>',preserve,'<=',baseR(end)});
assert(strcmp(original.radialMap,'uniform_core')&&isfield(original,'uniformCoreInfo'), ...
    'ipm_cr:RadialRedistributionMap','Enabled redistribution requires the original uniform_core analytic map.');
map=original.uniformCoreInfo;
assert(all(isfield(map,{'a','alpha','power','tailCoefficient'}))&& ...
    all(isfinite([map.a,map.alpha,map.power,map.tailCoefficient]))&& ...
    map.a>0&&map.alpha>=0&&map.power==3&&map.tailCoefficient>0, ...
    'ipm_cr:RadialRedistributionMap','Invalid original uniform_core map parameters.');
ua=inverse_radius(preserve);ub=inverse_radius(full);
assert(0<ua&&ua<ub&&ub<=N,'ipm_cr:RadialRedistributionMap');
integralN=N-(ua+ub)/2;
K=round((factor-1)*integralN);Nnew=N+K;
assert(isfinite(K)&&K>=0&&isfinite(Nnew)&&Nnew<flintmax, ...
    'ipm_cr:RadialRedistributionCount','The requested radial interval count is not representable.');
c=K/integralN;protectedCount=floor(ua)+1;
assert(protectedCount>=7,'ipm_cr:RadialRedistributionProtection', ...
    'The protected radial prefix must contain a complete seven-point stencil.');
u=(0:Nnew)';preservedCount=protectedCount;
if K==0
    % factor=1 or a rounded zero insertion is a literal geometry no-op.
    u=(0:N)';preservedCount=N+1;copiedMetricRows=(1:N+1)';
else
    j=u;protected=j<=floor(ua);
    atUpper=ub+c*(ub-ua)/2;
    far=~protected&j>=atUpper;
    u(far)=(j(far)+c*(ua+ub)/2)/(1+c);
    transition=~protected&~far;
    lower=ua+zeros(nnz(transition),1);upper=ub+zeros(nnz(transition),1);
    target=j(transition);
    for iteration=1:80
        middle=lower+(upper-lower)/2;
        below=middle+c*monitor_integral(middle,ua,ub)<target;
        lower(below)=middle(below);upper(~below)=middle(~below);
    end
    u(transition)=lower+(upper-lower)/2;
    u(1:protectedCount)=(0:protectedCount-1)';u(end)=N;
    r=analytic_radius(u);
    % Never reevaluate exact protected integer nodes through exp/log.
    r(1:protectedCount)=baseR(1:protectedCount);r(end)=baseR(end);
    assert(all(isfinite(r))&&all(diff(r)>0),'ipm_cr:RadialRedistributionAxis', ...
        'The redistributed analytic radial axis must remain strictly increasing.');
    drMetric=ipm_cr.mesh.fd((0:Nnew)',1,7)*r;
    rows=(1:Nnew+1)';stencilFirst=min(max(rows-3,1),Nnew+1-6);
    copiedMetricRows=find(stencilFirst+6<=protectedCount);
    drMetric(copiedMetricRows)=baseDr(copiedMetricRows);
end
assert(all(isfinite(drMetric))&&all(drMetric>0), ...
    'ipm_cr:RadialRedistributionMetric','The redistributed mapped radial metric must be finite and positive.');
inverseResidual=max(abs(u+c*monitor_integral(u,ua,ub)-(0:Nnew)'));
assert(inverseResidual<=256*eps(max(1,Nnew)), ...
    'ipm_cr:RadialRedistributionInverse','The smooth monitor inverse did not converge to machine brackets.');
gaps=diff(r);
info=struct('enabled',true,'kind','fixed_old_index_smooth_radial_refinement_v1', ...
    'monitorVersion','old_index_beta_8_8_integral_v1','geometryFixedInTime',true, ...
    'factor',factor,'preserveRadius',preserve,'fullRefinementRadius',full, ...
    'baseNr',N+1,'actualNr',Nnew+1,'baseIntervals',N,'actualIntervals',Nnew, ...
    'addedIntervals',K,'ua',ua,'ub',ub,'monitorAmplitude',c, ...
    'actualTailDensityFactor',1+c,'integratedMonitorAtOuter',integralN, ...
    'protectedPrefixCount',protectedCount,'preservedPrefixCount',preservedCount, ...
    'preservedLastRadius',r(preservedCount), ...
    'copiedMetricPrefixCount',numel(copiedMetricRows),'copiedMetricRows',copiedMetricRows, ...
    'metricDerivativeWidth',7,'oldIndexAtNewNodes',u, ...
    'baseR',baseR,'baseDrMetric',baseDr,'maximumInverseResidual',inverseResidual, ...
    'minimumRadialGap',min(gaps),'maximumRadialGap',max(gaps), ...
    'maximumAdjacentGapRatio',max([gaps(2:end)./gaps(1:end-1);gaps(1:end-1)./gaps(2:end)]), ...
    'analyticMap',struct('a',map.a,'alpha',map.alpha,'power',map.power, ...
        'tailCoefficient',map.tailCoefficient), ...
    'metricDefinition','FD7 on new integer index; original radialMetric/baseIntervals on complete protected stencils');

    function value=analytic_radius(oldIndex)
        s=oldIndex/N;
        if map.alpha==0
            tail=s.^map.power;
        else
            tail=zeros(size(s));positive=s>0;
            tail(positive)=exp(map.power*(log_sinh(map.alpha*s(positive))-log_sinh(map.alpha)));
        end
        value=map.a*s+map.tailCoefficient*tail;
        value(oldIndex==0)=0;value(oldIndex==N)=baseR(end);
    end

    function value=inverse_radius(radius)
        if radius==baseR(end),value=N;return;end
        lower=0;upper=N;
        for iteration=1:80
            middle=lower+(upper-lower)/2;sample=analytic_radius(middle);
            if sample==radius,lower=middle;upper=middle;break;end
            if sample<radius,lower=middle;else,upper=middle;end
        end
        value=lower+(upper-lower)/2;
    end
end

function value=monitor_integral(u,ua,ub)
value=zeros(size(u));upper=u>=ub;value(upper)=u(upper)-(ua+ub)/2;
inside=u>ua&u<ub;
if any(inside)
    t=(u(inside)-ua)/(ub-ua);
    value(inside)=(ub-ua)*(t.*betainc(t,8,8)-.5*betainc(t,9,8));
end
end

function value=log_sinh(argument)
value=argument+log(-expm1(-2*argument))-log(2);
end

function value=option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
