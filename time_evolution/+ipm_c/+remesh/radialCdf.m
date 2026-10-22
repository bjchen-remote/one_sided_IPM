function [r,info]=radialCdf(options)
%RADIALCDF Positive analytic monitor on the half-line, independent of N.
% m(r)=1/sqrt(a^2+r^2)+sum_k w_k/sigma_k*(phi((r-mu_k)/sigma_k)
%                                           +phi((r+mu_k)/sigma_k)).
% Each mirrored pair has nominal mass w_k on [0,infinity). The monitor is
% positive, even and analytic, so its inverse CDF has a smooth odd signed
% continuation across the physical wall origin. No PDE field is modified.
% Self-contained radial geometry proposal; no experiments runtime dependency.
% options: nr=129, rmax=1e10, coreScale=.375, featureMass=6,
% featureWeights=[.2 .3 .5], featureBounds=three [lo hi] radial projections,
% maximumAdjacentRatio=6. Missing bounds select an explicitly manufactured
% isotropic exp(-r^2) example at levels [.1 .5 .9], not the physical IVP.
if nargin<1,options=struct();end
assert(isstruct(options)&&isscalar(options),'ipm_c:CdfOptions','Candidate validation failed.');
nr=option(options,'nr',129);R=option(options,'rmax',1e10);
a=option(options,'coreScale',.375);mass=option(options,'featureMass',6);
weights=option(options,'featureWeights',[.2,.3,.5]);
ratioLimit=option(options,'maximumAdjacentRatio',6);
validateattributes(nr,{'numeric'},{'scalar','real','finite','integer','>=',9});
validateattributes(R,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(a,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(mass,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(weights,{'numeric'},{'vector','numel',3,'real','finite','positive'});
validateattributes(ratioLimit,{'numeric'},{'scalar','real','finite','>',1});
weights=reshape(weights,1,[]);
assert(abs(sum(weights)-1)<=64*eps,'ipm_c:CdfWeights', ...
    'Feature weights must sum to one; their fixed mass cannot depend on N.');
defaultBounds=~isfield(options,'featureBounds');
if defaultBounds
    bounds=[zeros(3,1),sqrt(-log([.1;.5;.9]))];
else
    bounds=options.featureBounds;
end
assert(isnumeric(bounds)&&isreal(bounds)&&isequal(size(bounds),[3,2])&& ...
    all(isfinite(bounds),'all')&&all(bounds(:,1)>=0)&& ...
    all(bounds(:,2)>bounds(:,1))&&all(bounds(:,2)<=R), ...
    'ipm_c:CdfBounds','Supply three positive-width radial feature intervals within [0,rmax].');
mu=mean(bounds,2)';sigma=.5*diff(bounds,1,2)';componentMass=mass*weights;
total=monitor_mass(R);targets=(0:nr-1)'*(total/(nr-1));
assert(isfinite(total)&&total>0,'ipm_c:CdfMass','Candidate validation failed.');
lower=zeros(nr,1);upper=R*ones(nr,1);
% Fixed iterations and identical target fractions preserve common dyadic
% quantiles across nested node counts, without using a tabulated old mesh.
for iteration=1:128
    middle=lower+(upper-lower)/2;
    below=monitor_mass(middle)<targets;
    lower(below)=middle(below);upper(~below)=middle(~below);
end
r=lower+(upper-lower)/2;r([1,end])=[0,R];
assert(all(isfinite(r))&&all(diff(r)>0),'ipm_c:CdfMonotonicity','Candidate validation failed.');
width=diff(r);adjacent=max([width(2:end)./width(1:end-1); ...
    width(1:end-1)./width(2:end)]);
assert(adjacent<=ratioLimit*(1+128*eps),'ipm_c:CdfAdjacentRatio', ...
    'The requested count/monitor exceeds its declared adjacent-cell ratio.');
observedMass=monitor_mass(r);quantileResidual=max(abs(observedMass-targets))/total;
analyticMetric=total./monitor_density(r);
assert(quantileResidual<=512*eps&&all(isfinite(analyticMetric))&& ...
    all(analyticMetric>0),'ipm_c:CdfQuantile','Candidate validation failed.');
info=struct('kind','positive_mirrored_gaussian_radial_cdf_v1', ...
    'coreScale',a,'featureBounds',bounds,'featureMass',mass, ...
    'featureWeights',weights,'componentNominalHalfLineMasses',componentMass, ...
    'componentMeans',mu,'componentSigma',sigma,'totalMonitorMass',total, ...
    'massStep',total/(nr-1),'cdfMassAtNodes',observedMass, ...
    'analyticMetricAtNodes',analyticMetric,'maximumQuantileResidual',quantileResidual, ...
    'maximumAdjacentRatio',adjacent,'maximumAdjacentRatioLimit',ratioLimit, ...
    'firstSpacing',r(2),'minimumSpacing',min(width), ...
    'monitorIndependentOfNodeCount',true,'signedMonitorEven',true, ...
    'signedInverseMapOddAndSmooth',true,'defaultManufacturedBounds',defaultBounds, ...
    'defaultBoundsDefinition','isotropic exp(-r^2), levels [.1 .5 .9]; not the original Gaussian IVP', ...
    'pdeFieldChanged',false,'sourceMaskUsed',false,'spatiallyQualified',false, ...
    'performanceTargetQualified',false);

    function value=monitor_mass(location)
        value=asinh(location/a);
        for k=1:3
            value=value+.5*componentMass(k)*( ...
                erf((location-mu(k))/(sigma(k)*sqrt(2)))+ ...
                erf((location+mu(k))/(sigma(k)*sqrt(2))));
        end
    end
    function value=monitor_density(location)
        value=1./hypot(location,a);
        for k=1:3
            minus=(location-mu(k))/sigma(k);plus=(location+mu(k))/sigma(k);
            value=value+(componentMass(k)/(sigma(k)*sqrt(2*pi)))*( ...
                exp(-.5*minus.^2)+exp(-.5*plus.^2));
        end
    end
end

function value=option(options,name,fallback)
if isfield(options,name),value=options.(name);else,value=fallback;end
end
