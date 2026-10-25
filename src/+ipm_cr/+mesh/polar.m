function grid = polar(options)
%POLAR_GRID Independent stretched polar tensor geometry on a half disk.
% Arrays have size [nr,ntheta].  The repeated r=0 row denotes ONE point.
% options: nr=129, ntheta=129, rmax=1e10, firstSpacing=.02,
% domain='semicircle', radialMap='sinh', stencilWidth=7.
% radialMap='uniform_core' adds a nearly uniform core.
% The default r(s)=R*sinh(alpha*s)/sinh(alpha) has a smooth ODD continuation
% across the signed physical wall.  Explicit radialMap='exp_gaps' retains
% the old increasing-gap geometry for historical diagnostic comparisons.
% This constructor returns axes, coordinates and the radial map metric;
% it does not build origin jets or transport derivative operators.
% No rectangle is replaced by a half disk silently.
% Requires the original ipm_cr.mesh.fd on the MATLAB search path.
if nargin<1,options=struct();end
validateattributes(options,{'struct'},{'scalar'},mfilename,'options');
domain=char(option(options,'domain','semicircle'));
if strcmp(domain,'rectangle')
    error('polar:RectangleDomainNotImplemented', ...
        'A rectangular polar mapping is not implemented; choose semicircle explicitly.');
end
assert(strcmp(domain,'semicircle'),'polar:Domain', ...
    'The only implemented polar domain is semicircle.');
nr=option(options,'nr',129);ntheta=option(options,'ntheta',129);
rmax=option(options,'rmax',1e10);firstSpacing=option(options,'firstSpacing',.02);
radialMap=char(option(options,'radialMap','sinh'));
assert(ismember(radialMap,{'sinh','exp_gaps','uniform_core'}),'polar:RadialMap', ...
    'radialMap must be sinh, uniform_core, or the explicit historical exp_gaps option.');
width=option(options,'stencilWidth',7);
assert(~any(isfield(options,{'originRings','originFitDegree'})), ...
    'ipm_cr:RetiredOption','originRings/originFitDegree no longer affect this origin-free solver.');
validateattributes(nr,{'numeric'},{'scalar','real','finite','integer','>=',7});
validateattributes(ntheta,{'numeric'},{'scalar','real','finite','integer','>=',7});
validateattributes(rmax,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(firstSpacing,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(width,{'numeric'},{'scalar','real','finite','integer','>=',3});
assert(width<=min(nr,ntheta),'polar:StencilWidth', ...
    'The derivative stencil must fit both tensor axes.');
intervals=nr-1;
assert(log(firstSpacing)+log(intervals)<=log(rmax)+32*eps(max(1,abs(log(rmax)))), ...
    'polar:RadialSpacing','firstSpacing must not exceed the average radial spacing.');
if strcmp(radialMap,'uniform_core')
    [r,uniformCoreInfo]=ipm_cr.mesh.radialAxis(options);
    firstSpacing=uniformCoreInfo.firstSpacing;
    mapParameter=uniformCoreInfo.alpha;
else
    logTarget=log(rmax)-log(firstSpacing);
    if abs(logTarget-log(intervals))<=32*eps(max(1,abs(logTarget)))
        mapParameter=0;
    else
        lower=0;upper=1;
        while radial_ratio_log(intervals,upper,radialMap)<logTarget,upper=2*upper;end
        for iteration=1:100
            middle=lower+(upper-lower)/2;
            if radial_ratio_log(intervals,middle,radialMap)<logTarget
                lower=middle;
            else
                upper=middle;
            end
        end
        mapParameter=lower+(upper-lower)/2;
    end
    r=zeros(nr,1);
    for k=1:intervals
        if strcmp(radialMap,'exp_gaps')
            r(k+1)=exp(log(firstSpacing)+geometric_logsum(k,mapParameter));
        elseif mapParameter==0
            r(k+1)=rmax*k/intervals;
        else
            r(k+1)=exp(log(rmax)+log_sinh(mapParameter*k/intervals)-log_sinh(mapParameter));
        end
    end
end
r(2)=firstSpacing;r(end)=rmax;
assert(all(isfinite(r))&&all(diff(r)>0),'polar:RadialRepresentability', ...
    'The requested radial axis must be finite and strictly increasing.');
theta=linspace(0,pi,ntheta);
cosTheta=cos(theta);sinTheta=sin(theta);
cosTheta([1,end])=[1,-1];sinTheta([1,end])=0;
if mod(ntheta,2)==1
    cosTheta((ntheta+1)/2)=0;sinTheta((ntheta+1)/2)=1;
end
X=r*cosTheta;Y=r*sinTheta;
logical=linspace(0,1,nr)';radialMetric=ipm_cr.mesh.fd(logical,1,7)*r;
assert(all(radialMetric>0),'polar:GridMetric', ...
    'The radial mapped WENO metric must remain positive.');
grid=struct('kind','independent_polar_half_disk_tensor_v1','domain',domain, ...
    'nr',nr,'ntheta',ntheta,'rmax',rmax,'firstSpacing',firstSpacing, ...
    'radialMap',radialMap,'radialMapParameter',mapParameter, ...
    'signedWallMapSmoothAtOrigin',ismember(radialMap,{'sinh','uniform_core'}), ...
    'r',r,'theta',theta, ...
    'X',X,'Y',Y,'cosTheta',cosTheta,'sinTheta',sinTheta, ...
    'radialMetric',radialMetric,'stencilWidth',width, ...
    'originIsSinglePhysicalPoint',true,'originalRectanglePreserved',false);
if strcmp(radialMap,'exp_gaps'),grid.radialGapRatio=exp(mapParameter);end
if strcmp(radialMap,'uniform_core'),grid.uniformCoreInfo=uniformCoreInfo;end
end

function value=option(options,name,fallback)
if isfield(options,name),value=options.(name);else,value=fallback;end
end

function value=geometric_logsum(number,logRatio)
% log(1+q+...+q^(number-1)), evaluated without a large intermediate sum.
if logRatio==0
    value=log(number);
else
    value=log_expm1(number*logRatio)-log_expm1(logRatio);
end
end

function value=radial_ratio_log(intervals,parameter,kind)
if strcmp(kind,'exp_gaps')
    value=geometric_logsum(intervals,parameter);
elseif parameter==0
    value=log(intervals);
else
    value=log_sinh(parameter)-log_sinh(parameter/intervals);
end
end

function value=log_sinh(argument)
% log(sinh(argument)), without overflowing sinh for a very large box.
value=argument+log(-expm1(-2*argument))-log(2);
end

function value=log_expm1(argument)
value=argument+log(-expm1(-argument));
end
