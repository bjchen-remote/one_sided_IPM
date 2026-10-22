function [measurement,rateEvaluator] = measure(field,r,wallOps,options)
%POLAR_WALL_MEASUREMENT Evaluate the original Cartesian wall I/P/C gauges.
% field is [numel(r),ntheta], with columns 1/end at theta=0/pi. Set
% options.radialDimension=2 for [ntheta,numel(r)]. The two rays are joined
% into an increasing signed X axis, with the origin represented once.
% wallOps.x, wallOps.Dx and wallOps.rescaling are the ORIGINAL measurement
% operators and persisted references; this routine never initializes them.
% options.rateField, when present, is a full polar RHS of the same size.
% Optional second output rateEvaluator(rhs,rateOptions) binds this exact
% density and geometry. It reuses the fixed I/P/C functionals, sampling only
% rhs and its paired rate offsets. No cross-state persistent cache is used.
% With paired options.fieldOffset/wallOffset, sample field-fieldOffset and
% add wallOffset on ORIGINAL x. Rate offsets are independently paired as
% rateFieldOffset/wallRateOffset (rateWallOffset is an accepted alias).
% Offsets do not redefine observables or references. Their rate values must
% be supplied consistently by the caller for every trial gauge/RHS.
%
% The default interpolation is LINEAR IN THE DATA. This is necessary for
% exact c_omega elimination: M(B+c_omega*R)=M(B)+c_omega*M(R). 'spline'
% is also data-linear; 'pchip' is an explicit diagnostic option and is not.
% A caller using pchip must independently audit the assembled complete RHS.
% No claim of equivalence to the old Cartesian WENO operator is made.

if nargin<4,options=struct();end
assert(isstruct(options) && isscalar(options), ...
    'polar:WallOptions','Measurement options must be a scalar structure.');
method='linear';radialDimension=1;
if isfield(options,'interpolation'),method=options.interpolation;end
if isfield(options,'radialDimension'),radialDimension=options.radialDimension;end
assert(ischar(method) || (isstring(method) && isscalar(method)), ...
    'polar:WallInterpolation','Specify linear, spline or pchip interpolation.');
method=char(method);
assert(ismember(method,{'linear','spline','pchip'}), ...
    'polar:WallInterpolation','Only linear, spline and pchip are supported.');
assert(isscalar(radialDimension) && ismember(radialDimension,[1,2]), ...
    'polar:WallOrientation','radialDimension must be 1 or 2.');
validateattributes(r,{'numeric'},{'real','finite','vector','increasing'});
r=r(:);
assert(numel(r)>=3 && r(1)==0 && r(end)>0, ...
    'polar:WallRadius','The radial grid must start at exactly zero.');
assert(isstruct(wallOps) && isscalar(wallOps) && ...
    all(isfield(wallOps,{'x','Dx','rescaling'})), ...
    'polar:WallOriginalOps','Inject the original Cartesian x, Dx and rescaling.');
x=wallOps.x(:)';n=numel(x);
assert(n>=3 && all(isfinite(x)) && all(diff(x)>0) && ...
    isequal(size(wallOps.Dx),[n,n]), ...
    'polar:WallOriginalOps','The original wall x/Dx must have compatible sizes.');
required={'wallIntegralRadius','referenceWallIntegral', ...
    'referenceWallGradientCentroid','peakPositionRestoringRate'};
assert(all(isfield(wallOps.rescaling,required)), ...
    'polar:WallReferences','Original I/C references and fixed gauge settings are required.');
assert(x(1)>=-r(end) && x(end)<=r(end), ...
    'polar:WallCoverage','Polar wall rays must cover the complete original wall domain.');

[wallValues,signedX,signedValues,originMismatch,fieldOffsetUsed]= ...
    sample_with_offset(field,'fieldOffset','wallOffset',options);
q=ipm_c.gauge.wallDerivative(wallValues,wallOps.Dx);
integral=ipm.evolve.wallGradientIntegral(q,wallOps,wallOps.rescaling);
peak=ipm.evolve.globalSplineMaximum(q,wallOps);
centroid=ipm.evolve.wallGradientCentroid(q,wallOps);
assert(isfinite(peak.value) && peak.value>0 && ~centroid.degenerateEnergy, ...
    'polar:WallDegenerate','The original I/P/C gauges require positive P and nonzero energy.');
measurement=struct('definition','original_Cartesian_wall_IPC_after_polar_sampling', ...
    'x',x,'wallValues',wallValues,'q',q, ...
    'I',integral.value,'P',peak.value,'C',centroid.value, ...
    'integral',integral,'peak',peak,'centroid',centroid, ...
    'interpolation',method,'dataLinearInterpolation',~strcmp(method,'pchip'), ...
    'radialDimension',radialDimension,'signedWallX',signedX, ...
    'signedWallValues',signedValues,'originMismatch',originMismatch, ...
    'referenceWallIntegral',wallOps.rescaling.referenceWallIntegral, ...
    'referenceWallGradientCentroid',wallOps.rescaling.referenceWallGradientCentroid, ...
    'relativeIntegralDrift',(integral.value-wallOps.rescaling.referenceWallIntegral)/ ...
        max(abs(wallOps.rescaling.referenceWallIntegral),realmin), ...
    'centroidDisplacement',centroid.value-wallOps.rescaling.referenceWallGradientCentroid, ...
    'fieldOffsetUsed',fieldOffsetUsed,'rate',[]);
if isfield(options,'rateWallOffset')
    if isfield(options,'wallRateOffset')
        assert(isequaln(options.rateWallOffset,options.wallRateOffset), ...
            'polar:WallRateOffsetAlias','wallRateOffset and rateWallOffset must agree.');
    else
        options.wallRateOffset=options.rateWallOffset;
    end
end
assert((~isfield(options,'rateFieldOffset') && ~isfield(options,'wallRateOffset')) || ...
    (isfield(options,'rateField') && ~isempty(options.rateField)), ...
    'polar:WallRateOffset','Rate offsets require an actual complete rateField.');
if isfield(options,'rateField') && ~isempty(options.rateField)
    measurement.rate=measure_rate(options.rateField,options);
end
rateEvaluator=[];
if nargout>1,rateEvaluator=@prepared_rate;end

    function value=prepared_rate(rhs,rateOptions)
        % The captured measurement is immutable for this invocation of
        % measure. The caller must create a new closure for a new density.
        if nargin<2,rateOptions=struct();end
        value=measurement;value.rate=measure_rate(rhs,rateOptions);
    end
    function rate=measure_rate(rhs,rateOptions)
        assert(isstruct(rateOptions)&&isscalar(rateOptions), ...
            'polar:WallOptions','Rate options must be a scalar structure.');
        assert(isequal(size(rhs),size(field)), ...
            'polar:WallRateSize','rateField must have the same polar shape as field.');
        if isfield(rateOptions,'rateWallOffset')
            if isfield(rateOptions,'wallRateOffset')
                assert(isequaln(rateOptions.rateWallOffset,rateOptions.wallRateOffset), ...
                    'polar:WallRateOffsetAlias','wallRateOffset and rateWallOffset must agree.');
            else
                rateOptions.wallRateOffset=rateOptions.rateWallOffset;
            end
        end
        [wallRate,~,~,rateOriginMismatch,rateOffsetUsed]= ...
            sample_with_offset(rhs,'rateFieldOffset','wallRateOffset',rateOptions);
        qDot=ipm_c.gauge.wallDerivative(wallRate,wallOps.Dx);
        iDot=integral.weights*qDot';
        pDot=peak.weights*qDot(peak.indices)';
        cDot=centroid.weights*qDot(centroid.indices)';
        rate=struct('wallValues',wallRate,'qDot',qDot, ...
            'I_dot',iDot,'P_dot',pDot,'C_dot',cDot, ...
            'normalized',[iDot/integral.value;pDot/peak.value;cDot], ...
            'originMismatch',rateOriginMismatch,'fieldOffsetUsed',rateOffsetUsed);
    end

    function [sampled,axis,values,mismatch,used]=sample_with_offset(polarValues,fieldKey,wallKey,sampleOptions)
        used=isfield(sampleOptions,fieldKey) || isfield(sampleOptions,wallKey);
        if ~used
            [sampled,axis,values,mismatch]=sample_wall(polarValues);
            return;
        end
        assert(isfield(sampleOptions,fieldKey) && isfield(sampleOptions,wallKey), ...
            'polar:WallOffsetPair','Polar and original-wall offsets must be supplied together.');
        offset=sampleOptions.(fieldKey);wallOffset=sampleOptions.(wallKey);
        assert(isnumeric(offset) && isreal(offset) && isequal(size(offset),size(polarValues)) && ...
            all(isfinite(offset),'all'),'polar:WallOffsetShape', ...
            'The polar offset must be finite and have the complete field shape.');
        assert(isnumeric(wallOffset) && isreal(wallOffset) && isvector(wallOffset) && ...
            numel(wallOffset)==n && all(isfinite(wallOffset)), ...
            'polar:WallOffsetShape','The exact wall offset must match the original x grid.');
        % Check original and offset origins separately; subtraction must not
        % hide two incompatible multi-valued physical origin rows.
        [~,axis,values,mismatch]=sample_wall(polarValues);
        sample_wall(offset);
        sampled=sample_wall(polarValues-offset)+reshape(wallOffset,1,[]);
        assert(all(isfinite(sampled)),'polar:WallOffsetFinite','Offset wall sampling must remain finite.');
    end

    function [sampled,axis,values,mismatch]=sample_wall(polarValues)
        validateattributes(polarValues,{'numeric'},{'real','finite','2d','nonempty'});
        if radialDimension==2,polarValues=polarValues';end
        assert(size(polarValues,1)==numel(r) && size(polarValues,2)>=2, ...
            'polar:WallFieldSize','Field must contain both theta=0/pi wall rays.');
        positive=polarValues(:,1);negative=polarValues(:,end);
        mismatch=abs(positive(1)-negative(1));
        tolerance=256*eps(max([1;abs(positive(1));abs(negative(1))]));
        assert(mismatch<=tolerance,'polar:WallOrigin', ...
            'The two wall rays disagree at the common polar origin.');
        origin=positive(1)+(negative(1)-positive(1))/2;
        axis=[-flipud(r(2:end));0;r(2:end)];
        values=[flipud(negative(2:end));origin;positive(2:end)];
        sampled=interp1(axis,values,x,method);
        sampled=reshape(sampled,1,[]);
        assert(all(isfinite(sampled)), ...
            'polar:WallSampling','Original-wall interpolation must remain finite.');
    end
end
