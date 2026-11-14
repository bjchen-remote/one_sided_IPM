classdef observe
    %OBSERVE Original wall traces, I/P/C observables and their time derivatives.
    methods (Static)
        function [measurement,rateEvaluator]=wall(field,r,wallOps,options)
            %WALLMEASURE Original I/P/C with a constant-shifted natural-spline peak.
            % Subtracting max(q) changes neither the mathematical natural spline nor
            % its maximizing location/functional rows. It prevents a P-scaled tie band
            % from selecting a nearby node instead of the stationary maximum. In the
            % shifted data every node is nonpositive; the local overshoot and its tie
            % tolerance then have the same small scale near a knot crossing.
            % I, C, all original references, wall sampling, and rate sampling remain
            % those of ipm_cr.gauge.measure. Its bound rate closure uses the same single
            % centered peak computation. Independent endpoint audits call wall afresh.
            if nargin<4,options=struct();end
            assert(isstruct(options)&&isscalar(options),'polar:WallOptions', ...
                'Measurement options must be a scalar structure.');
            options.centerPeak=true;
            [measurement,rateEvaluator]=ipm_cr.gauge.observe.measure(field,r,wallOps,options);
        end

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
            % rateFieldOffset/wallRateOffset.
            % Offsets do not redefine observables or references. Their rate values must
            % be supplied consistently by the caller for every trial gauge/RHS.
            %
            % The default interpolation is LINEAR IN THE DATA. This is necessary for
            % exact c_omega elimination: M(B+c_omega*R)=M(B)+c_omega*M(R). 'spline'
            % is also data-linear; 'pchip' is an explicit diagnostic option and is not.
            % The nonlinear pchip variant is retained for explicit diagnostic use.
            % options.centerPeak=true computes the original global spline maximum from
            % q-max(q) once, then restores its amplitude; I/C and fixed references stay
            % unchanged. The original I/C observables use finite-axis piecewise-linear
            % and trapezoidal quadrature, not exact continuum wall integrals.


            if nargin<4,options=struct();end
            assert(isstruct(options) && isscalar(options), ...
                'polar:WallOptions','Measurement options must be a scalar structure.');
            method='linear';radialDimension=1;centerPeak=false;
            if isfield(options,'interpolation'),method=options.interpolation;end
            if isfield(options,'radialDimension'),radialDimension=options.radialDimension;end
            if isfield(options,'centerPeak'),centerPeak=options.centerPeak;end
            validateattributes(centerPeak,{'logical'},{'scalar'});
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
            q=ipm_cr.gauge.observe.wallDiff(wallValues,wallOps.Dx);
            integral=ipm_cr.gauge.observe.integral(q,wallOps,wallOps.rescaling);
            if centerPeak
                anchor=max(q);peak=ipm_cr.gauge.peak(q-anchor,wallOps);
                peak.evaluationAnchor=anchor;
                peak.centeredMaximumCandidateValue=peak.maximumCandidateValue;
                peak.value=peak.value+anchor;
                peak.gridValue=peak.gridValue+anchor;
                peak.maximumCandidateValue=peak.maximumCandidateValue+anchor;
                peak.tieToleranceUsesCenteredValues=true;
            else
                peak=ipm_cr.gauge.peak(q,wallOps);
            end
            centroid=ipm_cr.gauge.observe.centroid(q,wallOps);
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
            if centerPeak
                measurement.peakEvaluation='constant_shifted_original_global_natural_cubic_spline';
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
                    [wallRate,~,~,rateOriginMismatch,rateOffsetUsed]= ...
                        sample_with_offset(rhs,'rateFieldOffset','wallRateOffset',rateOptions);
                    qDot=ipm_cr.gauge.observe.wallDiff(wallRate,wallOps.Dx);
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
                    assert(isnumeric(offset) && isreal(offset) && isequal(size(offset),size(polarValues)), ...
                        'polar:WallOffsetShape','The polar offset must have the complete field shape.');
                    assert(isnumeric(wallOffset) && isreal(wallOffset) && isvector(wallOffset) && ...
                        numel(wallOffset)==n,'polar:WallOffsetShape', ...
                        'The exact wall offset must match the original x grid.');
                    [~,axis,values,mismatch]=sample_wall(polarValues);
                    sampled=sample_wall(polarValues-offset)+reshape(wallOffset,1,[]);
                end

                function [sampled,axis,values,mismatch]=sample_wall(polarValues)
                    validateattributes(polarValues,{'numeric'},{'real','2d','nonempty'});
                    if radialDimension==2,polarValues=polarValues';end
                    assert(size(polarValues,1)==numel(r) && size(polarValues,2)>=2, ...
                        'polar:WallFieldSize','Field must contain both theta=0/pi wall rays.');
                    positive=polarValues(:,1);negative=polarValues(:,end);
                    mismatch=abs(positive(1)-negative(1));
                    origin=positive(1)+(negative(1)-positive(1))/2;
                    axis=[-flipud(r(2:end));0;r(2:end)];
                    values=[flipud(negative(2:end));origin;positive(2:end)];
                    sampled=interp1(axis,values,x,method);
                    sampled=reshape(sampled,1,[]);
                end
        end

        function integral = integral(wallGradient,ops,settings)
            %IPM.EVOLVE.WALLGRADIENTINTEGRAL Integrate R_X on the wall over [-L,L].
            %   L=wallIntegralRadius is fixed in rescaled X coordinates. The weights
            %   integrate the piecewise-linear wall-gradient interpolant exactly, including
            %   partial end cells when -L or L is not a mesh node. Weights are rebuilt
            %   after remeshing; the coordinate interval and initial reference are fixed.

            x = reshape(ops.x,1,[]);
            radius = settings.wallIntegralRadius;

            left = x(1:end-1);
            right = x(2:end);
            lower = max(left,-radius);
            upper = min(right,radius);
            active = upper > lower;
            indices = find(active);
            cellWidth = right(active)-left(active);
            clippedWidth = upper(active)-lower(active);
            averageCoordinate = ((lower(active)-left(active)) + ...
                (upper(active)-left(active)))./(2*cellWidth);
            weights = zeros(size(x));
            weights(indices) = clippedWidth.*(1-averageCoordinate);
            weights(indices+1) = weights(indices+1) + ...
                clippedWidth.*averageCoordinate;
            support = weights > 0;
            supportPoints = nnz(support);

            wallGradient = reshape(wallGradient,1,[]);
            value = sum(weights.*wallGradient);
            absoluteValue = sum(weights.*abs(wallGradient));
            condition = value/max(absoluteValue,realmin);
            integral = struct('weights',weights,'value',value, ...
                'radius',radius,'center',0,'support',support, ...
                'supportPoints',supportPoints,'condition',condition, ...
                'absoluteValue',absoluteValue,'meanValue',value/(2*radius));
        end

        function centroid=centroid(values,ops)
            %WALLGRADIENTCENTROID Positive-quadrature centroid of wall R_X squared.
            % On fixed [-L,L], q is sampled by piecewise-linear interpolation at the
            % clipped endpoints and unchanged interior nodes. Positive trapezoidal
            % weights w give E=sum(w*q^2), C=sum(w*x*q^2)/E, and the exact derivative of
            % this discrete observable Cdot=sum(2*w*(x-C)*q*qdot)/E. Endpoint sampling
            % is included in the returned global-node functional. No peak curvature,
            % local branch, time-step, or basin condition is used.
            x=ops.x(:)';q=values(:)';L=1;
            if isfield(ops,'rescaling') && isfield(ops.rescaling,'wallIntegralRadius')
                L=ops.rescaling.wallIntegralRadius;
            end
            inside=find(x>-L & x<L);z=[-L,x(inside),L];n=numel(z);
            B=sparse(n,numel(x));B(1,:)=sample_row(x,-L);B(end,:)=sample_row(x,L);
            for j=1:numel(inside),B(j+1,inside(j))=1;end
            h=diff(z);w=[h(1),h(1:end-1)+h(2:end),h(end)]/2;
            sampled=(B*q')';qScale=max(abs(sampled));
            scaled=sampled/qScale;
            normalizedEnergy=sum(w.*scaled.^2);
            value=sum(w.*z.*scaled.^2)/normalizedEnergy;
            weights=(2*(w.*(z-value).*scaled)/(qScale*normalizedEnergy))*B;
            energyWeights=(2*(w.*sampled))*B;
            indices=find(weights~=0 | energyWeights~=0);
            reference=value;referencePersisted=false;
            if isfield(ops,'rescaling') && isfield(ops.rescaling,'referenceWallGradientCentroid')
                reference=ops.rescaling.referenceWallGradientCentroid;
                referencePersisted=true;
            end
            centroid=struct('definition','fixed_window_wall_R_X_squared_centroid', ...
                'value',value,'energy',sum(sampled.*(w.*sampled)), ...
                'normalizedEnergy',normalizedEnergy,'normalizationScale',qScale, ...
                'degenerateEnergy',qScale==0,'radius',L,'window',[-L,L], ...
                'indices',indices,'weights',weights(indices), ...
                'energyWeights',energyWeights(indices), ...
                'samplePoints',z,'sampleValues',sampled, ...
                'quadratureWeights',w,'sampleMap',B, ...
                'referenceValue',reference,'referencePersisted',referencePersisted, ...
                'quadratureOrder',2);
        end

        function q=wallDiff(values,Dx)
            %WALLDERIVATIVE Apply the original wall derivative to local value differences.
            % Each stencil annihilates constants before multiplication; large platforms
            % cannot produce a derivative from the floating point row sum of Dx.
            persistent previous rows columns weights
            if isempty(previous)||~isequal(previous,Dx)
                [rows,columns,weights]=find(Dx);previous=Dx;
            end
            values=values(:);
            assert(numel(values)==size(Dx,2), ...
                'polar:WallDerivative','The wall data must match the original derivative axis.');
            q=accumarray(rows,weights.*(values(columns)-values(rows)), ...
                [size(Dx,1),1],@sum,0).';
        end
    end
end

function row=sample_row(x,query)
row=sparse(1,numel(x));index=find(x==query,1);
if ~isempty(index),row(index)=1;return;end
left=find(x<query,1,'last');
t=(query-x(left))/(x(left+1)-x(left));
row(left:left+1)=[1-t,t];
end
