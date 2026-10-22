function rho = initialDensity(ops,physics)
%IPM.FIELD.INITIALDENSITY Physical density data for IPM.

if isa(physics.initialCondition,'function_handle')
    rho = physics.initialCondition(ops.X,ops.Y);
    if ~isequal(size(rho), [ops.ny,ops.nx])
        error('ipm:InitialSize', 'Initial-condition function returned the wrong size.');
    end
    return;
end

switch lower(string(physics.initialCondition))
    case "unstable_front"
        envelope = exp(-(ops.X/physics.frontWidth).^4);
        front = physics.frontHeight + physics.frontAmplitude*envelope.* ...
            cos(pi*ops.X/physics.frontWidth);
        % Larger rho lies above smaller rho: the Rayleigh-Taylor unstable
        % orientation for gravity in the negative y direction.
        rho = 0.5*(1+tanh((ops.Y-front)/physics.frontThickness));
    case "smooth_blob"
        rho = exp(-((ops.X/1.1).^2+((ops.Y-1.4)/0.7).^2));
    case "degenerate"
        k = physics.degeneratePower;
        if strcmp(ops.symmetryMode,'double_odd_omega')
            % The upper trace is even in x_1. Its extension across x_2=0
            % is odd and has the explicitly accepted jump at that axis.
            absoluteX = abs(ops.X);
            rho = absoluteX.^k ./ ...
                (1+absoluteX.^k+ops.Y.^k);
        else
            rho = zeros(ops.ny,ops.nx);
            positiveX = ops.X > 0;
            rho(positiveX) = ops.X(positiveX).^k ./ ...
                (1+ops.X(positiveX).^k+ops.Y(positiveX).^k);
        end
    case "degenerate_primitive"
        k = physics.degeneratePower;
        % Analytic primitive of x^(k-1)/(1+x^k+y^k), anchored by
        % rho(0,y)=0.
        if strcmp(ops.symmetryMode,'double_odd_omega')
            absoluteX = abs(ops.X);
            rho = log1p(absoluteX.^k./(1+ops.Y.^k))/k;
        else
            rho = zeros(ops.ny,ops.nx);
            positiveX = ops.X > 0;
            rho(positiveX) = log1p(ops.X(positiveX).^k ./ ...
                (1+ops.Y(positiveX).^k))/k;
        end
    case "shifted_degenerate_primitive"
        % Selected 256-cell experiment. The single user-facing power k
        % lives in main.m. With a=(k-1)^(1/k) and z=(X+a)_+,
        %   R_0=log(1+z^k/(1+Y^k))/k,
        %   R_{0,X}(X,0)=z^(k-1)/(1+z^k).
        % The continuous wall-gradient maximum is at z=a, hence X=0.
        % This primitive is anchored to zero for X<=-a. Adding an arbitrary
        % constant is mathematically allowed but degrades floating accuracy.
        k = physics.degeneratePower;
        a = (k-1)^(1/k);
        z = max(ops.X+a,0);
        rho = log1p(z.^k./(1+ops.Y.^k))/k;
    case "shifted_degenerate_primitive_cutoff"
        % Preserve the shifted primitive exactly through z=B, then
        % saturate its coordinate.  The horizontal Poisson source R_X is
        % compactly supported even though R approaches a positive plateau.
        k = physics.degeneratePower;
        a = (k-1)^(1/k);
        z = max(ops.x(:)'+a,0);
        h = cutoff_coordinate(z,physics.shiftedCutoffB);
        rho = log1p(repmat(h.^k,ops.ny,1)./(1+ops.Y.^k))/k;
    case "localized_gaussian_primitive"
        % Prescribe a rapidly decaying positive R_X, then integrate in X.
        % R(-infinity,Y)=0. This fixes the Y-dependent primitive constant;
        % an odd-erf primitive would describe a different physical IVP.
        A=physics.gaussianAmplitude;a=physics.gaussianCenterX;
        sx=physics.gaussianWidthX;sy=physics.gaussianWidthY;
        w=physics.gaussianSkewWeight;z=(ops.X-a)/sx;
        % erfc avoids cancellation in 1+erf(z) on the negative tail.
        primitive=(1-w)*erfc(-z)+w*erfc(-(z-.8)/1.25);
        rho=A*sx*sqrt(pi)/2*primitive.*exp(-(ops.Y/sy).^2);
    case "ccf_heavy_tail"
        % A resolved approximation to the CCF boundary profile
        %   d_x rho = 1_{x>1}/sqrt(x-1).
        % For z=(|x|-a)_+, the primitive below has
        %   d_|x| rho = A*z/(z^2+epsilon^2)^(3/4),
        % which is finite initially and converges pointwise to the target
        % tail away from its onset as epsilon tends to zero.  The algebraic
        % wall-normal envelope keeps the source integrable in y without
        % changing its wall trace.
        if strcmp(ops.symmetryMode,'double_odd_omega')
            tailCoordinate = abs(ops.X);
        else
            tailCoordinate = ops.X;
        end
        z = max(tailCoordinate-physics.ccfTailOnset,0);
        unboundedPrimitive = 2*physics.ccfTailAmplitude * ...
            ((z.^2+physics.ccfTailRegularization^2).^(1/4) - ...
            sqrt(physics.ccfTailRegularization));
        if isfinite(physics.ccfTailCutoffScale)
            saturation = 2*physics.ccfTailAmplitude* ...
                sqrt(physics.ccfTailCutoffScale);
            radialPrimitive = saturation*tanh(unboundedPrimitive/saturation);
        else
            radialPrimitive = unboundedPrimitive;
        end
        verticalEnvelope = 1 ./ ...
            (1+(ops.Y/physics.ccfTailVerticalScale).^ ...
            physics.ccfTailVerticalPower);
        rho = radialPrimitive.*verticalEnvelope;
    otherwise
        error('ipm:UnknownInitialCondition', ...
            'Unknown initialCondition "%s".',physics.initialCondition);
end
end

function h = cutoff_coordinate(z,B)
% h'=1-I_t(8,8) on B<z<2B, with t=(z-B)/B.  The integral is exact:
% h=B+B*[t*(1-I_t(8,8))+(1/2)*I_t(9,8)].  The regularized beta
% polynomials give a C^8 coordinate without numerical quadrature;
% betainc(1-t,8,8) avoids subtracting two nearly equal numbers.
h = min(z,1.5*B);
transition = z>B & z<2*B;
t = (z(transition)-B)/B;
tailDerivative = betainc(1-t,8,8);
h(transition) = B*(1+t.*tailDerivative+ ...
    0.5*betainc(t,9,8));
end
