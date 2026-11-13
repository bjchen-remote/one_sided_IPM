function rho=initial(ops,physics)
%INITIAL Registered Gaussian datum with an optional smooth radial cutoff.
% The immutable wall authority admits this datum only. Keep its original
% arithmetic: erfc avoids cancellation in 1+erf on the negative tail.
assert(isfield(physics,'initialCondition')&& ...
    strcmp(physics.initialCondition,'localized_gaussian_primitive'), ...
    'ipm_cr:InitialCondition','The maintained solver uses the original Gaussian datum.');
A=physics.gaussianAmplitude;a=physics.gaussianCenterX;
sx=physics.gaussianWidthX;sy=physics.gaussianWidthY;
w=physics.gaussianSkewWeight;z=(ops.X-a)/sx;
primitive=(1-w)*erfc(-z)+w*erfc(-(z-.8)/1.25);
rho=A*sx*sqrt(pi)/2*primitive.*exp(-(ops.Y/sy).^2);
if isfield(physics,'gaussianCutoff')
    radii=physics.gaussianCutoff;
    validateattributes(radii,{'numeric'},{'real','finite','positive','increasing','vector','numel',2});
    rho=rho.*ipm_cr.mesh.compactCutoff(hypot(ops.X,ops.Y),radii(1),radii(2));
elseif isfield(physics,'gaussianSourceCutoff')
    radii=physics.gaussianSourceCutoff;
    validateattributes(radii,{'numeric'},{'real','finite','positive','increasing','vector','numel',2});
    outside=abs(ops.X)>radii(1)|abs(ops.Y)>radii(1);
    if any(outside,'all')
        model=struct('amplitude',A,'centerX',a,'widthX',sx,'widthY',sy,'skewWeight',w,'cutoff',radii);
        value=ipm_cr.model.compactGaussian(model,ops.X(outside),ops.Y(outside),false);
        rho(outside)=value.rho;
    end
end
end
