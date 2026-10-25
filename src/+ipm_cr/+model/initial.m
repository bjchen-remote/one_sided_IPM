function rho=initial(ops,physics)
%INITIAL Original Gaussian datum; historical cases belong to legacy/solvers.
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
end
