function result = test_ipm_c_gauge()
%TEST_IPM_C_GAUGE Synthetic wall-transfer and complete-RHS gauge checks.
% Requires the original runtime's +ipm package on the MATLAB path. No PDE,
% Darcy solve, checkpoint edit or reference initialization is performed.

r=linspace(0,2,161)';theta=linspace(0,pi,25);
X=r*cos(theta);Y=r*sin(theta);a=.02;
R=X-a*X.^3+.04*Y.^2;Q=1-3*a*X.^2;
x=linspace(-1.9,1.9,81);n=numel(x);h=x(2)-x(1);
Dx=spdiags([-ones(n,1),zeros(n,1),ones(n,1)],[-1,0,1],n,n)/(2*h);
Dx(1,:)=0;Dx(1,1:3)=[-3,4,-1]/(2*h);
Dx(end,:)=0;Dx(end,end-2:end)=[1,-4,3]/(2*h);
settings=struct('wallIntegralRadius',1,'cOmegaGaugeFloor',1e-14, ...
    'pointStrainConditionFloor',1e-12,'referenceWallIntegral',17.25, ...
    'referenceWallGradientCentroid',.17,'peakPositionRestoringRate',.5, ...
    'pinX',0,'gaugeRateGuess',[.1,-.2]);
ops=struct('x',x,'Dx',Dx,'rescaling',settings);saved=ops;
m=ipm_c.gauge.measure(R,r,ops);
assert(strcmp(m.interpolation,'linear') && m.dataLinearInterpolation, ...
    'polar:testDefaultInterpolation','Gauge interpolation must default to linear.');
signedX=[-flipud(r(2:end));0;r(2:end)];
signedR=[flipud(R(2:end,end));R(1,1);R(2:end,1)];
expectedR=interp1(signedX,signedR,x,'linear');expectedQ=expectedR*Dx';
assert(max(abs(m.q-expectedQ))<1e-13, ...
    'polar:testOriginalDerivative','Measurement must apply the original Dx after transfer.');
assert(isequaln(ops,saved) && m.referenceWallIntegral==17.25 && ...
    m.referenceWallGradientCentroid==.17 && abs(m.I-17.25)>.1, ...
    'polar:testReferences','Saved references must not be reset to the mapped state.');
mp=ipm_c.gauge.measure(R,r,ops,struct('interpolation','pchip'));
assert(~mp.dataLinearInterpolation && isequaln(ops,saved), ...
    'polar:testPchipDiagnostic','PCHIP must be explicit and identified as nonlinear in data.');
mt=ipm_c.gauge.measure(R',r,ops,struct('radialDimension',2));
assert(isequaln(mt.q,m.q),'polar:testOrientation','Both declared radial dimensions must agree.');

Uref=.31;guess=[.1,-.2];
callback=@(cl,delta)-((cl+.23)*X+(delta-.16)).*Q+.07*R;
options=struct('rootTolerance',1e-9,'auditTolerance',1e-8);
[F,rates,info]=ipm_c.gauge.solve(R,r,ops,callback,Uref,guess,options);
rateOptions=struct('rateField',F);mf=ipm_c.gauge.measure(R,r,ops,rateOptions);
target=[0;0;settings.peakPositionRestoringRate*(m.C-settings.referenceWallGradientCentroid)];
assert(max(abs(mf.rate.normalized+target))<1e-8, ...
    'polar:testFullRhs','The actual complete polar RHS must satisfy original I/P/C rates.');
assert(abs(rates.c_r-(rates.delta-Uref))<1e-14 && ...
    isequaln(info.originalGuess,guess) && info.Uref==Uref && ...
    isequaln(ops,saved),'polar:testGaugeCoordinates', ...
    'Physical translation, Uref, original guess and references must be retained.');
% Independent affine 3x3 reference using the same original wall functionals.
b0=callback(0,0);b1=callback(1,0)-b0;b2=callback(0,1)-b0;
v0=observe(b0);matrix=[observe(b1),observe(b2),observe(R)];
expected=matrix\(-v0-target);
actual=[rates.c_l;rates.delta;rates.c_omega];
assert(norm(actual-expected,inf)<1e-6, ...
    'polar:testAffineReference','The 2D elimination must match the affine full 3x3 reference.');

% The origin is one physical point; incompatible ray data must be rejected.
bad=R;bad(1,end)=bad(1,end)+1e-4;rejected=false;
try
    ipm_c.gauge.measure(bad,r,ops);
catch problem
    rejected=strcmp(problem.identifier,'polar:WallOrigin');
end
assert(rejected,'polar:testOriginRejection','Incompatible polar origins must not be averaged away.');
result=struct('passed',true,'I',m.I,'P',m.P,'C',m.C, ...
    'rates',rates,'auditResidual',info.auditResidual, ...
    'affineReferenceDifference',norm(actual-expected,inf), ...
    'originalReferencesPreserved',isequaln(ops,saved));

    function values=observe(field)
        opts=struct('rateField',field);
        item=ipm_c.gauge.measure(R,r,ops,opts);values=item.rate.normalized;
    end
end
