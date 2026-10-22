function [g,cart,rho,reference]=ipm_c_test_fixture(n)
%IPM_C_TEST_FIXTURE Analytic datum and fixed original wall operators.
% No experimental trajectory or runtime path is needed. Reference wall
% observables are initialized once from the direct Cartesian Gaussian datum.
if nargin<1,n=33;end
R=4;alpha=3;h=R*sinh(alpha/(n-1))/sinh(alpha);
g=ipm_c.mesh.build(struct('nr',n,'ntheta',n,'rmax',R,'firstSpacing',h));
p=struct('initialCondition','localized_gaussian_primitive','symmetryMode','half_plane', ...
    'gaussianAmplitude',.7,'gaussianCenterX',0,'gaussianWidthX',.75, ...
    'gaussianWidthY',1,'gaussianSkewWeight',.2);
scale=struct('logC_l',0,'logC_omega',0,'X_shift',0, ...
    'canonicalTime',0,'physicalTime',0);
x=linspace(-3.8,3.8,161);Dx=ipm.mesh.fdMatrix(x,1,7);
wallOps=struct('x',x,'Dx',Dx);
settings=struct('wallIntegralRadius',1,'cOmegaGaugeFloor',1e-20, ...
    'pointStrainConditionFloor',1e-6,'referenceWallIntegral',1, ...
    'referenceWallGradientCentroid',0,'peakPositionRestoringRate',.5,'pinX',0);
wallOps.rescaling=settings;
wallRho=ipm.field.initialDensity(struct('X',x,'Y',zeros(size(x)), ...
    'ny',1,'nx',numel(x)),p);
wallQ=wallRho*Dx';
I=ipm.evolve.wallGradientIntegral(wallQ,wallOps,settings);
P=ipm.evolve.globalSplineMaximum(wallQ,wallOps);
C=ipm.evolve.wallGradientCentroid(wallQ,wallOps);
settings.referenceWallIntegral=I.value;
settings.referenceWallGradientCentroid=C.value;
wallOps.rescaling=settings;
rho=ipm.field.initialDensity(struct('X',g.X,'Y',g.Y,'ny',g.nr,'nx',g.ntheta),p);
cart=struct('originalMetadata',struct('config',struct('physics',p)), ...
    'originalWallOps',wallOps,'scale',scale,'flow',struct('c_l',0,'c_r',0));
reference=struct('wallRho',wallRho,'wallQ',wallQ,'IPC',[I.value,P.value,C.value]);
end
