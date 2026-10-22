function report=test_bspline_poisson()
%TEST_BSPLINE_POISSON Exact polynomials and physical-domain refinement.
% Spatial accuracy is checked against analytic psi/u, independently of the
% Galerkin linear-system residual. No evolved trajectory is certified here.
g=ipm_c.mesh.build(struct('nr',65,'ntheta',65,'rmax',1,'firstSpacing',1/64));
o=struct('manufactured',true,'outerPsi',zeros(1,g.ntheta));
source=8*g.Y;exact=g.Y.*(1-g.X.^2-g.Y.^2);
[f,info]=ipm_c.spatial.poisson(source,g,o);
psiError=max(abs(f.psi-exact),[],'all');
uError=max(abs([f.u1+(1-g.X.^2-3*g.Y.^2),f.u2+2*g.X.*g.Y]),[],'all');
assert(info.bsplineOrder==4&&info.bsplineDegree==3&& ...
    psiError<1e-9&&uError<1e-8,'polar:BsplinePolynomial');
[warm,winfo]=ipm_c.spatial.poisson(source,g,o);
assert(~winfo.factorCacheRebuilt&&isequaln(warm.psi,f.psi)&& ...
    isequaln(warm.u1,f.u1)&&isequaln(warm.u2,f.u2),'polar:BsplineCacheReuse');
[zero,zinfo]=ipm_c.spatial.poisson(zeros(size(source)),g);
assert(all(zero.psi==0,'all')&&all(zero.u1==0,'all')&&all(zero.u2==0,'all'), ...
    'polar:BsplineZero');
forcing=-2+2*g.X.^2+14*g.Y.^2;
[wall,wallInfo]=ipm_c.spatial.poisson(forcing,g,o);
exactWall=g.Y.^2.*(1-g.X.^2-g.Y.^2);
wallPsiError=max(abs(wall.psi-exactWall),[],'all');
wallU1=-(2*g.Y-2*g.X.^2.*g.Y-4*g.Y.^3);wallU2=-2*g.X.*g.Y.^2;
wallUError=max(abs([wall.u1-wallU1,wall.u2-wallU2]),[],'all');
assert(wallPsiError<2e-8&&wallUError<2e-7,'polar:BsplineNonzeroWallSource');
single=zeros(size(source));single(22,19)=1;
% This underresolved single-node source checks only the boundary kernel,
% not spatial convergence. Production keeps the .02 strong-defect gate.
[point,pinfo]=ipm_c.spatial.poisson(single,g,struct('manufactured',true, ...
    'maxStrongResidual',.05,'greenSourceTolerance',1e-10,'greenMaxSources',2048));
strength=g.r(22)*(g.r(23)-g.r(21))/2*(g.theta(20)-g.theta(18))/2;
xs=g.X(22,19);ys=g.Y(22,19);
boundary=log1p(4*g.Y(end,:)*ys./((g.X(end,:)-xs).^2+(g.Y(end,:)-ys).^2))/(4*pi)*strength;
greenError=max(abs(point.psi(end,:)-boundary));
assert(greenError<1e-12,'polar:BsplineGreenBoundary');
% Regression for the former nonlocal source-spline failure. This deliberate
% legacy axis has a steep geometric tail: a global not-a-knot source spline
% spreads an isolated core ring to the far boundary. The single sine mode
% has zero wall source, isolating source-load locality from a sharply varying
% wall lift (covered separately by the nonzero-wall polynomial manufacture).
% The discontinuous radial ring tests locality, not smooth convergence;
% only this manufactured test uses the explicit finite strong-defect bound 10.
legacy=ipm_c.mesh.build(struct('nr',385,'ntheta',128,'rmax',1e10, ...
    'firstSpacing',.03,'radialMap','sinh'));
xi=linspace(0,1,legacy.nr)';
legacyAxis=11.52*xi+(1e10-11.52)*(sinh(100*xi)/sinh(100)).^3;
legacyAxis(1)=0;legacyAxis(2)=.03;legacyAxis(end)=1e10;
legacy=ipm_c.remesh.radialGeometry(legacy,legacyAxis);
[~,ringRow]=min(abs(legacy.r-10.2));
ringSource=zeros(legacy.nr,legacy.ntheta);ringSource(ringRow,:)=sin(legacy.theta);
ringSource(:,[1,end])=0;
ringOptions=struct('manufactured',true,'outerPsi',zeros(1,legacy.ntheta), ...
    'maxStrongResidual',10);
[ring,ringInfo]=ipm_c.spatial.poisson(ringSource,legacy,ringOptions);
far=legacy.R>1e6;
farVelocity=max(hypot(ring.u1(far),ring.u2(far)));
assert(farVelocity<1e-6&&~ringInfo.sourceReconstructionUsesGlobalInverse&& ...
    ringInfo.sourceReconstructionMaximumSupportNodes<=6, ...
    'polar:BsplineSourceLocality', ...
    'A compact source on a steep radial tail must not create a spurious far-field velocity.');
[ringWarm,ringWarmInfo]=ipm_c.spatial.poisson(ringSource,legacy,ringOptions);
assert(~ringWarmInfo.factorCacheRebuilt&& ...
    ringWarmInfo.sourceReconstructionGeometryCacheReused&& ...
    isequaln(ringWarm.psi,ring.psi)&&isequaln(ringWarm.u1,ring.u1)&& ...
    isequaln(ringWarm.u2,ring.u2),'polar:BsplineSourceLocalityCache');
legacyGaps=diff(legacy.r);
sourceLocality=struct('kind','compact_ring_on_legacy_steep_tail_v1', ...
    'sourceRadius',legacy.r(ringRow),'sourceRow',ringRow, ...
    'maximumAdjacentGapRatio',max(legacyGaps(2:end)./legacyGaps(1:end-1)), ...
    'farRadiusThreshold',1e6,'farVelocityInfinity',farVelocity, ...
    'velocityInfinity',max(hypot(ring.u1,ring.u2),[],'all'), ...
    'strongResidual',ringInfo.maximumOffGridStrongResidualScaledBySource, ...
    'manufacturedStrongResidualLimit',ringOptions.maxStrongResidual, ...
    'factorCacheReused',~ringWarmInfo.factorCacheRebuilt, ...
    'smoothSpatialAccuracyTest',false);
records=struct([]);
for n=[161,321,641]
    cfg=ipm_c.config();settings=cfg.grid;settings.nr=n;settings.ntheta=128;
    settings.firstSpacing=.06*160/(n-1);
    grid=ipm_c.mesh.build(settings);R2=grid.X.^2+grid.Y.^2;E=exp(-R2);
    exact=grid.Y.*E;source=(8-4*R2).*grid.Y.*E;
    options=struct('manufactured',true,'outerPsi',zeros(1,grid.ntheta));
    timer=tic;[v,details]=ipm_c.spatial.poisson(source,grid,options);seconds=toc(timer);
    expectedU1=-(1-2*grid.Y.^2).*E;expectedU2=-2*grid.X.*grid.Y.*E;
    core=grid.R<=2;
    errPsi=max(abs(v.psi(core)-exact(core)));
    errU=max([max(abs(v.u1(core)-expectedU1(core))),max(abs(v.u2(core)-expectedU2(core)))]);
    assert(all(isfinite(v.psi),'all')&&all(isfinite(v.u1),'all')&& ...
        details.linearSolveBackwardError<=1e-10&&grid.r(end)==1e10, ...
        'polar:BsplineFullDomain');
    row=struct('nr',n,'ntheta',grid.ntheta,'h',settings.firstSpacing, ...
        'rEnd',grid.r(end),'corePsiInfError',errPsi,'coreVelocityInfError',errU, ...
        'coldPoissonSeconds',seconds,'poisson',details);
    if isempty(records),records=row;else,records(end+1)=row;end %#ok<AGROW>
end
psiOrder=log2([records(1:2).corePsiInfError]./[records(2:3).corePsiInfError]);
velocityOrder=log2([records(1:2).coreVelocityInfError]./[records(2:3).coreVelocityInfError]);
assert(all(psiOrder>3.2)&&all(velocityOrder>2.5), ...
    'polar:BsplineRefinement','Cubic B-spline fields must converge against the analytic solution.');
report=struct('kind','polar_cubic_bspline_galerkin_exact_and_refined_v1','passed',true, ...
    'bsplineOrder',4,'bsplineDegree',3,'psiPolynomialError',psiError, ...
    'velocityPolynomialError',uError,'nonzeroWallPsiError',wallPsiError, ...
    'nonzeroWallVelocityError',wallUError,'imageGreenBoundaryError',greenError, ...
    'polynomialInfo',info,'wallInfo',wallInfo,'zeroInfo',zinfo,'greenInfo',pinfo, ...
    'sourceLocalityFaultRegression',sourceLocality, ...
    'refinement',records,'psiOrders',psiOrder,'velocityOrders',velocityOrder, ...
    'factorCacheReused',true,'timeEvolutionTested',false,'longTrajectoryQualified',false);
end
