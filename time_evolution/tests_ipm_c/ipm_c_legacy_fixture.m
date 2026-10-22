function varargout=ipm_c_legacy_fixture(action,varargin)
%IPM_C_LEGACY_FIXTURE Retired analytic/background kernels for test oracles.
% Production solves never call this file. It retains the original uncached
% Poisson and Gaussian split only for historical parity/manufactured checks.
switch char(action)
    case 'background',[varargout{1:nargout}]=background(varargin{:});
    case 'poisson',[varargout{1:nargout}]=poisson(varargin{:});
    case 'rhs',[varargout{1:nargout}]=rhs(varargin{:});
    otherwise,error('ipm_c:LegacyFixtureAction','Unknown legacy fixture action.');
end
end

function bg=background(physics,X,Y,scale)
%POLAR_GAUSSIAN_BACKGROUND Exact original left-zero Gaussian reference field.
% b=e^omega*rho_init(e^-ell*(X-S),e^-ell*Y), where ell=logC_l,
% omega=logC_omega and S=X_shift. Derivatives are with respect to X/Y.
% This is a representation of the same physical reference datum; no state,
% clock, gauge reference, source threshold or coordinate anchor is reset.
assert(isstruct(physics) && isscalar(physics) && ...
    isfield(physics,'initialCondition') && ...
    strcmp(physics.initialCondition,'localized_gaussian_primitive'), ...
    'polar:BackgroundDatum','Use the original localized_gaussian_primitive datum.');
names={'gaussianAmplitude','gaussianCenterX','gaussianWidthX', ...
    'gaussianWidthY','gaussianSkewWeight'};
assert(all(isfield(physics,names)), ...
    'polar:BackgroundParameters','The original Gaussian parameters must all be supplied.');
for k=1:numel(names)
    validateattributes(physics.(names{k}),{'numeric'},{'real','finite','scalar'});
end
A=physics.gaussianAmplitude;a=physics.gaussianCenterX;
sx=physics.gaussianWidthX;sy=physics.gaussianWidthY;w=physics.gaussianSkewWeight;
assert(A>0 && sx>0 && sy>0 && w>=0 && w<=1, ...
    'polar:BackgroundParameters','Gaussian widths/amplitude must be positive and 0<=w<=1.');
assert(isnumeric(X) && isreal(X) && isnumeric(Y) && isreal(Y) && ...
    isequal(size(X),size(Y)) && all(isfinite(X),'all') && all(isfinite(Y),'all'), ...
    'polar:BackgroundCoordinates','Supply matching finite real X/Y arrays.');
assert(isstruct(scale) && isscalar(scale) && ...
    all(isfield(scale,{'logC_l','logC_omega','X_shift'})), ...
    'polar:BackgroundScale','Use the actual logC_l/logC_omega/X_shift state.');
for key={'logC_l','logC_omega','X_shift'}
    validateattributes(scale.(key{1}),{'numeric'},{'real','finite','scalar'});
end
lengthInverse=exp(-scale.logC_l);amplitude=exp(scale.logC_omega);
xp=lengthInverse*(X-scale.X_shift);yp=lengthInverse*Y;
assert(isfinite(lengthInverse) && isfinite(amplitude) && ...
    all(isfinite(xp),'all') && all(isfinite(yp),'all'), ...
    'polar:BackgroundRepresentability','The background coordinate/scale is not representable.');
z=(xp-a)/sx;t=(z-.8)/1.25;s=yp/sy;
[e1,ze1]=gaussian_terms(z);[e2,te2]=gaussian_terms(t);
[ey,sey,s2ey]=gaussian_terms(s);
primitive=(1-w)*erfc(-z)+w*erfc(-t);
primitive= A*sx*sqrt(pi)/2*primitive;
qx=A*((1-w)*e1+(w/1.25)*e2);
qxx=(-2*A/sx)*((1-w)*ze1+(w/1.25^2)*te2);
eyPrime=(-2/sy)*sey;eySecond=(2/sy^2)*s2ey;
gradientScale=amplitude*lengthInverse;
hessianScale=gradientScale*lengthInverse;
bg=struct('rho',amplitude*(primitive.*ey), ...
    'rx',gradientScale*(qx.*ey),'ry',gradientScale*(primitive.*eyPrime), ...
    'rxx',hessianScale*(qxx.*ey),'rxy',hessianScale*(qx.*eyPrime), ...
    'ryy',hessianScale*(primitive.*eySecond), ...
    'definition','exp(omega)*left_zero_Gaussian_primitive(exp(-ell)*(X-S),exp(-ell)*Y)', ...
    'scale',scale,'primitiveConstant','rho_init(-infinity,y)=0');
for key={'rho','rx','ry','rxx','rxy','ryy'}
    assert(all(isfinite(bg.(key{1})),'all'), ...
        'polar:BackgroundNonfinite','An exact background component is not representable.');
end
end

function [e,ze,z2e]=gaussian_terms(z)
% Avoid Inf*0 after exact Gaussian underflow; no nonzero threshold is used.
e=exp(-z.^2);ze=zeros(size(e));z2e=ze;represented=e~=0;
ze(represented)=z(represented).*e(represented);
z2e(represented)=(2*z(represented).^2-1).*e(represented);
end


function [flow,info]=poisson(source,grid,options)
%POLAR_POISSON Half-plane Darcy law on a standard semicircular polar grid.
% -Delta psi=source, psi=0 on the physical wall. Outer values use the
% SAME half-plane image Green kernel, with polar area weights. A supplied
% outerPsi is restricted to manufactured/operator validation.
if nargin<3,options=struct();end
assert(strcmp(grid.domain,'semicircle'),'polar:PoissonDomain', ...
    'This separated Poisson implementation requires a standard r-by-theta grid.');
r=grid.r(:);theta=grid.theta(:)';[nr,nt]=size(source);
assert(isequal(size(source),[numel(r),numel(theta)]) && all(isfinite(source),'all'));
assert(max(abs(theta-linspace(0,pi,nt)))<128*eps(pi),'polar:UniformTheta','Use a uniform angular axis.');
if isfield(options,'outerPsi')
    assert(isfield(options,'manufactured') && options.manufactured, ...
        'polar:BoundaryOverride','Boundary overrides are only for manufactured tests.');
    outer=reshape(options.outerPsi,1,[]);
    greenInfo=struct('kind','manufactured_boundary_override','sourceCount',0);
else
    [outer,greenInfo]=green_outer(source,grid,options);
end
assert(numel(outer)==nt && all(isfinite(outer)));
outer([1,end])=0;
modes=1:nt-2;S=sin(theta(2:end-1)'*modes);
% A nonzero source trace at the physical wall has a discontinuous odd sine
% extension. Lift its endpoint values before the sine solve. This is an
% algebraic Poisson decomposition, not a source mask or changed boundary.
A=sin(theta).^2.*(1+cos(theta))/4;
B=sin(theta).^2.*(1-cos(theta))/4;
Ap=sin(2*theta)/4+(-sin(theta)+3*sin(3*theta))/16;
Bp=sin(2*theta)/4-(-sin(theta)+3*sin(3*theta))/16;
App=cos(2*theta)/2-cos(theta)/16+9*cos(3*theta)/16;
Bpp=cos(2*theta)/2+cos(theta)/16-9*cos(3*theta)/16;
A([1,end])=0;B([1,end])=0;Ap([1,end])=0;Bp([1,end])=0;
a=-r.^2.*source(:,1);b=-r.^2.*source(:,end);
lift=a*A+b*B;liftTheta=a*Ap+b*Bp;
radialA=grid.Drr*a;radialB=grid.Drr*b;
radialA(2:end)=radialA(2:end)+(grid.Dr(2:end,:)*a)./r(2:end);
radialB(2:end)=radialB(2:end)+(grid.Dr(2:end,:)*b)./r(2:end);
lapLift=radialA*A+radialB*B;
lapLift(2:end,:)=lapLift(2:end,:)+(a(2:end)*App+b(2:end)*Bpp)./r(2:end).^2;
liftedSource=source+lapLift;
sourceModes=liftedSource(:,2:end-1)*S*(2/(nt-1));
outerResidual=outer-lift(end,:);
outerModes=outerResidual(2:end-1)*S*(2/(nt-1));
I=2:nr-1;R=spdiags(r,0,nr,nr);
base=-R^2*grid.Drr-R*grid.Dr;
persistent cache
rebuild=isempty(cache) || ~isequal(cache.r,r) || cache.nt~=nt || ~isequal(cache.Dr,grid.Dr);
if rebuild
    cache=struct('r',r,'nt',nt,'Dr',grid.Dr,'factors',{{}},'rowScale',{{}});
    for k=modes
        A=base(I,I)+k*k*speye(nr-2);
        row=max(abs(A),[],2);row=max(full(row),realmin);
        cache.rowScale{k}=row;cache.factors{k}=decomposition(spdiags(1./row,0,nr-2,nr-2)*A,'lu');
    end
end
coeff=zeros(nr,nt-2);coeff(end,:)=outerModes;worst=0;physicalWorst=0;
for k=modes
    b=r(I).^2.*sourceModes(I,k)-base(I,end)*outerModes(k);
    coeff(I,k)=cache.factors{k}\(b./cache.rowScale{k});
    defect=(base(I,I)+k*k*speye(nr-2))*coeff(I,k)-b;
    worst=max(worst,norm(defect,inf)/max(norm(b,inf),realmin));
    physicalDefect=defect./r(I).^2;
    physicalWorst=max(physicalWorst,norm(physicalDefect,inf)/max(1,norm(sourceModes(I,k),inf)));
end
psi=coeff*sin(modes'*theta)+lift;psi(:,[1,end])=0;psi(1,:)=0;
psiTheta=coeff*(modes'.*cos(modes'*theta))+liftTheta;
psiR=grid.Dr*psi;
ur=zeros(nr,nt);ur(2:end,:)=-psiTheta(2:end,:)./r(2:end);
ut=psiR;
ux=ur.*cos(theta)-ut.*sin(theta);
uy=ur.*sin(theta)+ut.*cos(theta);
% Only the m=1 sine coefficient contributes to the Cartesian first jet.
% Use the same seven-point radial derivative at r=0; a quadratic 2-D fit
% would reduce this otherwise polynomial-exact manufacture to second order.
originUx=-(grid.Dr(1,:)*coeff(:,1));
ux(1,:)=originUx;uy(1,:)=0;
ur(1,:)=originUx*cos(theta);ut(1,:)=-originUx*sin(theta);
uy(:,[1,end])=0;
assert(all(isfinite(ux),'all') && all(isfinite(uy),'all'));
flow=struct('psi',psi,'u1',ux,'u2',uy,'ur',ur,'utheta',ut,'source',source);
info=struct('definition','minus_polar_laplacian_psi_equals_rho_X', ...
    'angularMethod','Dirichlet sine collocation','radialMethod','original seven-point physical-axis finite differences', ...
    'green',greenInfo,'maximumModalRelativeResidual',worst,'factorCacheRebuilt',rebuild, ...
    'maximumUnweightedModalResidualScaledBySource',physicalWorst, ...
    'nonzeroWallSourceLiftApplied',true,'liftIsExactPoissonEquationDecomposition',true, ...
    'wallNormalVelocityInfinity',max(abs(uy(:,[1,end])),[],'all'), ...
    'originIsSinglePhysicalPoint',true,'nativeCartesianOperatorUnchanged',true, ...
    'sameHalfDiscreteOperatorClaimed',false);
end

function [outer,info]=green_outer(q,g,options)
% Same signed half-plane image kernel as original greenBoundary.m:105-142.
% Only the geometry/quadrature and bin coordinates are changed.
tol=opt(options,'greenSourceTolerance',1e-10);cap=opt(options,'greenMaxSources',2048);
r=g.r(:);t=g.theta(:)';nr=numel(r);nt=numel(t);
wr=[r(2)-r(1);r(3:end)-r(1:end-2);r(end)-r(end-1)]/2;
wt=[t(2)-t(1),t(3:end)-t(1:end-2),t(end)-t(end-1)]/2;
area=(r.*wr)*wt;
interior=false(nr,nt);interior(2:end-1,2:end-1)=true;
mx=max(abs(q(interior)));keep=interior & abs(q)>=tol*mx;
if mx==0,outer=zeros(1,nt);info=struct('sourceCount',0,'candidateCount',0);return;end
[ri,ti]=find(keep);xs=g.X(keep);ys=g.Y(keep);strength=q(keep).*area(keep);count=numel(strength);
if count>cap
    binsR=max(1,floor(sqrt(floor(cap/2)*nr/nt)));binsT=max(1,floor(floor(cap/2)/binsR));
    br=min(floor((ri-1)*binsR/nr)+1,binsR);bt=min(floor((ti-1)*binsT/nt)+1,binsT);
    bins=br+(bt-1)*binsR;n=binsR*binsT;
    [xp,yp,sp]=aggregate(bins,xs,ys,max(strength,0),n,1);
    [xn,yn,sn]=aggregate(bins,xs,ys,max(-strength,0),n,-1);
    xs=[xp;xn];ys=[yp;yn];strength=[sp;sn];
end
outer=zeros(1,nt);tx=g.X(end,:)';ty=g.Y(end,:)';
for a=1:64:nt
    z=a:min(a+63,nt);distance=(tx(z)-xs').^2+(ty(z)-ys').^2;
    kernel=log1p(4*(ty(z)*ys')./distance)/(4*pi);
    outer(z)=(kernel*strength)';
end
outer([1,end])=0;
info=struct('candidateCount',count,'sourceCount',numel(strength),'sourceThreshold',tol, ...
    'maxSources',cap,'areaWeight','r dr dtheta positive trapezoidal', ...
    'imageKernel','log1p(4*y*ys/distanceSquared)/(4*pi)', ...
    'sourceMaskIsOnlyOriginalGreenCompressionPolicy',true);
end

function [x,y,s]=aggregate(b,xs,ys,w,n,sgn)
w0=accumarray(b,w,[n,1],@sum,0);active=w0>0;
x0=accumarray(b,w.*xs,[n,1],@sum,0);y0=accumarray(b,w.*ys,[n,1],@sum,0);
x=x0(active)./w0(active);y=y0(active)./w0(active);s=sgn*w0(active);
end

function v=opt(s,k,d)
v=d;if isfield(s,k),v=s.(k);end
end


function [F,flow,info]=rhs(rho,grid,cart,options)
%POLAR_BACKGROUND_RHS Complete IPM RHS using an exact Gaussian background.
% Input/output remain TOTAL rho and F. The representation rho=b+v computes
% source=b_X+D_X(v) and B=T(v)-Vx*b_X-Vy*b_Y, then F=B+c_omega*rho.
% The background is never removed from the physical source, Poisson problem,
% full gauge, amplitude reaction, or clocks. Original fixed Cartesian wall
% observables sample remainders plus their analytic wall background.
if nargin<4,options=struct();end
assert(isstruct(options)&&isscalar(options),'polar:BackgroundOptions');
assert(isfield(cart,'originalMetadata')&&isfield(cart.originalMetadata,'config')&& ...
    isfield(cart.originalMetadata.config,'physics')&&isfield(cart,'scale'), ...
    'polar:BackgroundConfiguration','The original Gaussian physics and current complete scale are required.');
assert(isequal(size(rho),[grid.nr,grid.ntheta])&&all(isfinite(rho),'all')&& ...
    all(rho(1,:)==rho(1,1)),'polar:BackgroundDensity', ...
    'Input must be one finite total density with a single shared origin.');
cfg=cart.originalMetadata.config;
bg=ipm_c_legacy_fixture('background',cfg.physics,grid.X,grid.Y,cart.scale);
v=rho-bg.rho;
primaryRemainder=false;
if isfield(options,'remainderState')
    v=options.remainderState;primaryRemainder=true;
    assert(isequal(size(v),size(rho))&&all(isfinite(v),'all')&& ...
        all(v(1,:)==v(1,1)), 'polar:PrimaryRemainder','Supply one finite single-origin remainder.');
    assert(isequaln(rho,bg.rho+v),'polar:PrimaryRemainderIdentity', ...
        'The complete state must be assembled from this same background and primary remainder.');
end
dv=ipm_c.field.derivatives(v,grid,options);
source=bg.rx+dv.rx;
assert(all(isfinite(source),'all'),'polar:BackgroundSource','The complete source must remain finite.');
poissonOptions=options;
if isfield(cfg,'elliptic')
    for key={'greenSourceTolerance','greenMaxSources'}
        if isfield(cfg.elliptic,key{1}),poissonOptions.(key{1})=cfg.elliptic.(key{1});end
    end
end
[physical,pinfo]=ipm_c_legacy_fixture('poisson',source,grid,poissonOptions);
method='linear';if isfield(options,'measurementInterpolation'),method=options.measurementInterpolation;end
wallX=cart.originalWallOps.x(:)';
bgwall=ipm_c_legacy_fixture('background',cfg.physics,wallX,zeros(size(wallX)),cart.scale);
measurementOptions=struct('interpolation',method,'fieldOffset',bg.rho,'wallOffset',bgwall.rho);
measurement=ipm_c.gauge.measure(rho,grid.r,cart.originalWallOps,measurementOptions);
signedX=[-flipud(grid.r(2:end));grid.r];
wallU=[flipud(physical.u1(2:end,end));physical.u1(:,1)];
Uref=interp1(signedX,wallU,measurement.C,'linear');
guess=[cart.flow.c_l,cart.flow.c_r];
if isfield(options,'gaugeGuess'),guess=options.gaugeGuess;end
gaugeOptions=options;gaugeOptions.measurementOptions=measurementOptions;
gaugeOptions.rateMeasurementOffsets=@rate_offsets;
transportOptions=options;
if isfield(cfg,'transport'),transportOptions.transport=cfg.transport;end
[F,rates,ginfo]=ipm_c.gauge.solve(rho,grid.r,cart.originalWallOps,@base,Uref,guess,gaugeOptions);
vx=(physical.u1-Uref)+rates.c_l*grid.X+rates.delta;
vy=physical.u2+rates.c_l*grid.Y;vy(:,[1,end])=0;
[backgroundRate,~]=rate_offsets(rates.c_l,rates.delta,rates.c_omega);
backgroundRateX=(rates.c_omega-rates.c_l)*bg.rx- ...
    (rates.c_l*grid.X+rates.c_r).*bg.rxx-rates.c_l*grid.Y.*bg.rxy;
% Compute the primary remainder equation directly. Recovering it by
% subtracting two total fields loses tiny increments and pollutes D_X F.
rateRemainder=ipm_c.field.transport(v,vx,vy,grid,transportOptions)+rates.c_omega*v- ...
    physical.u1.*bg.rx-physical.u2.*bg.ry;
directF=backgroundRate+rateRemainder;
assemblyDifference=max(abs(F-directF),[],'all');
assert(assemblyDifference<=512*eps(max([1;abs(F(:));abs(directF(:))])), ...
    'polar:BackgroundAssemblyIdentity','The direct remainder equation changed the complete PDE.');
F=directF;
auditOptions=measurementOptions;auditOptions.rateField=F;
[auditOptions.rateFieldOffset,auditOptions.wallRateOffset]= ...
    rate_offsets(rates.c_l,rates.delta,rates.c_omega);
directAudit=ipm_c.gauge.measure(rho,grid.r,cart.originalWallOps,auditOptions);
directResidual=directAudit.rate.normalized+ginfo.target;
assert(max(abs(directResidual))<=ginfo.auditTolerance, ...
    'polar:BackgroundDirectFullAudit','The returned complete RHS failed original I/P/C.');
ginfo.fullMeasurement=directAudit;ginfo.auditResidual=directResidual;
ginfo.directRemainderAssemblyDifferenceInfinity=assemblyDifference;
dRate=ipm_c.field.derivatives(rateRemainder,grid,options);
fullQ=dRate.rx+backgroundRateX;
assert(all(isfinite(fullQ),'all')&&all(isfinite(backgroundRateX),'all'), ...
    'polar:BackgroundFullQ','The full derivative of the complete F must remain finite.');

% Match the ordinary polar RHS diagnostics; these rates are not a proved
% stability bound and do not introduce any time reparameterization.
logical=linspace(0,1,grid.nr)';metric=ipm.mesh.fdMatrix(logical,1,7)*grid.r;
assert(all(metric>0),'polar:RadialMetric','Mapped radial WENO metric must stay positive.');
angularSpacing=pi/(grid.ntheta-1);positive=2:grid.nr;
vr=vx.*grid.cosTheta+vy.*grid.sinTheta;
vt=-vx.*grid.sinTheta+vy.*grid.cosTheta;
rateField=zeros(size(rho));
rateField(positive,:)=abs(vr(positive,:))./(metric(positive)/(grid.nr-1))+ ...
    abs(vt(positive,:))./(grid.r(positive)*angularSpacing);
signedWallX=signedX';wallLogical=linspace(-1,1,numel(signedWallX));
wallMetric=(ipm.mesh.fdMatrix(wallLogical,1,7)*signedWallX')';
assert(all(wallMetric>0),'polar:WallMetric','The joined physical wall metric must be positive.');
[~,wallNode]=min(abs(signedWallX));
signedV=[flipud(vx(2:end,end));vx(:,1)]';
wallRates=abs(signedV)./(wallMetric/(grid.nr-1));
rateField(:,1)=wallRates(grid.nr:end)';rateField(:,end)=flipud(wallRates(1:grid.nr)');
rateField(1,:)=wallRates(wallNode);
[cfl,index]=max(rateField(:));[ri,ti]=ind2sub(size(rho),index);
flow=physical;flow.source=source;
flow.c_l=rates.c_l;flow.c_omega=rates.c_omega;flow.c_r=rates.c_r;
flow.generatorTranslationRate=rates.delta;flow.generatorReferenceVelocity=Uref;
flow.transportU1=vx;flow.transportU2=vy;
flow.fullQ=fullQ;flow.backgroundRate=backgroundRate;flow.backgroundRateX=backgroundRateX;
flow.remainderRate=rateRemainder;flow.remainderState=v;
flow.backgroundSource=bg.rx;flow.remainderSource=dv.rx;
info=struct('poisson',pinfo,'gauge',ginfo,'maximumMappedPolarCflRate',cfl, ...
    'originWallMappedTransportRate',wallRates(wallNode), ...
    'diagnosticTransportRateIncludingOrigin',max(cfl,wallRates(wallNode)), ...
    'transportRateIsNotAProvedStabilityBound',true, ...
    'cflPosition',struct('r',grid.r(ri),'theta',grid.theta(ti),'X',grid.X(ri,ti),'Y',grid.Y(ri,ti)), ...
    'singleOriginCartesianEquation',true,'completePolarPdeEvaluated',true, ...
    'sameContinuousEquation',true,'originalWallFunctionalGridPreserved',true, ...
    'sameSemiDiscreteEquation',false,'referencesReset',false,'physicalTimeReparameterized',false, ...
    'poleRegularization',false, ...
    'exactGaussianBackground',true,'backgroundIsPhysicalSourceComponent',true, ...
    'backgroundRemovedFromPde',false,'fullAuditIncludesBackgroundTimeRate',true, ...
    'sourceDefinition','analytic b_X + polar D_X(rho-b)', ...
    'fullQDefinition','polar D_X(F-Fbg) + analytic (Fbg)_X', ...
    'inputAndOutputRemainTotalFields',true,'remainderInfinity',max(abs(v),[],'all'));
info.primaryRemainderSupplied=primaryRemainder;
info.remainderRateRecoveredBySubtraction=false;
info.completeRhsAssemblyDifferenceInfinity=assemblyDifference;

    function B=base(cl,delta)
        vx0=(physical.u1-Uref)+cl*grid.X+delta;
        vy0=physical.u2+cl*grid.Y;vy0(:,[1,end])=0;
        B=ipm_c.field.transport(v,vx0,vy0,grid,transportOptions)-vx0.*bg.rx-vy0.*bg.ry;
        assert(all(isfinite(B),'all'),'polar:NonfiniteRhs','The complete background-plus-remainder base must be finite.');
    end
    function [polarOffset,wallOffset]=rate_offsets(cl,delta,cw)
        cr=delta-Uref;
        polarOffset=-(cl*grid.X+cr).*bg.rx-cl*grid.Y.*bg.ry+cw*bg.rho;
        wallOffset=-(cl*wallX+cr).*bgwall.rx+cw*bgwall.rho;
    end
end
