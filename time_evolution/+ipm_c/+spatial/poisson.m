function [flow,info]=poisson(source,grid,options)
%POISSON Cached cubic B-spline Galerkin/FFT polar Darcy operator.
% -Delta psi=source on the half disk, with zero physical-wall trace and
% the original image-Green outer values. Radial trial functions are order-4
% (degree-3), C2 physical-r B-splines. No finite-difference Poisson operator
% or iterative linear solve is used. Velocity differentiates this SAME psi.
if nargin<3,options=struct();end
if nargin==1&&ischar(source)&&strcmp(source,'transform-test')
    flow=transform_test();info=struct('physicalPdeCalls',0);return;
end
assert(~isfield(options,'bsplineOrder')||isequal(options.bsplineOrder,4), ...
    'polar:BSplineOrder','This Poisson operator requires order 4 (degree 3) B-splines.');
assert(~isfield(options,'radialDiscretization')|| ...
    strcmp(options.radialDiscretization,'cubic_bspline_galerkin'), ...
    'polar:BSplineDiscretization','Use the physical-r cubic B-spline Galerkin discretization.');
assert(strcmp(grid.domain,'semicircle'),'polar:PoissonDomain', ...
    'The separated Poisson implementation requires a standard r-by-theta grid.');
r=grid.r(:);theta=grid.theta(:)';[nr,nt]=size(source);
assert(isequal(size(source),[numel(r),numel(theta)])&&all(isfinite(source),'all'));
assert(max(abs(theta-linspace(0,pi,nt)))<128*eps(pi),'polar:UniformTheta','Use a uniform angular axis.');
persistent cache
rebuild=isempty(cache)||~isequal(cache.r,r)||~isequal(cache.theta,theta);
if rebuild,cache=prepare_geometry(r,theta);end
if isfield(options,'outerPsi')
    assert(isfield(options,'manufactured')&&options.manufactured, ...
        'polar:BoundaryOverride','Boundary overrides are only for manufactured tests.');
    outer=reshape(options.outerPsi,1,[]);
    greenInfo=struct('kind','manufactured_boundary_override','sourceCount',0);
else
    [outer,greenInfo]=green_outer(source,grid,options,cache.green);
end
assert(numel(outer)==nt&&all(isfinite(outer)));outer([1,end])=0;
modes=cache.modes;nm=numel(modes);I=cache.interior;
% One local cubic Hermite reconstruction serves every modal source and both
% wall traces. Its nodal slopes use a cached five-point physical-r map. A
% local source perturbation therefore cannot create a tail across the box.
sourceModes=sine_analysis(source(:,2:end-1));
sourceData=[sourceModes,source(:,[1,end])];
sourceSlopes=cache.sourceDerivative*sourceData;
qg=cache.q.value*sourceData+cache.q.slope*sourceSlopes;
wallQ=qg(:,nm+(1:2));
wallD=cache.q.first*sourceData(:,nm+(1:2))+cache.q.firstSlope*sourceSlopes(:,nm+(1:2));
% L=a(r)A(theta)+b(r)B(theta), a=-r^2*q_wall. Differentiate that
% reconstruction analytically; do not spline the already multiplied a data.
aprime=-2*cache.q.r.*wallQ-cache.q.r.^2.*wallD;
aOverR=-cache.q.r.*wallQ;
load=cache.loadSource*qg(:,1:nm) ...
    -(cache.loadDerivative*aprime(:,1))*cache.Ahat ...
    -(cache.loadDerivative*aprime(:,2))*cache.Bhat ...
    +(cache.loadOverR*aOverR(:,1))*cache.Apphat ...
    +(cache.loadOverR*aOverR(:,2))*cache.Bpphat;
a=-r.^2.*source(:,1);b=-r.^2.*source(:,end);
lift=a*cache.A+b*cache.B;
outerModes=sine_analysis(outer(2:end-1)-lift(end,2:end-1));
coefficients=zeros(cache.basisCount-1,nm);coefficients(end,:)=outerModes;
rhs=load-cache.outerColumn.*outerModes;
for k=modes
    coefficients(I,k)=cache.solveScale(:,k).*(cache.factors{k}\(rhs(:,k).*cache.solveScale(:,k)));
end
% The weak algebraic residual checks the actual assembled Galerkin solve.
% It is explicitly distinct from the off-grid unweighted strong PDE audit.
c=coefficients(I,:);
defect=zeros(size(c));backwardDenominator=zeros(size(c));
for k=modes
    defect(:,k)=cache.operators{k}*c(:,k)-rhs(:,k);
    backwardDenominator(:,k)=cache.absoluteOperators{k}*abs(c(:,k))+abs(rhs(:,k));
end
backwardError=max(max(abs(defect),[],1)./max(max(backwardDenominator,[],1),realmin));
componentwiseBackwardError=max(abs(defect)./max(backwardDenominator,realmin),[],'all');
weakRelative=max(max(abs(defect),[],1)./max(max(abs(rhs),[],1),realmin));
assert(isfinite(backwardError)&&backwardError<=1e-10,'polar:BSplineLinearSolve', ...
    'The cached direct B-spline Galerkin solve failed its backward-error gate.');
nativeModes=cache.Bnodes*coefficients;
psi=sine_synthesis(nativeModes)+lift;psi(:,[1,end])=0;psi(1,:)=0;
wallNodeD=sourceSlopes(:,nm+(1:2));
wallAPrime=-2*r.*source(:,[1,end])-r.^2.*wallNodeD;
psiR=sine_synthesis(cache.Dnodes*coefficients)+wallAPrime(:,1)*cache.A+wallAPrime(:,2)*cache.B;
psiTheta=cosine_synthesis(nativeModes.*modes)+a*cache.Ap+b*cache.Bp;
ur=zeros(nr,nt);ur(2:end,:)=-psiTheta(2:end,:)./r(2:end);
ut=psiR;
ux=ur.*cache.cosTheta-ut.*cache.sinTheta;
uy=ur.*cache.sinTheta+ut.*cache.cosTheta;
% The Cartesian first jet is the m=1 radial slope. Higher sine modes need
% not have an exact zero endpoint slope in a finite Galerkin approximation;
% report that defect, and retain one Cartesian velocity at the physical pole.
originSlopes=cache.Dnodes(1,:)*coefficients;originUx=-originSlopes(1);
ux(1,:)=originUx;uy(1,:)=0;
ur(1,:)=originUx*cache.cosTheta;ut(1,:)=-originUx*cache.sinTheta;
uy(:,[1,end])=0;
[strong,coreStrong,nearStrong,strongInfinity]=strong_audit(coefficients,sourceData,sourceSlopes,cache);
assert(all(isfinite(ux),'all')&&all(isfinite(uy),'all')&&isfinite(strong), ...
    'polar:NonfiniteBSplinePoisson','The complete spline solution and independent strong residual must be finite.');
maxStrongResidual=opt(options,'maxStrongResidual',.02);
validateattributes(maxStrongResidual,{'numeric'},{'scalar','real','finite','positive'});
assert(strong<=maxStrongResidual,'polar:BSplineStrongResidual', ...
    'The off-grid strong residual %.6g exceeds the configured spatial-defect limit %.6g.', ...
    strong,maxStrongResidual);
flow=struct('psi',psi,'u1',ux,'u2',uy,'ur',ur,'utheta',ut,'source',source);
info=struct('definition','minus_polar_laplacian_psi_equals_rho_X', ...
    'discretizationContract','physical_r_cubic_bspline_galerkin_fft_v1', ...
    'angularMethod','Dirichlet sine collocation via equivalent FFT odd/even extensions', ...
    'radialMethod','physical-r order-4 degree-3 C2 B-spline Galerkin', ...
    'bsplineOrder',4,'bsplineDegree',3,'splineOrder',4,'splineDegree',3,'radialBasisCount',cache.basisCount, ...
    'radialInteriorUnknownCount',numel(I),'quadraturePointsPerSpan',8, ...
    'sourceReconstruction','physical-r local C1 cubic Hermite; cached five-point nodal derivative map', ...
    'sourceReconstructionStencilWidth',5,'sourceReconstructionMaximumSupportNodes',6, ...
    'sourceReconstructionUsesGlobalInverse',false, ...
    'sourceReconstructionInfinityAmplification',max(abs(qg),[],'all')/max(max(abs(sourceData),[],'all'),realmin), ...
    'green',greenInfo,'factorCacheRebuilt',rebuild,'sourceReconstructionGeometryCacheReused',~rebuild, ...
    'linearSolveBackwardError',backwardError,'weakGalerkinRelativeResidual',weakRelative, ...
    'linearSolveComponentwiseBackwardError',componentwiseBackwardError, ...
    'maximumModalRelativeResidual',weakRelative, ...
    'maximumOffGridStrongResidualScaledBySource',strong, ...
    'coreOffGridStrongResidualScaledBySource',coreStrong, ...
    'nearOriginOffGridStrongResidualScaledBySource',nearStrong, ...
    'offGridStrongResidualInfinity',strongInfinity,'strongAuditCoreRadius',cache.auditCoreRadius, ...
    'maximumAllowedOffGridStrongResidualScaledBySource',maxStrongResidual, ...
    'strongAuditSites','midpoint of every physical radial knot interval, complete native angular axis', ...
    'strongResidualDefinition','unweighted -Delta(psi)-local_Hermite_reconstructed_source; analytic trial/lift derivatives', ...
    'strongResidualIsAlgebraicResidual',false,'strongResidualIsAccuracyCertificate',false, ...
    'originHigherModeRadialSlopeInfinity',max(abs(originSlopes(2:end))), ...
    'originHigherModeSlopeRelativeToFirst',max(abs(originSlopes(2:end)))/max(abs(originSlopes(1)),realmin), ...
    'originVelocityQualification','single Cartesian value uses mode 1; unresolved higher-mode endpoint slopes remain spatial approximation errors', ...
    'nonzeroWallSourceLiftApplied',true,'liftIsExactPoissonEquationDecomposition',true, ...
    'liftRadialDerivativeMethod','analytic derivative of -r^2 times reconstructed wall source', ...
    'wallNormalVelocityInfinity',max(abs(uy(:,[1,end])),[],'all'), ...
    'originIsSinglePhysicalPoint',true,'nativeCartesianOperatorUnchanged',false, ...
    'velocityDerivativesConsistentWithSpline',true, ...
    'sameHalfDiscreteOperatorClaimed',false,'iterativeLinearSolverUsed',false);
end

function cache=prepare_geometry(r,theta)
% One retained geometry; rebuilding a remeshed axis releases its predecessor.
nr=numel(r);nt=numel(theta);modes=1:nt-2;
knots=[repmat(r(1),1,4),r(2:end-1)',repmat(r(end),1,4)];nb=numel(knots)-4;
active=2:nb;I=1:nb-2;
[Bnodes,Dnodes]=basis_matrices(knots,r);
[x,w]=gauss8();h=diff(r);
rq=reshape((r(1:end-1)+h*(1+x')/2)',[],1);
wq=reshape((h*(w'/2))',[],1);
[Bq,Dq]=basis_matrices(knots,rq);Bq=Bq(:,active);Dq=Dq(:,active);
Wr=spdiags(wq.*rq,0,numel(rq),numel(rq));
WOverR=spdiags(wq./rq,0,numel(rq),numel(rq));
K=Dq'*Wr*Dq;M=Bq'*WOverR*Bq;
K=(K+K')/2;M=(M+M')/2;
% The origin B-spline is excluded before assembling M: it has infinite
% m>=1 angular energy. Every retained basis vanishes at r=0.
cache=struct('r',r,'theta',theta,'h',h,'modes',modes,'interior',I,'basisCount',nb, ...
    'Bnodes',Bnodes(:,active),'Dnodes',Dnodes(:,active), ...
    'K',K(I,I),'M',M(I,I), ...
    'cosTheta',cos(theta),'sinTheta',sin(theta), ...
    'factors',{cell(1,numel(modes))},'operators',{cell(1,numel(modes))}, ...
    'absoluteOperators',{cell(1,numel(modes))}, ...
    'solveScale',zeros(numel(I),numel(modes)), ...
    'outerColumn',K(I,end)+M(I,end)*(modes.^2));
cache.A=sin(theta).^2.*(1+cos(theta))/4;
cache.B=sin(theta).^2.*(1-cos(theta))/4;
cache.Ap=sin(2*theta)/4+(-sin(theta)+3*sin(3*theta))/16;
cache.Bp=sin(2*theta)/4-(-sin(theta)+3*sin(3*theta))/16;
cache.App=cos(2*theta)/2-cos(theta)/16+9*cos(3*theta)/16;
cache.Bpp=cos(2*theta)/2+cos(theta)/16-9*cos(3*theta)/16;
cache.A([1,end])=0;cache.B([1,end])=0;cache.Ap([1,end])=0;cache.Bp([1,end])=0;
cache.Ahat=sine_analysis(cache.A(2:end-1));cache.Bhat=sine_analysis(cache.B(2:end-1));
cache.Apphat=sine_analysis(cache.App(2:end-1));cache.Bpphat=sine_analysis(cache.Bpp(2:end-1));
cache.loadSource=Bq(:,I)'*Wr;
cache.loadDerivative=Dq(:,I)'*Wr;
cache.loadOverR=Bq(:,I)'*spdiags(wq,0,numel(rq),numel(rq));
for k=modes
    matrix=cache.K+k*k*cache.M;
    cache.operators{k}=matrix;
    cache.absoluteOperators{k}=abs(matrix);
    scale=1./sqrt(full(diag(matrix)));
    assert(all(isfinite(scale))&&all(scale>0),'polar:BSplineEnergy','The modal spline energy must have positive diagonal.');
    cache.solveScale(:,k)=scale;
    diagonal=spdiags(scale,0,numel(I),numel(I));
    cache.factors{k}=decomposition(diagonal*matrix*diagonal,'lu');
end
% Source interpolation is local, independently of the global C2 B-spline
% trial space for psi. A global source spline on geometric tails can have
% tiny nodal values but a large r*q between nodes. Five-point slopes retain
% fourth-order source interpolation without a nonlocal curvature inverse.
cache.sourceDerivative=ipm.mesh.fdMatrix(r,1,5);
cache.q=source_maps(r,rq);cache.q.r=rq;
ra=r(1:end-1)+h/2;
[Ba,Da,D2a]=basis_matrices(knots,ra);
cache.audit=source_maps(r,ra);cache.audit.r=ra;
cache.audit.B=Ba(:,active);cache.audit.D=Da(:,active);cache.audit.D2=D2a(:,active);
cache.auditCoreRadius=min(2,r(end));
wr=[r(2)-r(1);r(3:end)-r(1:end-2);r(end)-r(end-1)]/2;
wt=[theta(2)-theta(1),theta(3:end)-theta(1:end-2),theta(end)-theta(end-1)]/2;
interior=false(nr,nt);interior(2:end-1,2:end-1)=true;
cache.green=struct('area',(r.*wr)*wt,'interior',interior);
end

function maps=source_maps(r,x)
% Sparse cubic Hermite value/derivative maps from nodal values and slopes.
n=numel(r);x=x(:);cellIndex=discretize(x,r);h=r(cellIndex+1)-r(cellIndex);
u=(x-r(cellIndex))./h;rows=(1:numel(x))';
row=[rows;rows];column=[cellIndex;cellIndex+1];
make=@(a,b)sparse(row,column,[a;b],numel(x),n);
maps=struct('value',make(2*u.^3-3*u.^2+1,-2*u.^3+3*u.^2), ...
    'slope',make(h.*(u.^3-2*u.^2+u),h.*(u.^3-u.^2)), ...
    'first',make((6*u.^2-6*u)./h,(-6*u.^2+6*u)./h), ...
    'firstSlope',make(3*u.^2-4*u+1,3*u.^2-2*u), ...
    'second',make((12*u-6)./h.^2,(-12*u+6)./h.^2), ...
    'secondSlope',make((6*u-4)./h,(6*u-2)./h));
end

function [B,D,D2]=basis_matrices(knots,x)
% Cox-de Boor recurrence and exact first/second derivatives. Sparse local
% supports avoid a toolbox dependency and dense geometry maps.
x=x(:);breaks=unique(knots);span=discretize(x,breaks)+3;
assert(all(isfinite(span)),'polar:BSplineSites','Spline sites must remain inside the physical radial domain.');
count=numel(x);nt=numel(knots);
B=sparse((1:count)',span,ones(count,1),count,nt-1);D=spalloc(count,nt-1,1);D2=D;
X=spdiags(x,0,count,count);
for degree=1:3
    n=nt-degree-1;
    dl=knots(1+degree:n+degree)-knots(1:n);
    dr=knots(2+degree:n+degree+1)-knots(2:n+1);
    il=zeros(1,n);ir=il;il(dl>0)=1./dl(dl>0);ir(dr>0)=1./dr(dr>0);
    L=spdiags(il',0,n,n);R=spdiags(ir',0,n,n);
    left=B(:,1:n);right=B(:,2:n+1);
    next=(X*left-left*spdiags(knots(1:n)',0,n,n))*L+ ...
        (right*spdiags(knots(2+degree:n+degree+1)',0,n,n)-X*right)*R;
    nextD=degree*(left*L-right*R);
    nextD2=degree*(D(:,1:n)*L-D(:,2:n+1)*R);
    B=next;D=nextD;D2=nextD2;
end
assert(max(abs(full(sum(B,2))-1))<1e-11,'polar:BSplinePartition','The local B-spline basis must form a partition of unity.');
end

function [worst,core,near,infinity]=strong_audit(coefficients,sourceData,sourceSlopes,cache)
% Independent sites/analytic B-spline derivatives: this is a sampled strong
% PDE residual, not the assembled weak residual or an accuracy certificate.
audit=cache.audit;r=audit.r;m=cache.modes;nm=numel(m);
q=audit.value*sourceData+audit.slope*sourceSlopes;
wall=q(:,nm+(1:2));
wallFirst=audit.first*sourceData(:,nm+(1:2))+audit.firstSlope*sourceSlopes(:,nm+(1:2));
wallSecond=audit.second*sourceData(:,nm+(1:2))+audit.secondSlope*sourceSlopes(:,nm+(1:2));
radialLift=-4*wall-5*r.*wallFirst-r.^2.*wallSecond;
lapLift=radialLift(:,1)*cache.A+radialLift(:,2)*cache.B- ...
    wall(:,1)*cache.App-wall(:,2)*cache.Bpp;
phi=audit.B*coefficients;phiR=audit.D*coefficients;phiRR=audit.D2*coefficients;
lapModes=phiRR+phiR./r-(phi.*(m.^2))./r.^2;
reconstructed=sine_synthesis(q(:,1:nm));
reconstructed(:,1)=wall(:,1);reconstructed(:,end)=wall(:,2);
defect=-sine_synthesis(lapModes)-lapLift-reconstructed;
infinity=max(abs(defect),[],'all');
worst=infinity/max(1,max(abs(reconstructed),[],'all'));
inside=r<=cache.auditCoreRadius;
if any(inside),core=max(abs(defect(inside,:)),[],'all')/max(1,max(abs(reconstructed(inside,:)),[],'all'));else,core=NaN;end
first=1:min(4,numel(r));near=max(abs(defect(first,:)),[],'all')/max(1,max(abs(reconstructed(first,:)),[],'all'));
end

function [x,w]=gauss8()
x=[-.9602898564975363;-.7966664774136267;-.5255324099163290;-.1834346424956498; ...
    .1834346424956498;.5255324099163290;.7966664774136267;.9602898564975363];
w=[.1012285362903763;.2223810344533745;.3137066458778873;.3626837833783620; ...
    .3626837833783620;.3137066458778873;.2223810344533745;.1012285362903763];
end

function coefficients=sine_analysis(values)
% Identical DST-I normalization: c=2/(N+1)*values*S.
n=size(values,2)+1;
extension=[zeros(size(values,1),1),values,zeros(size(values,1),1),-values(:,end:-1:1)];
transformed=fft(extension,[],2);
if isreal(values),coefficients=-imag(transformed(:,2:n))/n;
else,coefficients=(transformed(:,end:-1:n+2)-transformed(:,2:n))/(2i*n);end
end
function values=sine_synthesis(coefficients)
n=size(coefficients,2)+1;
extension=[zeros(size(coefficients,1),1),coefficients,zeros(size(coefficients,1),1),-coefficients(:,end:-1:1)];
transformed=fft(extension,[],2);
if isreal(coefficients),interior=-imag(transformed(:,2:n))/2;
else,interior=(transformed(:,end:-1:n+2)-transformed(:,2:n))/(4i);end
values=[zeros(size(coefficients,1),1),interior,zeros(size(coefficients,1),1)];
end
function values=cosine_synthesis(coefficients)
n=size(coefficients,2)+1;
extension=[zeros(size(coefficients,1),1),coefficients,zeros(size(coefficients,1),1),coefficients(:,end:-1:1)];
transformed=fft(extension,[],2);values=transformed(:,1:n+1)/2;
if isreal(coefficients),values=real(values);end
end
function report=transform_test()
rng(819);sizes=[17,65,129,1025];records=struct([]);
for nt=sizes
    theta=linspace(0,pi,nt);modes=1:nt-2;S=sin(theta(2:end-1)'*modes);
    v=randn(3,nt-2)+1i*randn(3,nt-2);c=randn(3,nt-2)+1i*randn(3,nt-2);a=sine_analysis(v);
    raw=[max(abs(a-v*S*(2/(nt-1))),[],'all'), ...
        max(abs(sine_synthesis(c)-c*sin(modes'*theta)),[],'all'), ...
        max(abs(cosine_synthesis(c.*modes)-c*(modes'.*cos(modes'*theta))),[],'all')];
    scaled=raw./[max(1,max(abs(v*S*(2/(nt-1))),[],'all')), ...
        max(1,max(abs(c*sin(modes'*theta)),[],'all')), ...
        max(1,max(abs(c*(modes'.*cos(modes'*theta))),[],'all'))];
    assert(all(scaled<2e-12),'polar:FFTTransformIdentity');
    row=struct('ntheta',nt,'absoluteErrors',raw,'scaledErrors',scaled, ...
        'analysisRoundtrip',max(abs(sine_synthesis(a)-[zeros(3,1),v,zeros(3,1)]),[],'all'));
    assert(row.analysisRoundtrip<1e-13,'polar:FFTRoundtrip');
    records=[records,row]; %#ok<AGROW>
end
report=struct('kind','FFT_DST_I_and_cosine_collocation_identity_zero_PDE_test_v1', ...
    'passed',true,'records',records,'physicalPdeCalls',0,'sameSineCollocation',true);
end

function [outer,info]=green_outer(q,g,options,geometry)
% Same signed half-plane image kernel as original greenBoundary.m:105-142.
% Only the geometry/quadrature and bin coordinates are changed.
tol=opt(options,'greenSourceTolerance',1e-10);cap=opt(options,'greenMaxSources',2048);
r=g.r(:);t=g.theta(:)';nr=numel(r);nt=numel(t);
area=geometry.area;interior=geometry.interior;
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
