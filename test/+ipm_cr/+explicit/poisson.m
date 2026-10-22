function [flow,info]=poisson(source,g,solverOptions,analyticSource,includePsi)
%POISSON Cached spectral hybrid solve for uniform-core/annular transport data.
% SOURCE is the complete already differentiated physical source on the two
% transport grids. Its unique active-core owner supplies each spectral
% query, including spectral annular points lying in the Cartesian owner.
% Return psi and both Cartesian derivatives of the SAME spectral solution
% on BOTH full transport grids. This function performs no transport step.
if nargin<3,solverOptions=struct();end
if nargin<4,analyticSource=[];end
if nargin<5,includePsi=true;end
validateattributes(includePsi,{'logical'},{'scalar'});
assert(isempty(analyticSource)||isa(analyticSource,'function_handle'), ...
    'ipm_cr:ExplicitAnalyticSource','An analytic source offset must be a function handle.');
timer=tic;
assert(isstruct(g)&&isscalar(g)&&all(isfield(g,{'core','polar'})), ...
    'ipm_cr:ExplicitPoissonGrid','Supply the uniform Cartesian/annular transport geometry.');
assert(isstruct(source)&&isscalar(source)&&all(isfield(source,{'core','polar'}))&& ...
    isnumeric(source.core)&&isreal(source.core)&&isequal(size(source.core),size(g.core.X))&& ...
    isnumeric(source.polar)&&isreal(source.polar)&&isequal(size(source.polar),size(g.polar.X))&& ...
    all(isfinite(source.core),'all')&&all(isfinite(source.polar),'all'), ...
    'ipm_cr:ExplicitPoissonSource','Both complete native source arrays must be finite and real.');
assert(isstruct(solverOptions)&&isscalar(solverOptions), ...
    'ipm_cr:ExplicitPoissonOptions','Use a scalar spectral Solver options structure.');
defaults=struct('coreDegree',[48,48],'angularDegree',64,'polarDegree',16, ...
    'rInner',.75,'rOuter',g.polar.r(end));
names=fieldnames(solverOptions);
for k=1:numel(names),defaults.(names{k})=solverOptions.(names{k});end
cfg=ipm_cr.config(defaults);
assert(abs(cfg.rOuter-g.polar.r(end))<=64*eps(g.polar.r(end)), ...
    'ipm_cr:ExplicitPoissonCoverage','The spectral and transport exterior radii must agree.');
geometry=geometry_key(g);persistent cache
rebuild=isempty(cache)||~isequaln(cache.geometry,geometry)||~isequaln(cache.cfg,cfg);
if rebuild
    next=prepare_geometry(g,geometry,cfg);cache=next;
end
phase=tic;
qc=apply_map(cache.sourceCore,source);qp=apply_map(cache.sourcePolar,source);
if ~isempty(analyticSource)
    qc=qc+analyticSource(cache.solver.Core.X,cache.solver.Core.Y);
    qp=qp+analyticSource(cache.solver.Polar.X,cache.solver.Polar.Y);
    assert(all(isfinite(qc),'all')&&all(isfinite(qp),'all'), ...
        'ipm_cr:ExplicitAnalyticSource','The complete analytic-offset spectral source must be finite.');
end
sourceSeconds=toc(phase);
[solution,info]=cache.solver.solve(struct('core',qc,'polar',qp));
phase=tic;
core=return_core(cache,solution,g,includePsi);
polar=return_polar(cache,solution,g,includePsi);
returnSeconds=toc(phase);
flow=struct('core',core,'polar',polar,'source',source,'solution',solution);
info.definition='minus_laplacian_psi_equals_physical_source_on_unique_active_transport_owner';
info.discretizationContract='uniform_core_origin_free_annular_fifth_order_with_spectral_overlap_Poisson_v2';
info.linearSolveBackwardError=info.linearBackwardError;
info.linearSolveComponentwiseBackwardError=info.linearBackwardError;
info.relativeResidual=info.linearBackwardError;
info.algebraicResidual=info.linearBackwardError;
info.factorCacheRebuilt=rebuild;info.sourceMapCacheReused=~rebuild;
info.returnMapCacheReused=~rebuild;info.explicitGeometryBuildSeconds=cache.buildSeconds;
info.sourceInterpolationSeconds=sourceSeconds;info.returnInterpolationSeconds=returnSeconds;
info.totalAdapterSolveSeconds=sourceSeconds+info.solveSeconds+returnSeconds;
info.totalCallSeconds=toc(timer);
info.sourceReconstruction='local seven-point physical-axis tensor polynomial with unique active Cartesian ownership';
info.sourceReconstructionStencilWidth=7;info.sourceReconstructionMaximumSupportNodes=49;
info.sourceReconstructionUsesGlobalInverse=false;
info.sourceInterpolationOperatorInfinityNorm=max(cache.sourceCore.operatorInfinityNorm,cache.sourcePolar.operatorInfinityNorm);
sourceInfinity=max([max(abs(source.core),[],'all'),max(abs(source.polar),[],'all'),realmin]);
info.sourceInterpolationInfinityAmplification=max([max(abs(qc),[],'all'),max(abs(qp),[],'all')])/sourceInfinity;
info.sourceReconstructionInfinityAmplification=info.sourceInterpolationInfinityAmplification;
info.sourceActiveCoreBounds=cache.sourceCore.activeCoreBounds;
info.sourceCoreOwnerQueries=nnz(cache.sourceCore.coreOwner)+nnz(cache.sourcePolar.coreOwner);
info.sourcePolarOwnerQueries=nnz(cache.sourceCore.polarOwner)+nnz(cache.sourcePolar.polarOwner);
info.polarOriginPresent=false;info.densityInterpolatedOrDifferentiatedByPoisson=false;
info.sameSemiDiscreteOperatorClaimed=false;info.iterativeLinearSolverUsed=false;
info.velocityDerivativesConsistentWithSpectralInterpolants=true;
info.coreReturnMethod=cache.coreReturnMethod;
info.returnIncludesPsi=includePsi;
info.maximumOffGridStrongResidualScaledBySource=NaN;
info.strongResidualIsAlgebraicResidual=false;info.strongResidualIsAccuracyCertificate=false;
info.sourceInterpolationAccuracyCertified=false;info.outerSourceTailCertified=false;
info.analyticBackgroundSource=~isempty(analyticSource);
if ~isempty(analyticSource)
    info.definition='minus_laplacian_psi_equals_analytic_background_x_plus_reconstructed_residual_x';
    info.sourceReconstruction='analytic background derivative at spectral nodes plus local residual source interpolation';
    info.sourceInterpolationInfinityAmplification=NaN;
    info.sourceReconstructionInfinityAmplification=NaN;
end
info.wallNormalVelocityInfinity=max([max(abs(core.u2(g.core.Y==0))), ...
    max(abs(polar.u2(:,[1,end])),[],'all')]);
end

function c=prepare_geometry(g,geometry,cfg)
timer=tic;s=ipm_cr.Solver(cfg);sc=s.Core;p=s.Polar;
zerosNative=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
[~,sourceCore]=ipm_cr.explicit.interpolate(zerosNative,g,sc.X,sc.Y);
[~,sourcePolar]=ipm_cr.explicit.interpolate(zerosNative,g,p.X,p.Y);
assert(sourceCore.allValid&&sourcePolar.allValid, ...
    'ipm_cr:ExplicitPoissonCoverage','Every spectral source point needs a valid unique native owner.');
coreSeparable=all(abs(g.core.x)<=cfg.coreHalfWidth*(1+16*eps))&& ...
    all(g.core.y>=0&g.core.y<=cfg.coreHeight*(1+16*eps));
if coreSeparable
    returnCoreX=ipm_cr.internal.bary(sc.x,g.core.x);
    returnCoreY=ipm_cr.internal.bary(sc.y,g.core.y);
    coreReturnMethod='separable exact spectral interpolation onto uniform Cartesian transport tensor';
else
    returnCoreX=[];returnCoreY=[];
    coreReturnMethod='generic exact spectral queries onto complete Cartesian transport geometry';
end
X=g.polar.X;Y=g.polar.Y;r=g.polar.r(:);theta=g.polar.theta(:)';
owner=abs(X)<=cfg.coreHalfWidth&Y>=0&Y<=cfg.coreHeight;
coreIndex=find(owner);
polarCoreX=ipm_cr.internal.bary(sc.x,X(coreIndex));
polarCoreY=ipm_cr.internal.bary(sc.y,Y(coreIndex));
returnAngular=ipm_cr.internal.bary(p.theta,theta);
returnAngularDerivative=returnAngular*p.Dtheta;
t=log(max(p.rInner,min(p.rOuter,r)));
eligible=r>=p.rInner*(1-16*eps)&r<=p.rOuter*(1+16*eps);
panels=cell(p.nPanels,1);covered=false(size(r));
for k=1:p.nPanels
    if k==p.nPanels,rows=find(eligible&t>=p.logEdges(k)&t<=p.logEdges(k+1));
    else,rows=find(eligible&t>=p.logEdges(k)&t<p.logEdges(k+1));end
    covered(rows)=true;
    panels{k}=struct('rows',rows,'radial',ipm_cr.internal.bary(p.panelNodes{k},t(rows)));
end
assert(all(covered(any(~owner,2))), ...
    'ipm_cr:ExplicitPoissonReturnCoverage','Every transport query outside the spectral core needs an annular panel.');
c=struct('geometry',geometry,'cfg',cfg,'solver',s,'sourceCore',sourceCore,'sourcePolar',sourcePolar, ...
    'coreSeparable',coreSeparable,'returnCoreX',returnCoreX,'returnCoreY',returnCoreY, ...
    'coreReturnMethod',coreReturnMethod,'polarCoreIndex',coreIndex, ...
    'polarCoreX',polarCoreX,'polarCoreY',polarCoreY, ...
    'returnAngular',returnAngular,'returnAngularDerivative',returnAngularDerivative, ...
    'returnPanels',{panels},'buildSeconds',toc(timer));
end

function values=apply_map(map,source)
values=reshape(map.core*source.core(:)+map.polar*source.polar(:),map.shape);
assert(all(isfinite(values),'all'),'ipm_cr:ExplicitPoissonSourceMap', ...
    'A complete spectral source interpolation became nonfinite.');
end

function value=return_core(c,solution,g,includePsi)
sc=c.solver.Core;
if c.coreSeparable
    By=c.returnCoreY;Bx=c.returnCoreX;
    psi=[];if includePsi,psi=By*solution.core.psi*Bx';end
    ux=-By*(sc.Dy*solution.core.psi)*Bx';
    uy=By*(solution.core.psi*sc.Dx')*Bx';
else
    v=c.solver.evaluate(solution,g.core.X,g.core.Y);psi=v.psi;ux=v.u1;uy=v.u2;
end
wall=g.core.Y==0;if includePsi,psi(wall)=0;end;uy(wall)=0;
assert(all(isfinite(psi),'all')&&all(isfinite(ux),'all')&&all(isfinite(uy),'all'), ...
    'ipm_cr:ExplicitPoissonReturn','The complete Cartesian return must be finite.');
value=struct('psi',psi,'u1',ux,'u2',uy);
end

function value=return_polar(c,solution,g,includePsi)
p=c.solver.Polar;sc=c.solver.Core;
r=g.polar.r(:);theta=g.polar.theta(:)';shape=size(g.polar.X);
ux=nan(shape);uy=ux;psi=[];if includePsi,psi=ux;end
co=cos(theta);si=sin(theta);si([1,end])=0;
for k=1:p.nPanels
    m=c.returnPanels{k};rows=m.rows;if isempty(rows),continue;end
    U=solution.polar.psi(p.panelIndices{k},:);sampled=m.radial*U;
    if includePsi,psi(rows,:)=sampled*c.returnAngular';end
    dt=(m.radial*(p.panelD{k}*U))*c.returnAngular';
    da=sampled*c.returnAngularDerivative';
    ux(rows,:)=-(si.*dt+co.*da)./r(rows);
    uy(rows,:)=(co.*dt-si.*da)./r(rows);
end
contract=@(U)sum((c.polarCoreY*U).*c.polarCoreX,2);
if includePsi,psi(c.polarCoreIndex)=contract(solution.core.psi);end
ux(c.polarCoreIndex)=-contract(sc.Dy*solution.core.psi);
uy(c.polarCoreIndex)=contract(solution.core.psi*sc.Dx');
if includePsi,psi(:,[1,end])=0;end;uy(:,[1,end])=0;
assert(all(isfinite(psi),'all')&&all(isfinite(ux),'all')&&all(isfinite(uy),'all'), ...
    'ipm_cr:ExplicitPoissonReturn','The complete annular return must be finite.');
value=struct('psi',psi,'u1',ux,'u2',uy);
end

function key=geometry_key(g)
key=struct('x',g.core.x,'y',g.core.y,'r',g.polar.r,'theta',g.polar.theta, ...
    'activeCoreRows',g.core.activeRows,'activeCoreColumns',g.core.activeColumns, ...
    'coreX',g.core.X,'coreY',g.core.Y,'polarX',g.polar.X,'polarY',g.polar.Y);
end
