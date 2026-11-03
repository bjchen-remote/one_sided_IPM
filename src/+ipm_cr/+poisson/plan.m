function c=plan(g,solverOptions)
%PLAN Build one simulation's spectral factors and fixed transfer maps.
if nargin<2,solverOptions=struct();end
defaults=struct('coreDegree',[48,48],'angularDegree',64,'polarDegree',16, ...
    'rInner',.75,'rOuter',g.polar.r(end));
names=fieldnames(solverOptions);
for k=1:numel(names),defaults.(names{k})=solverOptions.(names{k});end
cfg=ipm_cr.pconfig(defaults);
assert(abs(cfg.rOuter-g.polar.r(end))<=64*eps(g.polar.r(end)), ...
    'ipm_cr:ExplicitPoissonCoverage','The spectral and transport exterior radii must agree.');
geometry=geometry_key(g);
timer=tic;s=ipm_cr.psolver(cfg);sc=s.Core;p=s.Polar;
zerosNative=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
[~,sourceCore]=ipm_cr.mesh.sample(zerosNative,g,sc.X,sc.Y);
[~,sourcePolar]=ipm_cr.mesh.sample(zerosNative,g,p.X,p.Y);
assert(sourceCore.allValid&&sourcePolar.allValid, ...
    'ipm_cr:ExplicitPoissonCoverage','Every spectral source point needs a valid unique native owner.');
audit=s.AuditQueries;
[~,auditCore]=ipm_cr.mesh.sample(zerosNative,g,audit.core.x,audit.core.y);
[~,auditPolar]=ipm_cr.mesh.sample(zerosNative,g,audit.polar.x,audit.polar.y);
assert(auditCore.allValid&&auditPolar.allValid, ...
    'ipm_cr:ExplicitPoissonCoverage','Every fixed native-source audit point needs a valid owner.');
coreSeparable=all(abs(g.core.x)<=cfg.coreHalfWidth*(1+16*eps))&& ...
    all(g.core.y>=0&g.core.y<=cfg.coreHeight*(1+16*eps));
if coreSeparable
    returnCoreX=ipm_cr.internal.bary(sc.x,g.core.x);
    returnCoreY=ipm_cr.internal.bary(sc.y,g.core.y);
    returnCorePlan=[];
    coreReturnMethod='separable exact spectral interpolation onto Cartesian transport tensor';
else
    returnCoreX=[];returnCoreY=[];
    returnCorePlan=core_return_plan(g,sc,p,cfg);
    coreReturnMethod='cached tensor subrectangle and exact annular panel queries onto complete Cartesian transport geometry';
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
    'auditSourceCore',auditCore,'auditSourcePolar',auditPolar, ...
    'coreSeparable',coreSeparable,'returnCoreX',returnCoreX,'returnCoreY',returnCoreY, ...
    'coreReturnMethod',coreReturnMethod,'returnCorePlan',returnCorePlan,'polarCoreIndex',coreIndex, ...
    'polarCoreX',polarCoreX,'polarCoreY',polarCoreY, ...
    'returnAngular',returnAngular,'returnAngularDerivative',returnAngularDerivative, ...
    'returnPanels',{panels},'buildSeconds',toc(timer));
end


function plan=core_return_plan(g,sc,p,cfg)
% The spectral core stays fixed when the TRANSPORT tensor grows or maps.
% Preserve tensor contractions on its interior subrectangle. Outside it,
% cache exactly the unique-owner/panel queries used by psolver.evaluate.
x=g.core.x(:)';y=g.core.y(:);nx=numel(x);ny=numel(y);
assert(isequal(g.core.X,repmat(x,ny,1))&&isequal(g.core.Y,repmat(y,1,nx)), ...
    'ipm_cr:ExplicitPoissonTensor','The Cartesian transport return must be a tensor grid.');
rows=find(y>=0&y<=cfg.coreHeight);columns=find(abs(x)<=cfg.coreHalfWidth);
owner=false(ny,nx);owner(rows,columns)=true;
X=g.core.X(:);Y=g.core.Y(:);r=hypot(X,Y);
valid=isfinite(X)&isfinite(Y)&Y>=0&r<=cfg.rOuter*(1+16*eps);
assert(all(valid),'ipm_cr:ExplicitPoissonReturnCoverage', ...
    'All Cartesian transport queries must lie in the physical half disk.');
outer=~owner(:);
assert(all(r(outer)>=p.rInner-128*eps(p.rInner)), ...
    'ipm_cr:ExplicitPoissonReturnCoverage', ...
    'Cartesian queries outside the fixed spectral core need annular support.');
theta=min(pi,max(0,atan2(max(Y,0),X)));
t=log(max(p.rInner,min(p.rOuter,r)));
panels=cell(p.nPanels,1);covered=owner(:);
for k=1:p.nPanels
    if k==p.nPanels,selected=outer&t>=p.logEdges(k)&t<=p.logEdges(k+1);
    else,selected=outer&t>=p.logEdges(k)&t<p.logEdges(k+1);end
    indices=find(selected);covered(indices)=true;
    Br=ipm_cr.internal.bary(p.panelNodes{k},t(indices));
    Ba=ipm_cr.internal.bary(p.theta,theta(indices));
    panels{k}=struct('indices',indices,'radial',Br,'angular',Ba, ...
        'angularDerivative',Ba*p.Dtheta,'radius',r(indices), ...
        'cosTheta',cos(theta(indices)),'sinTheta',sin(theta(indices)));
end
assert(all(covered),'ipm_cr:ExplicitPoissonReturnCoverage', ...
    'Every Cartesian query requires one fixed spectral owner and panel.');
plan=struct('kind','fixed_spectral_tensor_subset_annular_query_cache_v1', ...
    'coreRows',rows,'coreColumns',columns,'coreOwner',owner, ...
    'coreX',ipm_cr.internal.bary(sc.x,x(columns)), ...
    'coreY',ipm_cr.internal.bary(sc.y,y(rows)), ...
    'panels',{panels},'tensorQueryCount',nnz(owner), ...
    'annularQueryCount',nnz(outer),'allQueriesCovered',true, ...
    'outsideCoreUsesPolar',true,'polynomialExtrapolationUsed',false, ...
    'fieldValuesCached',false,'sameSolutionAnalyticDerivatives',true);
end

function key=geometry_key(g)
key=struct('x',g.core.x,'y',g.core.y,'r',g.polar.r,'theta',g.polar.theta, ...
    'activeCoreRows',g.core.activeRows,'activeCoreColumns',g.core.activeColumns, ...
    'coreX',g.core.X,'coreY',g.core.Y,'polarX',g.polar.X,'polarY',g.polar.Y);
end
