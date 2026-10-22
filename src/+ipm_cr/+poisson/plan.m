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


function key=geometry_key(g)
key=struct('x',g.core.x,'y',g.core.y,'r',g.polar.r,'theta',g.polar.theta, ...
    'activeCoreRows',g.core.activeRows,'activeCoreColumns',g.core.activeColumns, ...
    'coreX',g.core.X,'coreY',g.core.Y,'polarX',g.polar.X,'polarY',g.polar.Y);
end
