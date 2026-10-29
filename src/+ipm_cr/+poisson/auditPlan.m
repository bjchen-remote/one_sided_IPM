function plan=auditPlan(c,p,g)
%AUDITPLAN Fixed queries and spectral interpolation, containing no RHS data.
theta=pi*((0:g.auditPoints-1)+.37)/g.auditPoints;
radius=linspace(g.rInner,.98*min(g.coreHalfWidth,g.coreHeight),5)';
x=radius*cos(theta);y=radius*sin(theta);x=x(:);y=y(:);
nx=x./hypot(x,y);ny=y./hypot(x,y);
z=linspace(0,1,g.auditPoints+2)';z=z(2:end-1);
x=[x;-g.coreHalfWidth*ones(size(z));g.coreHalfWidth*ones(size(z));g.coreHalfWidth*(2*z-1)];
y=[y;g.coreHeight*z;g.coreHeight*z;g.coreHeight*ones(size(z))];
nx=[nx;-ones(size(z));ones(size(z));zeros(size(z))];
ny=[ny;zeros(size(z));zeros(size(z));ones(size(z))];
keep=hypot(x,y)>=g.rInner;x=x(keep);y=y(keep);nx=nx(keep);ny=ny(keep);
z=((0:24)+.29)/25;
coreX=g.coreHalfWidth*(2*z-1);coreY=g.coreHeight*z;
[xx,yy]=meshgrid(coreX,coreY);
coreQueries=struct('x',xx(:),'y',yy(:));
mid=(p.logEdges(1:end-1)+p.logEdges(2:end))/2;
half=diff(p.logEdges)/2;
ts=reshape(mid+[-sqrt(3/5);0;sqrt(3/5)]*half,1,[]);
[tt,aa]=ndgrid(ts,theta);
polarQueries=struct('x',exp(tt).*cos(aa),'y',exp(tt).*sin(aa));
queries=struct('core',coreQueries,'polar',polarQueries, ...
    'polarPanel',repmat(repelem((1:p.nPanels)',3),1,numel(theta)));
% Native reconstruction audits use these tensor coordinates directly. The
% generic plans below remain for the unchanged continuous-callback gate.
coreTensor=struct('Bx',ipm_cr.internal.bary(c.x,coreX), ...
    'By',ipm_cr.internal.bary(c.y,coreY));
angular=ipm_cr.internal.bary(p.theta,theta);
polarPanels=cell(p.nPanels,1);
for j=1:p.nPanels
    rows=(j-1)*3+(1:3);
    polarPanels{j}=struct('rows',rows, ...
        'Br',ipm_cr.internal.bary(p.panelNodes{j},ts(rows)),'r',exp(ts(rows))');
end
polarTensor=struct('angular',angular,'angularSecond',angular*p.Dtheta2, ...
    'panels',{polarPanels},'shape',size(polarQueries.x));
ix=p.X(1,p.interior)';iy=p.Y(1,p.interior)';
plan=struct('queries',queries,'x',x,'y',y,'nx',nx,'ny',ny, ...
    'traceInnerCore',ipm_cr.poisson.core.prepare(c,ix,iy), ...
    'traceOuterPolar',ipm_cr.poisson.polar.prepare(p,c.traceX,c.traceY,{'psi'}), ...
    'overlapCore',ipm_cr.poisson.core.prepare(c,x,y), ...
    'overlapPolar',ipm_cr.poisson.polar.prepare(p,x,y), ...
    'originCore',ipm_cr.poisson.core.prepare(c,0,0), ...
    'nativeCoreTensor',coreTensor,'nativePolarTensor',polarTensor, ...
    'nativeCore',ipm_cr.poisson.core.prepare(c,coreQueries.x,coreQueries.y), ...
    'nativePolar',ipm_cr.poisson.polar.prepare(p,polarQueries.x,polarQueries.y,{'laplacian'}));
end
