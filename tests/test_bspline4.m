function report=test_bspline4(output,n,width)
%TEST_BSPLINE4 Manufactured two-patch solve and spline transport, at fixed time.
if nargin<1,output=tempname;end
if nargin<2,n=[40,80,120];end
if nargin<3,width=7;end
main('config');
errors=zeros(numel(n),9);coupled=zeros(numel(n),3);data=cell(numel(n),1);
for j=1:numel(n)
    k=n(j);h=1.8/k;
    r=exp(linspace(log(.6),log(5),3*k+1)).';
    axes=struct('x',(-k:k)*h,'y',(0:k)'*h,'r',r, ...
        'theta',linspace(0,pi,2*k+1),'dxMetric',h*ones(1,2*k+1), ...
        'dyMetric',h*ones(k+1,1),'drMetric',r*log(5/.6)/(3*k));
    g=ipm_cr.mesh.build(struct('axes',axes,'staircase',true,'interpolationWidth',width));
    solver=ipm_cr.space.poisson(g);
    axes=struct('x',solver.c.xAxis,'y',solver.c.yAxis, ...
        't',solver.p.tAxis,'theta',solver.p.thetaAxis);
    assert(isequal(solver.c.X,g.core.X)&&isequal(solver.c.Y,g.core.Y)&& ...
        isequal(solver.p.X,g.polar.X)&&isequal(solver.p.Y,g.polar.Y));
    % psi = y exp(-r^2). Supply its exact forced DtN trace on r=R.
    R=g.polar.r(end);
    imposed=struct('kind','forced_dtn','exteriorStatus','manufactured_boundary', ...
        'value',@(x,y) 2*(1-R^2)*y.*exp(-R^2));
    [v,info]=solver.solve(@source,imposed);
    [f,vx,vy,reference,q,qRate]=manufactured(g);
    if j==1,check_wall_probe(g,axes,f);end
    d=ipm_cr.space.basis.gradient(f,g,axes);
    op=ipm_cr.space.transport(f,g,axes);
    a=op(vx,vy);
    da=ipm_cr.space.basis.gradient(ipm_cr.mesh.exchange(a.B,g),g,axes);
    errors(j,:)=[poisson_error(v,g),compare(a.B,reference,g),compare(d,q,g),compare(da,qRate,g), ...
        max(abs(v.core.u2(1,:))),max(abs(v.polar.u2(:,[1,end])),[],'all')];
    [coupled(j,:),where]=coupled_error(g,axes,solver,imposed);
    data{j}=struct('buildSeconds',solver.buildSeconds,'solveSeconds',info.solveSeconds, ...
        'weakBackwardError',info.weakBackwardError,'traceAbsolute',info.traceAbsolute, ...
        'coreWeakTraceMismatch',info.coreWeakTraceMismatch, ...
        'nativeDof',info.coefficientDof,'freeDof',info.freeWeakDof, ...
        'traceQuadratureCount',numel(solver.c.traceX),'angularTraceDof',numel(solver.p.jj), ...
        'cholesky',info.linearSolverStorage);
    data{j}.worstCoupledTransportDx=where;
    polar=g.polar;i=polar.active;
    x=polar.X(i);y=polar.Y(i);e=exp(-x.^2-y.^2);
    defect=[abs(v.polar.u1(i)-(2*y.^2-1).*e),abs(v.polar.u2(i)+2*x.*y.*e)];
    [~,index]=max(defect(:));[row,component]=ind2sub(size(defect),index);
    data{j}.worstPolarVelocity=struct('r',hypot(x(row),y(row)), ...
        'theta',atan2(y(row),x(row)),'component',component);
    fprintf('n=%d errors %s setup %.3gs\n',k,sprintf(' %.4g',errors(j,:)),solver.buildSeconds);
end
refinement=log(n(2:end)./n(1:end-1)).';
report=struct('intervals',n,'columns',{{'corePsi','coreVelocity','polarPsi','polarVelocity', ...
    'transport','densityDx','transportDx','coreWallNormal','polarWallNormal'}}, ...
    'errors',errors,'orders',log(errors(1:end-1,1:7)./errors(2:end,1:7))./refinement, ...
    'meshes',{data},'transferWidth',width,'timeAdvanced',0);
report.coupledColumns={'velocity','transport','transportDx'};
report.coupledErrors=coupled;
report.coupledOrders=log(coupled(1:end-1,:)./coupled(2:end,:))./refinement;
report.minimumFinePairOrder=3.8;
report.componentFinePairPass=report.orders(end,:)>report.minimumFinePairOrder;
report.coupledPass=all(report.coupledOrders(end,:)>report.minimumFinePairOrder);
report.passed=all(errors(:,8:9)<5e-12,'all')&&report.coupledPass&&all(report.componentFinePairPass);
if ~isfolder(output),mkdir(output);end
ipm_cr.io.identity('json',fullfile(output,'bspline4.json'),report);
assert(report.passed,'ipm_cr:SplineManufacturedOrder', ...
    'The coupled or individual fine-pair order missed its target.');
end

function check_wall_probe(g,axes,f)
wallX=linspace(-2.5,2.5,141);
[~,map]=ipm_cr.mesh.interpolate.sample(f,g,wallX,zeros(size(wallX)));
wall=ipm_cr.mesh.interpolate.wallPlan(g,map);
assert(any([wall.core.targets;wall.polar.targets]),'ipm_cr:WallProbeGhostCoverage');
[full,local]=ipm_cr.space.transport(f,g,axes,wall.nodes);
vx=struct();vy=struct();localX=struct();localY=struct();
for patch={'core','polar'}
    key=patch{1};p=g.(key);idx=wall.nodes.(key);
    vx.(key)=sin(p.X+.23);vy.(key)=.31*cos(p.X-p.Y);
    localX.(key)=vx.(key)(idx);localY.(key)=vy.(key)(idx);
end
wide=full(vx,vy,struct('cfl',false,'partials',true));
[b,partials]=local(localX,localY);
reference={wide.B,wide.partials.c_l,wide.partials.delta};
compact={b,partials.c_l,partials.delta};
for k=1:numel(reference)
    old=ipm_cr.mesh.interpolate.sample(ipm_cr.mesh.exchange(reference{k},g),map);
    new=ipm_cr.mesh.interpolate.sample(ipm_cr.mesh.interpolate.wallExchange(compact{k},wall),wall.map);
    assert(isequal(old,new),'ipm_cr:WallProbeEquivalence');
end
end

function [errors,where]=coupled_error(g,axes,solver,boundary)
% Feed the spline derivative of a manufactured density into Poisson,
% then use the resulting numerical velocity in the spline transport.
f=struct();b=f;db=f;
for label={'core','polar'}
    k=label{1};x=g.(k).X;y=g.(k).Y;r2=x.^2+y.^2;e=exp(-r2);
    H=sqrt(pi)/2*erf(x);
    a=2*x.*exp(-x.^2)+(6-4*y.^2).*H;
    f.(k)=y.*exp(-y.^2).*a;
    fy=exp(-y.^2).*((1-2*y.^2).*a-8*y.^2.*H);
    q=y.*(8-4*r2).*e;qx=x.*y.*(8*r2-24).*e;
    qy=e.*((8-4*r2).*(1-2*y.^2)-8*y.^2);
    u=(2*y.^2-1).*e;v=-2*x.*y.*e;
    b.(k)=-u.*q-v.*fy;
    db.(k)=2*x.*u.*q-u.*qx-(-2+4*x.^2).*y.*e.*fy-v.*qy;
end
f=ipm_cr.mesh.exchange(f,g);
d=ipm_cr.space.basis.gradient(f,g,axes);
q=struct('core',d.core.x,'polar',d.polar.x);
flow=solver.solve(q,boundary);
op=ipm_cr.space.transport(f,g,axes);
rhs=op(struct('core',flow.core.u1,'polar',flow.polar.u1), ...
    struct('core',flow.core.u2,'polar',flow.polar.u2));
grad=ipm_cr.space.basis.gradient(ipm_cr.mesh.exchange(rhs.B,g),g,axes);
err=poisson_error(flow,g);
[qError,where]=compare(grad,db,g);
errors=[max(err([2,4])),compare(rhs.B,b,g),qError];
end

function s=source(x,y)
r2=x.^2+y.^2;s=y.*(8-4*r2).*exp(-r2);
end

function e=poisson_error(v,g)
e=zeros(1,4);j=0;
for label={'core','polar'}
    key=label{1};p=g.(key);i=p.active;u=v.(key);
    x=p.X(i);y=p.Y(i);b=exp(-x.^2-y.^2);
    e(j+1)=max(abs(u.psi(i)-y.*b));
    e(j+2)=max(abs([u.u1(i)-(2*y.^2-1).*b;u.u2(i)+2*x.*y.*b]));
    j=j+2;
end
end

function [f,vx,vy,reference,q,qRate]=manufactured(g)
f=struct();vx=f;vy=f;reference=f;q=f;qRate=f;
for label={'core','polar'}
    key=label{1};p=g.(key);x=p.X;y=p.Y;
    f.(key)=exp(-(.4*(x-.7).^2+.6*(y-.8).^2));
    fx=-.8*(x-.7).*f.(key);fy=-1.2*(y-.8).*f.(key);
    vx.(key)=.6+.3*sin(x);vy.(key)=.2*sin(y);
    reference.(key)=-vx.(key).*fx-vy.(key).*fy;q.(key)=fx;
    fxx=-.8*f.(key)-.8*(x-.7).*fx;fxy=-1.2*(y-.8).*fx;
    qRate.(key)=-.3*cos(x).*fx-vx.(key).*fxx-vy.(key).*fxy;
end
f=ipm_cr.mesh.exchange(f,g);
end

function [e,where]=compare(a,b,g)
e=-inf;where=struct();
for label={'core','polar'}
    k=label{1};i=g.(k).active;
    if isstruct(a.(k)),values=a.(k).x;else,values=a.(k);end
    ids=find(i);[local,j]=max(abs(values(ids)-b.(k)(ids)));
    if local>e
        e=local;j=ids(j);
        where=struct('chart',k,'x',g.(k).X(j),'y',g.(k).Y(j));
    end
end
end
