function report=test_staircase(outputDir)
%TEST_STAIRCASE Independent geometry, accuracy, transfer and driver checks.
if nargin<1,outputDir=tempname;end
c=main('config');
assert(c.method.geometry.staircase&&strcmp(c.method.geometry.preset,'medium'));
angleCounts=struct('coarse',91,'medium',181,'fine',361);
for name={'coarse','medium','fine'}
    c.method.geometry.preset=name{1};
    g=ipm_cr.mesh.build(c.method.geometry);
    a=g.core;stair=a.stair;
    layers=stair.activeRadialLayers;
    assert(numel(g.polar.theta)==angleCounts.(name{1}));
    assert(g.polar.theta(1)==0&&g.polar.theta(end)==pi);
    if strcmp(name{1},'medium'),assert(isequal(coarse.polar.theta,g.polar.theta(1:2:end)));end
    assert_thin_overlap(layers);
    assert(nnz(a.active)<numel(a.activeRows)*numel(a.activeColumns));
    assert(all(stair.support(a.active)));
    allocation=g.axisSource.coreAllocation;
    budgets=[(a.nx-1)/2,a.ny-1];
    cut=allocation.coreBounds;
    counts=[nnz(a.x>0&a.x<=cut(1)),nnz(a.y>0&a.y<=cut(2))];
    assert(all(abs(counts-allocation.targetAxisShare.*budgets)<=1));
    values=struct('core',ones(size(a.X)),'polar',ones(size(g.polar.X)));
    values.core(a.ghost)=NaN;values.polar(g.polar.ghost)=NaN;
    filled=ipm_cr.mesh.fill(values,g);
    assert(all(filled.core==1,'all')&&all(filled.polar==1,'all'));
    linear=struct('core',a.X+.3*a.Y,'polar',g.polar.X+.3*g.polar.Y);
    linear.core(a.ghost)=NaN;
    d=ipm_cr.mesh.diff(linear,g);
    assert(max(abs(d.core.x(a.active)-1))<2e-12&&max(abs(d.core.y(a.active)-.3))<2e-12);
    report.geometry.(name{1})=struct('sizeCore',size(a.X),'sizePolar',size(g.polar.X), ...
        'activeCore',nnz(a.active),'layers',[min(layers),median(layers),max(layers)], ...
        'observedRelativeCoreShare',prod(counts./budgets)/allocation.referenceAxisShare^2, ...
        'coreFractionBounds',cut,'mapNorms',[g.exchange.coreMapInfinityNorm,g.exchange.polarMapInfinityNorm]);
    if strcmp(name{1},'coarse'),coarse=g;end
end

solver=ipm_cr.poisson.fem.solver(coarse);
assert(isequal(solver.c.X,coarse.core.X)&&isequal(solver.p.X,coarse.polar.X));
boundary=struct('kind','green_dirichlet','value',@(x,y)0*x+0*y, ...
    'exteriorStatus','declared_compact_source');
[u,info]=solver.solve(@source,boundary);
for name={'core','polar'}
    k=name{1};m=coarse.(k);x=m.X(m.active);y=m.Y(m.active);
    e=exp(-(x.^2+y.^2));
    defect=[u.(k).psi(m.active)-y.*e;u.(k).u1(m.active)+(1-2*y.^2).*e; ...
        u.(k).u2(m.active)+2*x.*y.*e];
    report.manufactured.(k)=max(abs(defect));
end
assert(max([report.manufactured.core,report.manufactured.polar])<.005&&info.traceAbsolute<1e-8);
report.manufacturedTraceError=info.traceAbsolute;

c.method.geometry.preset='coarse';c.run.outputDir=outputDir;
c.run.verbose=false;c.run.tauEnd=.001;c.run.wallSeconds=120;
[summary,s,ctx]=main(c);
assert(summary.targetReached&&isfield(ctx.mesh.core,'stair'));
c.run.tauEnd=.002;
[summary,s,ctx]=main(c);
assert(summary.targetReached&&summary.resumed&&s.scale.canonicalTime==.002);
assert(isfinite(s.rhs.rate)&&all(isfinite(s.rhs.Q.core),'all'));
report.driver=struct('tau',s.scale.canonicalTime,'step',s.step,'resumed',summary.resumed, ...
    'QInfinity',max([abs(s.rhs.Q.core(:));abs(s.rhs.Q.polar(:))]));

raw=jsondecode(fileread(fullfile(ipm_cr.root(),'data','mesh.json')));
old=ctx.mesh;
ax=struct('x',1.12*old.core.x,'y',1.06*old.core.y,'theta',old.polar.theta, ...
    'r',raw.axes.r(1:2:end),'drMetric',2*raw.axes.drMetric(1:2:end), ...
    'dxMetric',1.12*old.core.dxMetric,'dyMetric',1.06*old.core.dyMetric);
settings=struct('candidateAxes',ax);
[moved,new,change]=ipm_cr.remesh(s,ctx,settings,true);
assert(change.remeshed&&isfield(new.mesh.core,'stair')&&~change.transfer.ghostValuesRead);
assert(isequal(new.mesh.polar.theta,coarse.polar.theta));
assert(isequaln(moved.scale,s.scale)&&isequaln(moved.originalAuthority,s.originalAuthority)&& ...
    moved.step==s.step&&all(isfinite(moved.rhs.Q.core),'all'));
layers=new.mesh.core.stair.activeRadialLayers;
report.remesh=struct('layers',[min(layers),median(layers),max(layers)], ...
    'sameScaleAndClock',true,'QInfinity',max([abs(moved.rhs.Q.core(:));abs(moved.rhs.Q.polar(:))]));
assert_thin_overlap(layers);
report.passed=true;
ipm_cr.io.json(fullfile(outputDir,'staircase_test.json'),report);
end

function assert_thin_overlap(layers)
% A minority of rays cross an extra staircase corner near the wall.
assert(ismember(median(layers),[1,2])&&min(layers)>=1&&max(layers)<=4);
assert(mean(layers<=2)>=.7);
end

function q=source(x,y)
r2=x.^2+y.^2;
q=(8-4*r2).*y.*exp(-r2);
end
