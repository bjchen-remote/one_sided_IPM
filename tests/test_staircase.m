function report=test_staircase(outputDir,staticOnly)
%TEST_STAIRCASE Independent geometry, accuracy, transfer and driver checks.
if nargin<1,outputDir=tempname;end
if nargin<2,staticOnly=false;end
c=main('config');c.method.geometry.interpolationWidth=7;
assert(c.method.geometry.staircase&&strcmp(c.method.geometry.preset,'medium'));
angleCounts=struct('coarse',91,'medium',181,'fine',361);
for name={'coarse','medium','fine'}
    c.method.geometry.preset=name{1};
    [~,method]=ipm_cr.config(c,'run');
    g=ipm_cr.mesh.build(method.geometry);
    a=g.core;stair=a.stair;
    layers=stair.activeRadialLayers;
    assert(numel(g.polar.theta)==angleCounts.(name{1}));
    assert(g.polar.theta(1)==0&&g.polar.theta(end)==pi);
    if strcmp(name{1},'medium'),assert(isequal(coarse.polar.theta,g.polar.theta(1:2:end)));end
    assert_thin_overlap(layers);
    assert_radial_match(g);
    assert_core_spacing(g);
    assert(nnz(a.active)<numel(a.activeRows)*numel(a.activeColumns));
    assert(all(stair.support(a.active)));
    e=g.exchange;
    assert(g.interpolationWidth==7&&e.maximumSupportNodes==49);
    assert(max(sum(spones(e.coreFromPolar),2))==49&&max(sum(spones(e.polarFromCore),2))==49);
    allocation=g.axisSource.coreAllocation;
    reference=g.axisSource.interfaceCore.referenceAxes;
    budgets=[(numel(reference.x)-1)/2,numel(reference.y)-1];
    cut=allocation.coreBounds;
    counts=g.axisSource.interfaceCore.coreIntervals;
    assert(isequal(counts,[nnz(reference.x>0&reference.x<=cut(1)), ...
        nnz(reference.y>0&reference.y<=cut(2))]));
    values=struct('core',ones(size(a.X)),'polar',ones(size(g.polar.X)));
    values.core(a.ghost)=NaN;values.polar(g.polar.ghost)=NaN;
    filled=ipm_cr.mesh.exchange(values,g);
    assert(all(filled.core==1,'all')&&all(filled.polar==1,'all'));
    ghostPolynomialError=assert_ghost_polynomials(g);
    report.geometry.(name{1})=struct('sizeCore',size(a.X),'sizePolar',size(g.polar.X), ...
        'activeCore',nnz(a.active),'layers',[min(layers),median(layers),max(layers)], ...
        'observedRelativeCoreShare',prod(counts./budgets)/allocation.referenceAxisShare^2, ...
        'coreFractionBounds',cut,'interfaceRadial',g.axisSource.interfaceRadial, ...
        'interfaceCore',g.geometryKey.staircase.interfaceCore, ...
        'ghostPolynomialError',ghostPolynomialError, ...
        'mapNorms',[g.exchange.coreMapInfinityNorm,g.exchange.polarMapInfinityNorm]);
    if strcmp(name{1},'coarse'),coarse=g;end
end

solver=ipm_cr.space.poisson(coarse);
assert(isequal(solver.c.X,coarse.core.X)&&isequal(solver.p.X,coarse.polar.X));
[u,info]=solver.solve(@source);
for name={'core','polar'}
    k=name{1};m=coarse.(k);x=m.X(m.active);y=m.Y(m.active);
    e=exp(-(x.^2+y.^2));
    defect=[u.(k).psi(m.active)-y.*e;u.(k).u1(m.active)+(1-2*y.^2).*e; ...
        u.(k).u2(m.active)+2*x.*y.*e];
    report.manufactured.(k)=max(abs(defect));
end
assert(max([report.manufactured.core,report.manufactured.polar])<.005&&info.traceAbsolute<1e-8);
report.manufacturedTraceError=info.traceAbsolute;
c.method.geometry.preset='coarse';
if staticOnly
    [s,ctx]=ipm_cr.init(c.method);
    report.timeAdvanced=0;
    if ~isfolder(outputDir),mkdir(outputDir);end
else
    c.run.outputDir=outputDir;
    c.run.verbose=false;c.run.tauEnd=.001;c.run.wallSeconds=120;
    [summary,s,ctx]=main(c);
    assert(summary.targetReached&&isfield(ctx.mesh.core,'stair'));
    c.run.tauEnd=.002;
    [summary,s,ctx]=main(c);
    assert(summary.targetReached&&summary.resumed&&s.scale.canonicalTime==.002);
    report.driver=struct('tau',s.scale.canonicalTime,'step',s.step,'resumed',summary.resumed, ...
        'QInfinity',max([abs(s.rhs.Q.core(:));abs(s.rhs.Q.polar(:))]));
end
assert(isfinite(s.rhs.rate)&&all(isfinite(s.rhs.Q.core),'all'));

initial=ipm_cr.mesh.axes.reference(ctx.mesh);
memory=struct('initialMesh',initial,'initialExtents',[1,1]);
ax=ipm_cr.mesh.axes.propose([1.12,1.06],memory);
settings=struct('candidateAxes',ax);
[moved,new,change]=ipm_cr.remesh(s,ctx,settings,true);
assert(change.remeshed&&isfield(new.mesh.core,'stair')&&~change.transfer.ghostValuesRead);
assert(change.transfer.polynomialDegree==5&&change.transfer.interpolationWidth==6);
assert(new.mesh.interpolationWidth==7);
assert(isequal(new.mesh.polar.theta,coarse.polar.theta));
assert(isequaln(moved.scale,s.scale)&&isequaln(moved.originalAuthority,s.originalAuthority)&& ...
    moved.step==s.step&&all(isfinite(moved.rhs.Q.core),'all'));
layers=new.mesh.core.stair.activeRadialLayers;
report.remesh=struct('layers',[min(layers),median(layers),max(layers)], ...
    'sameScaleAndClock',true,'QInfinity',max([abs(moved.rhs.Q.core(:));abs(moved.rhs.Q.polar(:))]));
assert_thin_overlap(layers);
assert_radial_match(new.mesh);
assert_core_spacing(new.mesh);
report.passed=true;
ipm_cr.io.identity('json',fullfile(outputDir,'staircase_test.json'),report);
end

function result=assert_ghost_polynomials(g)
c=g.core;p=g.polar;
cartesian=@(x,y) (x/c.x(end)).^4.*(1+(y/c.y(end)).^4);
% This polynomial checks the closest active radial donors without using
% the 1e10 far endpoint as the scale for the inner-ring data.
radialScale=p.r(g.ghostWidth+6);
polar=@(r,t) (r/radialScale).^4.*(1+(t/pi).^4);
f=struct('core',cartesian(c.X,c.Y),'polar',polar(p.r,p.theta));
f.core(c.ghost)=NaN;f.polar(p.ghost)=NaN;
filled=ipm_cr.mesh.exchange(f,g);
x=c.X(c.ghost);y=c.Y(c.ghost);
reference=polar(hypot(x,y),atan2(y,x));
coreError=max(abs(filled.core(c.ghost)-reference)./max(1,abs(reference)));
reference=cartesian(p.X(p.ghost),p.Y(p.ghost));
polarError=max(abs(filled.polar(p.ghost)-reference)./max(1,abs(reference)));
result=[coreError,polarError];
assert(all(result<2e-11),'ipm_cr:StairGhostPolynomial');
end

function assert_core_spacing(g)
spacing=g.geometryKey.staircase.interfaceCore;
limit=g.axisSource.interfaceRadial.adjacentStepRatioLimit;
axes={g.core.x(g.core.originIndex:end),g.core.y};
metrics={g.core.dxMetric(g.core.originIndex:end),g.core.dyMetric};
for direction=1:2
    axis=axes{direction};metric=metrics{direction};
    if direction==1
        band=spacing.x;
        outer=g.core.activeXBounds(2);
    else
        band=spacing.y;
        outer=g.core.activeYBounds(2);
    end
    step=diff(axis(:));lo=band.seamCells(1);hi=band.seamCells(2);
    assert(axis(hi+1)>=outer&&axis(lo)<outer);
    assert(find(axis==outer,1)+g.ghostWidth<=hi+1);
    assert(max(abs(step(1:band.centralCells)/band.centralStep-1))<1e-11);
    assert(max(abs(step(lo:hi)/band.seamStep-1))<1e-11);
    assert(all(metric(lo:hi+1)==metric(lo)));
    adjacent=step(2:end)./step(1:end-1);
    assert(max([adjacent;1./adjacent])<limit);
end
end

function assert_radial_match(g)
radial=g.geometryKey.staircase.interfaceRadial;
rows=1:g.ghostWidth+1+radial.uniformActiveCells;
h=radial.radialStep;
assert(max(abs(diff(g.polar.r(rows))/h-1))<2e-12);
assert(all(g.polar.drMetric(rows)==h));
assert(h<=radial.baseRadialStep&&h<=radial.medianCoreRadialStep);
limit=1.1;
if numel(g.polar.theta)==91,limit=1.2;end
last=find(g.polar.r==radial.transitionEndRadius,1);
step=diff(g.polar.r(1:min(last+1,end)));
ratio=step(2:end)./step(1:end-1);
actual=max([ratio;1./ratio]);
assert(radial.adjacentStepRatioLimit==limit&&actual<limit);
assert(abs(actual-radial.maxAdjacentStepRatio)<2e-12);
fullSteps=diff(g.polar.r);
fullRatio=fullSteps(2:end)./fullSteps(1:end-1);
assert(max([fullRatio;1./fullRatio])<limit);
expected=min(radial.baseRadialStep,median(radial.coreStepAtAngle));
assert(abs(h/expected-1)<2e-12);
assert(isequal(radial.radialToCoreStepRange, ...
    [min(h./radial.coreStepAtAngle),max(h./radial.coreStepAtAngle)]));
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
