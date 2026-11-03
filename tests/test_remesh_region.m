function report=test_remesh_region(outputDir,controls)
%TEST_REMESH_REGION Pure geometry tests of the positive-q main component.
% Only mesh.build and synthetic q are used; no context, RHS or PDE solve.
if nargin<2,controls=struct();end
assert(isstruct(controls)&&isscalar(controls)&& ...
    all(ismember(fieldnames(controls),{'wallSeconds'})),'ipm_cr:RegionTestOptions');
budget=90;if isfield(controls,'wallSeconds'),budget=controls.wallSeconds;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','remesh_region');end
if ~isfolder(outputDir),mkdir(outputDir);end
file=fullfile(outputDir,'remesh_region.json');
assert(~isfile(file),'ipm_cr:RegionTestFresh','Use a fresh output directory.');
timer=tic;report=struct('kind','static_positive_q_region_geometry_v1', ...
    'passed',false,'timeAdvanced',0,'contextConstructed',false,'newRHS',false, ...
    'newPoissonSolve',false,'continuumAccuracyCertified',false);
try
    g=ipm_cr.mesh.build(geometry());p=ipm_cr.remesh('config');
    assert(p.level==.25,'ipm_cr:RegionTestDefaultLevel');
    widths=[.65,.42];k=sqrt(log(1/p.level));
    base=state(g,widths,1);shape=region(base,g,p);
    report.base=check_ellipse(shape,widths,k);
    report.geometry=struct('coreSize',size(g.core.X),'polarSize',size(g.polar.X), ...
        'activeDofCount',g.activeDofCount,'level',p.level);

    % A relative superlevel region must not change under positive amplitude
    % scaling; its peak and threshold scale by exactly the same factor.
    scaled=region(state(g,widths,7.5),g,p);
    report.amplitude=struct('factor',7.5, ...
        'extentDifference',max(abs(scaled.extents-shape.extents)), ...
        'peakScalingError',abs(scaled.peak/shape.peak-7.5));
    assert(report.amplitude.extentDifference<1e-9&& ...
        report.amplitude.peakScalingError<1e-10,'ipm_cr:RegionTestAmplitude');

    % Below-threshold changes accumulate against the retained reference,
    % rather than being discarded by comparing only consecutive checks.
    first=region(state(g,widths.*[1.10,1],1),g,p);
    second=region(state(g,widths.*[1.21,1],1),g,p);
    height=region(state(g,widths.*[1,1.25],1),g,p);
    shrinking=region(state(g,widths.*[1,.78],1),g,p);
    report.widthFirst=check_ellipse(first,widths.*[1.10,1],k);
    report.widthSecond=check_ellipse(second,widths.*[1.21,1],k);
    report.height=check_ellipse(height,widths.*[1,1.25],k);
    report.shrinking=check_ellipse(shrinking,widths.*[1,.78],k);
    d1=abs(first.extents./shape.extents-1);
    step2=abs(second.extents./first.extents-1);
    cumulative=abs(second.extents./shape.extents-1);
    dh=abs(height.extents./shape.extents-1);ds=abs(shrinking.extents./shape.extents-1);
    report.cumulative=struct('firstChange',d1,'secondIncrement',step2, ...
        'secondVsOriginal',cumulative,'heightChange',dh,'shrinkChange',ds, ...
        'threshold',p.threshold,'controllerInvoked',false);
    assert(max(d1)<p.threshold&&max(step2)<p.threshold&& ...
        cumulative(1)>p.threshold&&dh(2)>p.threshold&&ds(2)>p.threshold, ...
        'ipm_cr:RegionTestCumulativeGeometry');

    % Ghost data are not part of the mathematical observation.
    polluted=base;polluted.rhs.q.core(g.core.ghost)=NaN;
    polluted.rhs.q.polar(g.polar.ghost)=Inf;cleaned=region(polluted,g,p);
    report.ghostDifference=max(abs([cleaned.extents,cleaned.xBounds,cleaned.yBounds,cleaned.peak]- ...
        [shape.extents,shape.xBounds,shape.yBounds,shape.peak]));
    assert(report.ghostDifference==0&&~cleaned.ghostDonorsUsed,'ipm_cr:RegionTestGhosts');

    % Put a resolved, below-global-maximum positive island at an actual
    % outer native polar point. It must be monitored but not merged into
    % the main wall-connected component's bounds.
    ir=find(g.polar.r>=3,1);it=max(2,round(numel(g.polar.theta)/5));
    center=[g.polar.X(ir,it),g.polar.Y(ir,it)];sigma=.08*g.polar.r(ir);
    island=base;
    for patch={'core','polar'}
        key=patch{1};z=g.(key);
        island.rhs.q.(key)=island.rhs.q.(key)+.6* ...
            exp(-((z.X-center(1)).^2+(z.Y-center(2)).^2)/sigma^2);
    end
    separated=region(island,g,p);
    report.island=struct('center',center,'width',sigma,'peak',separated.peak, ...
        'otherPositiveNodes',separated.otherPositiveNodes, ...
        'extentDifference',max(abs(separated.extents-shape.extents)), ...
        'extents',separated.extents,'wallConnected',separated.wallConnected);
    assert(separated.otherPositiveNodes>0&&separated.wallConnected&& ...
        report.island.extentDifference<.01*max(shape.extents), ...
        'ipm_cr:RegionTestSeparateIsland');
    assert(toc(timer)<=budget,'ipm_cr:RegionTestBudget');
    report.passed=true;report.status='passed_static_region_geometry';
catch problem
    report.status='failed';report.failure=struct('identifier',problem.identifier, ...
        'message',problem.message,'stack',problem.stack);persist();rethrow(problem)
end
persist();
    function persist()
        report.elapsedSeconds=toc(timer);ipm_cr.io.json(file,report);
    end
end

function v=region(s,g,p)
v=ipm_cr.remesh('region',s,g,p);
assert(v.valid,'ipm_cr:RegionTestValid','Synthetic region is invalid: %s',v.reason);
assert(v.wallConnected,'ipm_cr:RegionTestWallConnected');
end
function s=state(g,widths,amplitude)
q=struct();
for patch={'core','polar'}
    key=patch{1};z=g.(key);
    q.(key)=amplitude*exp(-(z.X/widths(1)).^2-(z.Y/widths(2)).^2);
end
s=struct('rhs',struct('q',q));
end
function r=check_ellipse(shape,widths,k)
expected=widths*k;
r=struct('expectedExtents',expected,'measuredExtents',shape.extents, ...
    'relativeExtentError',abs(shape.extents./expected-1), ...
    'xBounds',shape.xBounds,'yBounds',shape.yBounds,'peak',shape.peak);
assert(all(r.relativeExtentError<.02)&& ...
    max(abs(shape.xBounds-[-expected(1),expected(1)]))/expected(1)<.02&& ...
    abs(shape.yBounds(1))<1e-12,'ipm_cr:RegionTestEllipse', ...
    'Piecewise-linear subcell contour must resolve analytic ellipse extents to two percent.');
end
function geo=geometry()
L=40;a=1.2;
geo=struct('coreHalfWidth',a,'coreHeight',a,'coreSpacing',a/L, ...
    'polarStart',.65,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',4*L+1,'ntheta',2*L+1,'rmax',1e10,'firstSpacing',a/(2*L), ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01,'stencilWidth',7), ...
    'coreMap',struct('enabled',true,'intervals',[2*L,L],'minSpacing',[a/L,a/L]/sinh(1)));
end
