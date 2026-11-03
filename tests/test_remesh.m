function report=test_remesh(outputDir,controls)
%TEST_REMESH Active transfer, same-time redistribution and cold replay.
% Uses only the public initial datum and a small mapped production geometry.
% No time evolution or private saved trajectory is required.
if nargin<2,controls=struct();end
assert(isstruct(controls)&&isscalar(controls)&& ...
    all(ismember(fieldnames(controls),{'wallSeconds'})),'ipm_cr:RemeshTestOptions');
budget=180;if isfield(controls,'wallSeconds'),budget=controls.wallSeconds;end
validateattributes(budget,{'numeric'},{'scalar','real','finite','positive'});
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','remesh');end
outputDir=char(java.io.File(char(outputDir)).getCanonicalPath());
if ~isfolder(outputDir),mkdir(outputDir);end
path=fullfile(outputDir,'remesh.json');
assert(~isfile(path),'ipm_cr:RemeshTestFresh','Use a fresh output directory.');
timer=tic;report=struct('kind','portable_same_time_remesh_regression_v3', ...
    'status','running','passed',false,'timeAdvanced',0, ...
    'continuumAccuracyCertified',false,'adaptiveTrajectoryQualified',false);
try
    p=ipm_cr.remesh('config');
    assert(p.enabled&&p.level==.25&&p.threshold==.2&&p.checkEveryTau==.5&& ...
        p.minCheckEveryTau==.005&&numel(fieldnames(p))==9&&~isfield(p,'intervalTau')&& ...
        isequaln(p,ipm_cr.remesh('config',p)),'ipm_cr:RemeshTestDefaults');
    reject(@()ipm_cr.remesh('config',struct('typo',1)),'ipm_cr:RemeshOptions');
    reject(@()ipm_cr.remesh('config',struct('blend',2)),'ipm_cr:RemeshOptions');
    geo=geometry(16);g=ipm_cr.mesh.build(geo);a=axes_of(g);
    a.x=a.x+1e-4*(a.x/1.2).*(1-(a.x/1.2).^2);
    a.y=a.y+1e-4*(a.y/1.2).*(1-a.y/1.2);
    a.dxMetric=metric(a.x,2)';a.dyMetric=metric(a.y,1);
    moved=ipm_cr.mesh.build(struct('axes',a));

    % Endpoint continuation is explicit and still reproduces degree six.
    % Ordinary sampling/source maps must continue to reject extrapolation.
    nodes=linspace(-1,1,13)';queries=[-1.02;.17;1.02];
    reject(@()ipm_cr.mesh.localMap(nodes,queries,7),'ipm_cr:CarrierInterpolation');
    M=ipm_cr.mesh.localMap(nodes,queries,7,true);
    report.endpointPolynomialError=max(abs(M*poly(nodes,0)-poly(queries,0)));
    assert(report.endpointPolynomialError<2e-12,'ipm_cr:RemeshTestPolynomial');

    % The two evolved patches deliberately disagree. A transfer may not
    % average/project them or select one owner throughout their overlap.
    fixture=synthetic(g);
    fixture.W.core(:)=.125;fixture.W.polar(:)=.375;
    fixture.W.core(g.core.ghost)=NaN;fixture.W.polar(g.polar.ghost)=Inf;
    [transferred,audit]=ipm_cr.remesh('transfer',fixture,g,moved);
    inside=moved.core.active& ...
        moved.core.X>=g.core.activeXBounds(1)&moved.core.X<=g.core.activeXBounds(2)& ...
        moved.core.Y>=g.core.activeYBounds(1)&moved.core.Y<=g.core.activeYBounds(2);
    exposed=moved.core.active&~inside;
    assert(any(exposed,'all')&&all(transferred.W.core(moved.core.active)==.125)&& ...
        all(transferred.W.polar(moved.polar.active)==.375)&&~audit.ghostValuesRead&& ...
        ~audit.activeOverlapProjected&&audit.routes.core.crossPatch==0&& ...
        audit.routes.core.samePatchBoundaryExtension==nnz(exposed)&& ...
        audit.routes.polar.crossPatch==0, ...
        'ipm_cr:RemeshTestIndependentPatches');
    assert(isempty(transferred.ctxId)&&isempty(transferred.ctxKey)&& ...
        isequaln(transferred.scale,fixture.scale)&&transferred.step==fixture.step&& ...
        isequaln(transferred.originalAuthority,fixture.originalAuthority), ...
        'ipm_cr:RemeshTestTransferIdentity');
    report.activeTransfer=audit;

    polynomial=synthetic(g);polynomial.W.core=poly(g.core.X,g.core.Y);
    polynomial.W.polar=.25+poly(g.polar.X/1e10,g.polar.Y/1e10);
    polynomial.W.core(g.core.ghost)=NaN;polynomial.W.polar(g.polar.ghost)=Inf;
    pt=ipm_cr.remesh('transfer',polynomial,g,moved);
    exact=poly(moved.core.X,moved.core.Y);
    report.samePatchPolynomialError=max(abs(pt.W.core(moved.core.active)-exact(moved.core.active)));
    assert(report.samePatchPolynomialError<2e-12,'ipm_cr:RemeshTestPolynomial');
    % A resolved nonpolynomial field must continue through the old active
    % boundary from its own chart, even when the other chart differs.
    resolved=synthetic(g);resolved.W.core=profile(g.core.X,g.core.Y);
    resolved.W.polar=.25+profile(g.polar.X,g.polar.Y);
    resolved.W.core(g.core.ghost)=NaN;resolved.W.polar(g.polar.ghost)=Inf;
    rt=ipm_cr.remesh('transfer',resolved,g,moved);
    exact=profile(moved.core.X,moved.core.Y);error=abs(rt.W.core-exact);
    report.resolvedNonpolynomialError=max(error(moved.core.active));
    report.newlyExposedNonpolynomialError=max(error(exposed));
    assert(report.resolvedNonpolynomialError<1e-8, ...
        'ipm_cr:RemeshTestResolvedContinuation');
    polynomial.rhs.q=struct('core',ones(size(g.core.X)),'polar',ones(size(g.polar.X)));
    constantMonitor=ipm_cr.remesh('monitor',polynomial,g);
    polynomial.rhs.q.core=exp(-((g.core.X+.25)/.08).^2-(g.core.Y/.07).^2);
    localizedMonitor=ipm_cr.remesh('monitor',polynomial,g);
    assert(constantMonitor.indicator<32*eps&&localizedMonitor.indicator>1e-2, ...
        'ipm_cr:RemeshTestMonitor','The monitor must distinguish a constant from an unresolved source.');
    report.monitor=struct('constant',constantMonitor.indicator,'localized',localizedMonitor.indicator);
    % Check cadence does not cause interpolation. Changes accumulate against
    % the last actual remesh reference, not against the preceding check.
    shapeState=region_state(g,1);shapeContext=struct('mesh',g);
    [shapeState,~,calibration]=ipm_cr.remesh(shapeState,shapeContext,struct());
    assert(strcmp(calibration.reason,'initial_calibration')&&~calibration.remeshed&& ...
        shapeState.remesh.remeshCount==0&&shapeState.remesh.checks==1&& ...
        shapeState.remesh.nextCheckTau==p.checkEveryTau&& ...
        isequaln(shapeState.remesh.lastShape,shapeState.remesh.referenceShape), ...
        'ipm_cr:RemeshTestCalibration');
    fixedShape=shapeState.remesh.referenceShape;
    for tau=[.1,.5,1,2,10,100]
        unchanged=shapeState;unchanged.scale.canonicalTime=tau;
        [checked,sameContext,event]=ipm_cr.remesh(unchanged,shapeContext,struct());
        assert(~event.remeshed&&isequaln(sameContext,shapeContext)&& ...
            isequaln(checked.W,unchanged.W)&&isequaln(checked.rhs,unchanged.rhs)&& ...
            isequaln(checked.remesh.referenceShape,fixedShape)&&checked.remesh.remeshCount==0, ...
            'ipm_cr:RemeshTestConstantShape','Unchanged shape must not cause a time-driven transfer.');
        if ~event.checked
            assert(isequaln(checked,unchanged),'ipm_cr:RemeshTestCadenceState');
        else
            assert(checked.remesh.nextCheckTau==tau+p.checkEveryTau&& ...
                isequaln(checked.remesh.lastShape,fixedShape), ...
                'ipm_cr:RemeshTestFixedCadence','A fixed shape uses the maximum check interval.');
        end
        shapeState=checked;
    end
    near=region_state(g,1.18);near.remesh=shapeState.remesh;
    near.scale.canonicalTime=near.remesh.lastCheckTau+p.checkEveryTau;
    [nearChecked,~,nearEvent]=ipm_cr.remesh(near,shapeContext,struct());
    expectedCadence=max(p.minCheckEveryTau,min(p.checkEveryTau, ...
        .8*max(0,p.threshold-max(nearEvent.directionalChange))/max(nearEvent.directionalRateEstimate)));
    assert(nearEvent.checked&&~nearEvent.triggered&&~nearEvent.remeshed&& ...
        nearEvent.checkCadence<p.checkEveryTau&&nearEvent.checkCadence>=p.minCheckEveryTau&& ...
        nearEvent.checkCadence==expectedCadence&& ...
        isequaln(nearChecked.remesh.referenceShape,fixedShape)&& ...
        nearChecked.remesh.nextCheckTau==near.scale.canonicalTime+expectedCadence, ...
        'ipm_cr:RemeshTestAdaptiveCadence','Approaching the threshold shortens checks without interpolation.');
    premature=nearChecked;
    premature.scale.canonicalTime=premature.remesh.nextCheckTau-p.minCheckEveryTau/2;
    [early,~,earlyEvent]=ipm_cr.remesh(premature,shapeContext,struct());
    assert(~earlyEvent.checked&&isequaln(early,premature),'ipm_cr:RemeshTestCadenceState');
    report.adaptiveCadence=struct('movingChange',nearEvent.directionalChange, ...
        'movingRate',nearEvent.directionalRateEstimate,'cadence',nearEvent.checkCadence, ...
        'nextCheckTau',nearChecked.remesh.nextCheckTau,'remeshed',false);
    changes=zeros(3,2);shapeScales=[1.10,1.15,1.25];
    for k=1:3
        test=region_state(g,shapeScales(k));
        test.remesh=shapeState.remesh;test.scale.canonicalTime=100+k*p.checkEveryTau;
        [checked,~,event]=ipm_cr.remesh(test,shapeContext,struct('candidateAxes',axes_of(g)));
        changes(k,:)=event.directionalChange;
        assert(event.triggered==(k==3)&&~event.remeshed&& ...
            isequaln(checked.remesh.referenceShape,fixedShape)&&checked.remesh.remeshCount==0, ...
            'ipm_cr:RemeshTestAccumulatedShape');
        shapeState=checked;
    end
    assert(max(changes(1,:))<p.threshold&&max(changes(2,:))<p.threshold&& ...
        max(changes(3,:))>=p.threshold,'ipm_cr:RemeshTestAccumulatedShape');
    report.shapeTrigger=struct('initialInterpolation',false,'unchangedTimes',[.1,.5,1,2,10,100], ...
        'directionalChanges',changes,'referenceUpdatedOnCheck',false);
    assert(toc(timer)<budget,'ipm_cr:RemeshTestBudget');

    options=struct('geometry',geo,'poissonBackend','fem', ...
        'solverOptions',struct('rOuter',1e10),'remesh',struct('enabled',false));
    [s,c]=ipm_cr.init(options);
    beforeStats=c.stats();
    [disabled,disabledC,off]=ipm_cr.remesh(s,c,struct('enabled',false));
    assert(isequaln(disabled,s)&&disabledC==c&&strcmp(off.reason,'disabled'), ...
        'ipm_cr:RemeshTestDisabled');

    % Geometry input errors are tested at their own interface, not through
    % remesh rejection/recovery machinery.
    unordered=a;unordered.x(2)=unordered.x(3)+.001;
    reject(@()ipm_cr.mesh.build(struct('axes',unordered)),'MATLAB:expectedIncreasing');
    [s,sameC,calibration]=ipm_cr.remesh(s,c,struct());
    assert(sameC==c&&strcmp(calibration.reason,'initial_calibration')&& ...
        ~calibration.remeshed&&isequaln(beforeStats,c.stats())&& ...
        s.remesh.remeshCount==0,'ipm_cr:RemeshTestCalibration');
    % Model the same current Gaussian after expansion from a smaller earlier
    % calibrated region. The ordinary threshold path must really transfer.
    expanded=s;expanded.scale.canonicalTime=s.scale.canonicalTime+p.checkEveryTau;
    for name={'xBounds','yBounds','widths','extents','peakXY','ownedPeakXY','monitorRadius'}
        key=name{1};expanded.remesh.referenceShape.(key)=expanded.remesh.referenceShape.(key)/1.25;
    end
    [automatic,autoC,autoEvent]=ipm_cr.remesh(expanded,c,struct());
    autoStats=autoC.stats();actualRegion=ipm_cr.remesh('region',automatic,autoC.mesh);
    actualRegion=rmfield(actualRegion,'crossings');
    assert(autoEvent.triggered&&autoEvent.remeshed&&~autoEvent.forced&&autoC~=c&& ...
        autoStats.rhsPoissonCalls==1&& ...
        isequaln(automatic.scale,expanded.scale)&& ...
        isequaln(automatic.originalAuthority,expanded.originalAuthority)&& ...
        isequaln(automatic.remesh.referenceShape,actualRegion)&& ...
        automatic.remesh.remeshCount==1&&isequaln(beforeStats,c.stats()), ...
        'ipm_cr:RemeshTestAutomatic','A resolved 25%% accumulated shape change must execute one fresh-RHS remesh.');
    report.automaticThresholdEvent=autoEvent;
    clear automatic autoC
    [after,newC,event]=ipm_cr.remesh(s,c,struct(),true);
    report.event=event;persist();
    assert(event.remeshed&&strcmp(event.status,'remeshed')&&event.forced&&newC~=c&& ...
        strcmp(after.ctxId,newC.id)&&isequaln(after.ctxKey,newC.key)&& ...
        isequaln(newC.opts,c.opts),'ipm_cr:RemeshTestMove','The explicit forced pass did not redistribute.');
    assert(isequaln(after.scale,s.scale)&&after.step==s.step&& ...
        isequaln(after.originalAuthority,s.originalAuthority)&& ...
        isequaln(newC.wall,c.wall)&&after.remesh.remeshCount==1&& ...
        after.remesh.lastRemeshTau==s.scale.canonicalTime&& ...
        isequaln(newC.mesh.polar.r,g.polar.r)&&isequaln(newC.mesh.polar.theta,g.polar.theta)&& ...
        isequaln(newC.mesh.polar.drMetric,g.polar.drMetric)&& ...
        isequal(newC.mesh.core.x([1,end]),g.core.x([1,end]))&& ...
        isequal(newC.mesh.core.y([1,end]),g.core.y([1,end]))&& ...
        newC.mesh.core.x(newC.mesh.core.originIndex)==0, ...
        'ipm_cr:RemeshTestFixedState');
    stats=newC.stats();
    assert(event.timeAdvanced==0&&after.rhs.diagnosticsComplete&& ...
        stats.rhsPoissonCalls==1&&after.rhs.poisson.sameNativeGrid&& ...
        ~after.rhs.poisson.independentPoissonGridUsed,'ipm_cr:RemeshTestFreshRhs');
    assert(isequaln(beforeStats,c.stats()),'ipm_cr:RemeshTestOldContextUntouched');
    [notYet,sameC,skip]=ipm_cr.remesh(after,newC,struct('enabled',true));
    assert(isequaln(notYet,after)&&sameC==newC&&strcmp(skip.reason,'check_cadence'), ...
        'ipm_cr:RemeshTestInterval','No due shape check may alter the accepted state.');
    unchangedAxes=axes_of(newC.mesh);unchangedAxes.dxMetric(1)=unchangedAxes.dxMetric(1)+eps(unchangedAxes.dxMetric(1));
    priorStats=newC.stats();
    [unchanged,sameC,same]=ipm_cr.remesh(after,newC,struct('candidateAxes',unchangedAxes),true);
    unchangedPrimary=unchanged;unchangedPrimary.remesh=after.remesh;
    assert(isequaln(unchangedPrimary,after)&&sameC==newC&&strcmp(same.reason,'geometry_unchanged')&& ...
        isequaln(unchanged.remesh.referenceShape,after.remesh.referenceShape)&& ...
        unchanged.remesh.checks==after.remesh.checks+1&& ...
        isequaln(priorStats,newC.stats()),'ipm_cr:RemeshTestExactNodes', ...
        'Identical nodes must retain the original metrics and avoid a new RHS.');
    [after,newC,followup]=ipm_cr.remesh(after,newC,struct(),true);
    assert(followup.remeshed&&followup.forced&&after.remesh.remeshCount==2&& ...
        after.remesh.lastRemeshTau==s.scale.canonicalTime&&isequaln(after.scale,s.scale), ...
        'ipm_cr:RemeshTestForced','An explicit forced pass must adapt at the same physical time.');
    report.controller=struct('firstEndpointCalibratedOnly',true,'checkEveryTau',p.checkEveryTau, ...
        'threshold',p.threshold,'forcedFollowupAdapted',true,'exactNodesRetainMetrics',true);
    assert(toc(timer)<budget,'ipm_cr:RemeshTestBudget');

    % Repeated same-time adaptation of the exact analytic Gaussian must
    % approach one reference-based density, not multiply the current mesh
    % density on every call. No runtime minimum-spacing gate is involved.
    reference=after.remesh.referenceCore;referenceMin=[min(diff(reference.x)),min(diff(reference.y))];
    motion=zeros(20,1);minimum=zeros(20,2);statuses=cell(20,1);
    for k=1:numel(motion)
        oldG=newC.mesh;
        [next,nextC,move]=ipm_cr.remesh(after,newC,struct(),true);
        newG=nextC.mesh;motion(k)=max([max(abs(newG.core.x-oldG.core.x)), ...
            max(abs(newG.core.y-oldG.core.y))]);
        minimum(k,:)=[min(diff(newG.core.x)),min(diff(newG.core.y))];statuses{k}=move.status;
        assert(isequaln(next.scale,s.scale)&&next.step==s.step&& ...
            isequaln(next.originalAuthority,s.originalAuthority)&& ...
            isequaln(next.referencePeak,s.referencePeak)&& ...
            isequaln(next.remesh.referenceCore,reference)&& ...
            isequal(size(next.W.core),size(s.W.core))&& ...
            isequal(size(next.W.polar),size(s.W.polar))&& ...
            all(next.eta.core==0,'all')&&all(next.eta.polar==0,'all'), ...
            'ipm_cr:RemeshTestRepeatedFixedState','Remeshing changed the same-time analytic state.');
        after=next;newC=nextC;
    end
    assert(all(minimum>=referenceMin/(1+p.monitorStrength),'all')&& ...
        (motion(end)<motion(1)||motion(end)<=64*eps(1.2)), ...
        'ipm_cr:RemeshTestReferenceDensity','Repeated fixed-source remeshing drifted or did not settle.');
    report.fixedSource=struct('iterations',numel(motion),'nodeMotion',motion, ...
        'minimumSpacing',minimum,'referenceMinimumSpacing',referenceMin, ...
        'statuses',{statuses},'timeAdvanced',0,'runtimeSpacingGateUsed',false);
    assert(toc(timer)<budget,'ipm_cr:RemeshTestBudget');

    % Serialize only data. A cold context must reconstruct every new mesh
    % plan from its actual axes while retaining the initial mesh prescription.
    payload=struct('state',after,'geometry',newC.mesh,'wall',newC.wall,'options',newC.opts);
    checkpoint=fullfile(outputDir,'remeshed_state.mat');save(checkpoint,'payload','-v7.3');
    loaded=load(checkpoint,'payload');cold=loaded.payload;
    assert(isequaln(payload,cold),'ipm_cr:RemeshTestPersistence');
    coldC=ipm_cr.context(cold.geometry,cold.wall,cold.options);
    replay=ipm_cr.refresh(cold.state,coldC);
    errors=struct('W',distance(after.W,replay.W),'q',distance(after.rhs.q,replay.rhs.q), ...
        'F',distance(after.rhs.F,replay.rhs.F),'Q',distance(after.rhs.Q,replay.rhs.Q));
    assert(errors.W==0&&errors.q==0&&max(errors.F,errors.Q)<1e-10&& ...
        isequaln(after.scale,replay.scale)&&isequaln(after.remesh,replay.remesh), ...
        'ipm_cr:RemeshTestColdRestore');
    report.cold=struct('checkpointSha256',ipm_cr.io.sha(checkpoint),'differences',errors, ...
        'controllerBitwise',true,'allScaleClockWordsBitwise',true);
    replay.scale.canonicalTime=replay.scale.canonicalTime+p.checkEveryTau;
    [checked,sameC,event]=ipm_cr.remesh(replay,coldC,struct());
    assert(~event.remeshed&&event.checked&&sameC==coldC&& ...
        isequaln(checked.remesh.referenceShape,replay.remesh.referenceShape)&& ...
        checked.remesh.lastCheckTau==replay.scale.canonicalTime&& ...
        checked.remesh.checks==replay.remesh.checks+1&& ...
        checked.remesh.nextCheckTau==replay.scale.canonicalTime+p.checkEveryTau&& ...
        isequaln(checked.remesh.lastShape,replay.remesh.lastShape), ...
        'ipm_cr:RemeshTestColdController');
    assert(toc(timer)<budget,'ipm_cr:RemeshTestBudget');
    report.status='passed';report.passed=true;
catch problem
    report.status='failed';report.failure=struct('identifier',problem.identifier, ...
        'message',problem.message,'stack',problem.stack);persist();rethrow(problem)
end
persist();
    function persist()
        report.elapsedSeconds=toc(timer);ipm_cr.io.json(path,report);
    end
end

function geo=geometry(L)
a=1.2;geo=struct('coreHalfWidth',a,'coreHeight',a,'coreSpacing',a/L, ...
    'polarStart',.65,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',4*L+1,'ntheta',2*L+1,'rmax',1e10,'firstSpacing',a/(2*L), ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01,'stencilWidth',7), ...
    'coreMap',struct('enabled',true,'intervals',[2*L,L],'minSpacing',[a/L,a/L]/sinh(1)));
end

function a=axes_of(g)
a=struct('x',g.core.x,'y',g.core.y,'r',g.polar.r,'theta',g.polar.theta, ...
    'dxMetric',g.core.dxMetric,'dyMetric',g.core.dyMetric,'drMetric',g.polar.drMetric);
end

function value=metric(axis,logicalLength)
value=ipm_cr.mesh.fd(linspace(0,logicalLength,numel(axis)),1,7)*axis(:)*logicalLength/(numel(axis)-1);
end

function s=synthetic(g)
s=struct('scale',struct('canonicalTime',0,'physicalTime',1,'physicalTimeLow',1e-20, ...
    'physicalTimeTail',[1e-40,1e-60],'densityOrigin',0), ...
    'step',4,'originalAuthority',struct('token',23),'referencePeak',.6, ...
    'W',struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X))), ...
    'ctxId','synthetic','ctxKey',struct('token',7),'rhs',struct('c_l',.2,'c_r',.3));
end

function s=region_state(g,scale)
s=synthetic(g);
for patch={'core','polar'}
    key=patch{1};p=g.(key);
    s.rhs.q.(key)=exp(-(p.X/(.28*scale)).^2-(p.Y/(.22*scale)).^2);
end
end

function v=poly(x,y)
v=.125+.01*x-.02*y+.003*x.*y+.001*x.^2+.0002*x.^6-.0003*y.^5+.0001*(x.*y).^3;
end

function v=profile(x,y)
v=.125+.02*sin(.8*x+.5*y).*exp(-.15*(x.^2+y.^2));
end

function value=distance(a,b)
value=max([max(abs(a.core-b.core),[],'all'),max(abs(a.polar-b.polar),[],'all')]);
end

function reject(callback,identifier)
try
    callback();
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:RemeshTestRejection', ...
        'Expected %s, received %s.',identifier,problem.identifier);return
end
error('ipm_cr:RemeshTestRejection','Expected rejection %s.',identifier);
end
