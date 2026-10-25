function report=test_ipm_cr_gauge_acceleration(outputDir,options)
%TEST_IPM_CR_GAUGE_ACCELERATION Frozen-state branch Jacobian and legacy gauge.
% Default execution uses only the analytic initial state. Historical tau200
% checks require includeTau200=true and an explicit checkpointFolder.
% No time is advanced. Both gauge solvers use exactly the same synchronized
% density, already computed Poisson field, wall samples and transport
% callback. Central differences test the complete transport arrays. Nodes
% whose upwind branch changes are recorded separately, never called smooth.
if nargin<2,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir)
    outputDir=fullfile(root,'results','gauge_acceleration');
end
o=resolve(options);outputDir=char(outputDir);
assert(java.io.File(outputDir).isAbsolute(),'ipm_cr:GaugeTestOutput','Use an absolute output directory.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'gauge_acceleration.json')), ...
    'ipm_cr:GaugeTestPreserve','Choose a fresh output directory.');
timer=tic;
report=struct('kind','frozen_density_transport_partials_and_legacy_IPC_comparison_v1', ...
    'status','running','options',o,'timeAdvanced',0,'readOnlyCheckpointInspection',false, ...
    'oldTrajectorySourceIdentityClaimed',false,'records',struct([]), ...
    'checks',struct('initial',struct('requested',o.includeInitial,'status','skipped','reason','Not requested'), ...
        'history',struct('requested',o.includeTau200,'status','skipped','reason','Historical checkpoint was not requested')), ...
    'requestedChecksSatisfied',false, ...
    'sourceInventory',source_inventory(),'sourceHashesUnchanged',false, ...
    'transportBranchDerivativeTolerance',1e-7,'legacyRootTolerance',1e-10, ...
    'legacyAuditTolerance',1e-8,'frozenEquivalenceTolerance',1e-8,'allChecksPassed',false,'elapsedSeconds',0, ...
    'limitations',{{'A central difference crossing a zero-speed upwind kink does not test a classical derivative; those nodes are listed separately.', ...
        'This is a frozen-state equivalence test, not a trajectory, speed benchmark, or fourth-order certificate.', ...
        'The old tau200 density is reevaluated with the current operators; its historical flow is not used as an exact reference.', ...
        'Both gauge solvers share one actual Poisson field per state; any common source or Poisson error is outside this comparison.'}});
if o.includeInitial,report.checks.initial.status='pending';report.checks.initial.reason='';end
if o.includeTau200,report.checks.history.status='pending';report.checks.history.reason='';end
persist(outputDir,report);
try
    if o.includeInitial
        [state,ctx]=ipm_cr.init(o.explicitOptions);
        g=ctx.mesh;wallOps=ctx.wall;method=ctx.opts;
        row=check_state('analytic_initial',state,g,wallOps,method,o);
        row.checkpoint=struct('used',false);report.records=append(report.records,row);
        report.checks.initial=struct('requested',true,'status','completed','reason','');
        report.elapsedSeconds=toc(timer);persist(outputDir,report);
    end
    historyAvailable=~isempty(o.checkpointFolder)&&isfile(fullfile(o.checkpointFolder,'latest.json'));
    if o.includeTau200&&historyAvailable
        [cp,pointer]=load_endpoint(o.checkpointFolder);
        report.readOnlyCheckpointInspection=true;
        beforeRho=cp.state.W;
        row=check_state('accepted_tau200_density',cp.state,cp.geometry,cp.wallOps,cp.explicitOptions,o);
        assert(isequaln(beforeRho,cp.state.W),'ipm_cr:GaugeTestMutation','The loaded accepted density changed.');
        assert_unchanged(o.checkpointFolder,pointer);
        row.checkpoint=struct('used',true,'pointer',pointer,'shaVerifiedBeforeAndAfter',true, ...
            'densityChanged',false,'timeAdvanced',0);
        report.records=append(report.records,row);
        report.checks.history=struct('requested',true,'status','completed','reason','');
    elseif o.includeTau200
        report.checks.history=struct('requested',true,'status','skipped', ...
            'reason','An explicit accepted checkpointFolder with latest.json is required; no historical test was executed');
    end
    report.sourceHashesUnchanged=isequaln(report.sourceInventory,source_inventory());
    assert(report.sourceHashesUnchanged,'ipm_cr:GaugeTestSources','Recorded numerical sources changed during testing.');
    report.allChecksPassed=~isempty(report.records)&&all([report.records.transportPartialChecksPassed])&& ...
        all([report.records.frozenGaugeEquivalencePassed]);
    report.requestedChecksSatisfied=(~o.includeInitial||strcmp(report.checks.initial.status,'completed'))&& ...
        (~o.includeTau200||strcmp(report.checks.history.status,'completed'));
    if isempty(report.records),report.status='skipped_no_executed_checks';
    elseif ~report.requestedChecksSatisfied,report.status='completed_with_skipped_requested_checks';
    elseif report.allChecksPassed,report.status='passed_frozen_state_comparisons';
    else,report.status='completed_measurements_with_failed_checks';end
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function o=resolve(options)
o=struct('includeInitial',true,'includeTau200',false,'explicitOptions',struct(), ...
    'checkpointFolder','', ...
    'partialSteps',[1e-4,1e-5,1e-6], ...
    'rateOffsets',[0,0;.13,.21;-.11,-.17]);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(o))), ...
    'ipm_cr:GaugeTestOptions','Unknown gauge qualification option.');
names=fieldnames(options);for j=1:numel(names),o.(names{j})=options.(names{j});end
validateattributes(o.includeInitial,{'logical'},{'scalar'});
validateattributes(o.includeTau200,{'logical'},{'scalar'});
validateattributes(o.partialSteps,{'numeric'},{'vector','positive','finite','real'});
o.partialSteps=reshape(o.partialSteps,1,[]);
assert(isnumeric(o.rateOffsets)&&size(o.rateOffsets,2)==2&&all(isfinite(o.rateOffsets),'all'), ...
    'ipm_cr:GaugeTestOptions','rateOffsets must have finite [c_l,delta] rows.');
assert(isstruct(o.explicitOptions)&&isscalar(o.explicitOptions),'ipm_cr:GaugeTestOptions');
assert(ischar(o.checkpointFolder)||(isstring(o.checkpointFolder)&&isscalar(o.checkpointFolder)), ...
    'ipm_cr:GaugeTestOptions','checkpointFolder must be an explicit path.');
o.checkpointFolder=char(o.checkpointFolder);
assert(isempty(o.checkpointFolder)||java.io.File(o.checkpointFolder).isAbsolute(), ...
    'ipm_cr:GaugeTestOptions','Use an absolute checkpointFolder.');
end

function row=check_state(name,state,g,wallOps,method,o)
timer=tic;rho=ipm_cr.mesh.fill(state.W,g);
guess=[state.rhs.c_l,state.rhs.c_r];
% Compute one current, complete physical state. The legacy comparison below
% consumes this exact q and Poisson result; it never recomputes either one.
ctx=ipm_cr.context(g,wallOps,method);
% This diagnostic intentionally tests the unsplit density operator even
% when its input accepted state also contains an analytic residual.
trial=state;if isfield(trial,'eta'),trial=rmfield(trial,'eta');end
trial.W=rho;trial.ctxId=ctx.id;trial.ctxKey=ctx.key;
fresh=ipm_cr.rhs(trial,ctx,true);
physical=fresh.physical;
[r,map]=wall_map(g,wallOps.x);
R=sample_wall(rho,map);R(1,:)=0;
[m,observe]=ipm_cr.gauge.measure(R,r,wallOps,struct('interpolation','linear'));
Uref=ipm_cr.mesh.sample(struct('core',physical.core.u1,'polar',physical.polar.u1),g,m.C,0);
prepared=ipm_cr.op.prepare(rho,g,struct('transport',method.transport));
[fast,fastInfo]=ipm_cr.gauge.rates(R,r,wallOps,@base,Uref,guess,m,observe);
[~,legacy,legacyInfo]=ipm_cr.gauge.seed(R,r,wallOps,@base,Uref,guess, ...
    struct('measurementOptions',struct('interpolation','linear'), ...
    'rootTolerance',1e-10,'auditTolerance',1e-8));
fastFull=complete(fast);legacyFull=complete(legacy);
partials=struct([]);
for j=1:size(o.rateOffsets,1)
    point=[fast.c_l,fast.delta]+o.rateOffsets(j,:);
    [vx,vy]=velocities(point(1),point(2));
    value=prepared(vx,vy,struct('cfl',false,'partials',true));analytic=value.partials;
    for k=1:2
        key='c_l';if k==2,key='delta';end
        for epsilon=o.partialSteps
            h=epsilon*max(1,abs(point(k)));plus=point;minus=point;
            plus(k)=plus(k)+h;minus(k)=minus(k)-h;
            [px,py]=velocities(plus(1),plus(2));[mx,my]=velocities(minus(1),minus(2));
            bp=prepared(px,py,struct('cfl',false));bm=prepared(mx,my,struct('cfl',false));
            bp=bp.B;bm=bm.B;
            fd=struct('core',(bp.core-bm.core)/(2*h),'polar',(bp.polar-bm.polar)/(2*h));
            masks=branch_masks(mx,my,px,py,g);
            difference=subtract(analytic.(key),fd);
            metrics=partial_metrics(difference,analytic.(key),masks,g);
            partials=append(partials,struct('point',point,'parameter',key,'relativeStep',epsilon, ...
                'actualStep',h,'metrics',metrics,'branchDiagnostics',masks.summary, ...
                'densityFrozen',true,'WenoWeightsFrozen',true, ...
                'classicalDerivativeClaimedAtKink',false));
        end
    end
end
names={'c_l','delta','c_r','c_omega','Uref'};a=zeros(1,numel(names));b=a;
for j=1:numel(names),a(j)=fast.(names{j});b(j)=legacy.(names{j});end
freshValues=[fresh.c_l,fresh.delta,fresh.c_r,fresh.c_omega,fresh.Uref];
row=struct('name',name,'tau',state.scale.canonicalTime,'step',state.step, ...
    'gridShapes',struct('core',size(g.core.X),'polar',size(g.polar.X)), ...
    'sameDensitySourcePoissonWallOps',true,'freshPoissonDiagnostics',fresh.poisson, ...
    'rateNames',{names},'fastRates',a,'legacyRates',b,'rateSignedDifference',a-b, ...
    'rateMixedRelativeDifference',abs(a-b)./max(1,abs(b)), ...
    'explicitEvaluateRateSignedDifference',a-freshValues, ...
    'fastInformation',fastInfo,'legacyInformation',compact_legacy(legacyInfo), ...
    'FComparison',field_metrics(subtract(fastFull.F,legacyFull.F),g), ...
    'explicitEvaluateFComparison',field_metrics(subtract(fastFull.F,fresh.F),g), ...
    'wallFunctionals',struct('I',m.I,'P',m.P,'C',m.C), ...
    'fastIPCResidual',fastFull.residual,'legacyIPCResidual',legacyFull.residual, ...
    'IPCResidualSignedDifference',fastFull.residual-legacyFull.residual, ...
    'fastNormalizedIPCRates',fastFull.normalizedRates,'legacyNormalizedIPCRates',legacyFull.normalizedRates, ...
    'cfl',struct('fastTransportRate',fastFull.rate,'legacyTransportRate',legacyFull.rate, ...
        'fastFullRate',fastFull.rate+abs(fast.c_omega),'legacyFullRate',legacyFull.rate+abs(legacy.c_omega), ...
        'fastPosition',fastFull.position,'legacyPosition',legacyFull.position, ...
        'bothActiveOverlapInteriorsIncluded',true, ...
        'signedFullRateDifference',fastFull.rate+abs(fast.c_omega)-legacyFull.rate-abs(legacy.c_omega)), ...
    'transportPartials',partials,'timeAdvanced',0,'seconds',toc(timer));
row.transportPartialChecksPassed=all(arrayfun(@(p)p.metrics.maximumSafeMixedError<=1e-7,partials));
Fscale=max([1;abs(fastFull.F.core(g.core.active));abs(fastFull.F.polar(g.polar.active))]);
row.FComparison.referenceScale=Fscale;
row.FComparison.mixedRelativeInfinity=row.FComparison.allActiveInfinity/Fscale;
row.frozenGaugeEquivalencePassed=all(abs(a-b)<=1e-8*max(1,abs(a)))&& ...
    row.FComparison.mixedRelativeInfinity<=1e-8&&max(abs(fastFull.residual))<=1e-8&& ...
    max(abs(legacyFull.residual))<=1e-8&&abs(row.cfl.signedFullRateDifference)<= ...
    1e-8*max(1,fastFull.rate+abs(fast.c_omega))&& ...
    all(abs(a-freshValues)<=1e-8*max(1,abs(a)))&& ...
    row.explicitEvaluateFComparison.allActiveInfinity<=1e-8*Fscale;
fprintf('GAUGE %s tau=%.6g rates %.3e F %.3e IPC fast/legacy %.3e/%.3e CFL-rate %.3e %.3fs\n', ...
    name,state.scale.canonicalTime,max(abs(a-b)),row.FComparison.allActiveInfinity, ...
    max(abs(fastFull.residual)),max(abs(legacyFull.residual)),row.cfl.signedFullRateDifference,row.seconds);

    function [vx,vy]=velocities(cl,delta)
        vx=struct('core',(physical.core.u1-Uref)+cl*g.core.X+delta, ...
            'polar',(physical.polar.u1-Uref)+cl*g.polar.X+delta);
        vy=struct('core',physical.core.u2+cl*g.core.Y, ...
            'polar',physical.polar.u2+cl*g.polar.Y);
    end
    function [wall,jac]=base(cl,delta)
        [vx,vy]=velocities(cl,delta);
        value=prepared(vx,vy,struct('cfl',false,'partials',nargout>1));
        B=value.B;jac=value.partials;
        wall=sample_wall(ipm_cr.mesh.fill(B,g),map);wall(1,:)=wall(1,1);
        if nargout>1
            for field={'c_l','delta'}
                key=field{1};v=ipm_cr.mesh.fill(jac.(key),g);
                jac.(key)=sample_wall(v,map);jac.(key)(1,:)=jac.(key)(1,1);
            end
        end
    end
    function result=complete(rates)
        [vx,vy]=velocities(rates.c_l,rates.delta);
        value=prepared(vx,vy);B=value.B;rate=value.rate;position=value.position;
        origin=B.core(1,g.core.originIndex);
        F=struct('core',B.core+rates.c_omega*rho.core-origin, ...
            'polar',B.polar+rates.c_omega*rho.polar-origin);
        F=ipm_cr.mesh.fill(F,g);F.core(1,g.core.originIndex)=0;
        wall=sample_wall(F,map);wall(1,:)=0;
        observed=ipm_cr.gauge.measure(R,r,wallOps,struct('interpolation','linear','rateField',wall));
        target=[0;0;wallOps.rescaling.peakPositionRestoringRate* ...
            (m.C-wallOps.rescaling.referenceWallGradientCentroid)];
        result=struct('F',F,'rate',rate,'position',position, ...
            'normalizedRates',observed.rate.normalized,'residual',observed.rate.normalized+target);
    end
end

function masks=branch_masks(mx,my,px,py,g)
my.core(1,:)=0;py.core(1,:)=0;my.polar(:,[1,end])=0;py.polar(:,[1,end])=0;
co=g.polar.cosTheta;si=g.polar.sinTheta;
mr=mx.polar.*co+my.polar.*si;pr=px.polar.*co+py.polar.*si;
ma=(-mx.polar.*si+my.polar.*co)./g.polar.r;
pa=(-px.polar.*si+py.polar.*co)./g.polar.r;
ma(:,[1,end])=0;pa(:,[1,end])=0;
masks.core=same_branch(mx.core,px.core)&same_branch(my.core,py.core);
masks.polar=same_branch(mr,pr)&same_branch(ma,pa);
outerChanged=~same_branch(mr(end,:),pr(end,:));
masks.polar(:,outerChanged)=false;
masks.summary=struct('coreSafeNodeCount',nnz(masks.core),'coreKinkOrBranchChangeCount',nnz(~masks.core), ...
    'polarSafeNodeCount',nnz(masks.polar),'polarKinkOrBranchChangeCount',nnz(~masks.polar), ...
    'outerRadialBranchChangedColumnCount',nnz(outerChanged), ...
    'unchangedZeroSpeedsTreatedAsConstantInThisParameterDirection',true, ...
    'criterion','Both directional component speeds retain sign at +/-h; unchanged exact zeros are safe for this direction. Changed exterior inflow columns are excluded in full.');
end

function safe=same_branch(a,b)
safe=(a>0&b>0)|(a<0&b<0)|(a==0&b==0);
end

function info=partial_metrics(error,analytic,masks,g)
info=struct('patches',struct(),'maximumSafeMixedError',0);
for patch={'core','polar'}
    key=patch{1};p=g.(key);e=error.(key);safe=masks.(key);
    assert(all(isfinite(e),'all'),'ipm_cr:GaugeTestFinite');
    row=struct('allNodes',metric(e,true(size(safe)),p.X,p.Y), ...
        'allActive',metric(e,p.active,p.X,p.Y), ...
        'safeAllNodes',metric(e,safe,p.X,p.Y), ...
        'kinkOrBranchChangeNodes',metric(e,~safe,p.X,p.Y));
    scale=max([1;abs(analytic.(key)(safe))]);
    row.safeAnalyticScale=scale;row.safeMixedError=row.safeAllNodes.absoluteInfinity/scale;
    info.maximumSafeMixedError=max(info.maximumSafeMixedError,row.safeMixedError);
    info.patches.(key)=row;
end
end

function info=field_metrics(field,g)
info=struct('patches',struct(),'allActiveInfinity',0);
for patch={'core','polar'}
    key=patch{1};p=g.(key);row=struct('allNodes',metric(field.(key),true(size(p.active)),p.X,p.Y), ...
        'allActive',metric(field.(key),p.active,p.X,p.Y));
    info.patches.(key)=row;info.allActiveInfinity=max(info.allActiveInfinity,row.allActive.absoluteInfinity);
end
end

function info=metric(value,mask,x,y)
indices=find(mask);info=struct('sampleCount',numel(indices),'absoluteInfinity',0,'rms',0,'maximumLocation',struct());
if isempty(indices),return;end
v=value(indices);[maximum,k]=max(abs(v));index=indices(k);
info.absoluteInfinity=maximum;info.rms=sqrt(mean(v.^2));
info.maximumLocation=struct('X',x(index),'Y',y(index),'r',hypot(x(index),y(index)));
end

function [r,map]=wall_map(g,x)
r=unique([0;abs(x(:))]);X=[r,-r];Y=zeros(size(X));
empty=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
[~,map]=ipm_cr.mesh.sample(empty,g,X,Y);
assert(map.allValid,'ipm_cr:GaugeTestWallMap','Every original wall node needs a native owner.');
end

function values=sample_wall(field,map)
values=reshape(map.core*field.core(:)+map.polar*field.polar(:),map.shape);
end

function [cp,pointer]=load_endpoint(folder)
assert(java.io.File(folder).isAbsolute(),'ipm_cr:GaugeTestCheckpoint','Use an absolute checkpoint folder.');
recorded=jsondecode(fileread(fullfile(folder,'latest.json')));
if strcmp(recorded.kind,'ipm_cr_explicit_accepted_pointer_v1')
    [cp,pointer]=ipm_cr.io.checkpoint('import',folder);
else
    [cp,pointer]=ipm_cr.io.checkpoint('load',folder);
end
assert(pointer.tau==200&&cp.state.scale.canonicalTime==200, ...
    'ipm_cr:GaugeTestCheckpoint','The frozen-state comparison requires an accepted tau200 state.');
end

function assert_unchanged(folder,pointer)
current=jsondecode(fileread(fullfile(folder,'latest.json')));
assert(current.sequence==pointer.sequence&&current.tau==pointer.tau&& ...
    strcmp(current.path,pointer.path)&&strcmp(current.sha256,pointer.sha256)&& ...
    strcmp(ipm_cr.io.sha(pointer_path(folder,pointer)),pointer.sha256), ...
    'ipm_cr:GaugeTestCheckpointChanged','Read-only inspection found a changed accepted checkpoint.');
end

function path=pointer_path(folder,pointer)
if strcmp(pointer.kind,'ipm_cr_accepted_pointer_v2'),path=fullfile(folder,pointer.path);
else,path=pointer.path;end
end

function info=compact_legacy(info)
if isfield(info,'baseRhs'),info=rmfield(info,'baseRhs');end
if isfield(info,'initialMeasurement'),info=rmfield(info,'initialMeasurement');end
if isfield(info,'fullMeasurement'),info=rmfield(info,'fullMeasurement');end
end

function d=subtract(a,b)
d=struct('core',a.core-b.core,'polar',a.polar-b.polar);
end

function inventory=source_inventory()
paths={mfilename('fullpath'),which('ipm_cr.op.prepare'),which('ipm_cr.op.apply'),which('ipm_cr.gauge.rates'), ...
    which('ipm_cr.rhs'),which('ipm_cr.mesh.fill'), ...
    which('ipm_cr.mesh.sample'),which('ipm_cr.mesh.diff'), ...
    which('ipm_cr.poisson.apply'),which('ipm_cr.gauge.seed'),which('ipm_cr.gauge.measure'), ...
    which('ipm_cr.io.sha'),which('ipm_cr.io.json')};
if ~endsWith(paths{1},'.m'),paths{1}=[paths{1},'.m'];end
privateDir=fullfile(fileparts(mfilename('fullpath')),'private');
paths{end+1}=fullfile(privateDir,'test_report.m');
inventory=test_inventory(paths);
end

function persist(folder,report)
test_report(folder,'gauge_acceleration',report);
end

function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
