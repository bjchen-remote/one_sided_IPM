function report=test_step_api(outputDir,options)
%TEST_STEP_API Explicit context/time-control ownership, optional frozen oracle.
% OPTIONS.baselineFile explicitly requests a historical MAT comparison;
% optional baselineSHA256 binds that input. The default has no old-data need.
% The four failureState/geometryCheckpoint file/SHA options together request
% the real secant-zero regression on frozen checkpoint spatial matrices.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','step_api');end
if nargin<2,options=struct();end
assert(isstruct(options)&&isscalar(options)&& ...
    all(ismember(fieldnames(options),{'baselineFile','baselineSHA256', ...
    'failureStateFile','failureStateSHA256','geometryCheckpointFile','geometryCheckpointSHA256'})), ...
    'ipm_cr:StepApiOptions','Use baseline options or all four failure-state/geometry-checkpoint file/SHA options.');
failureFields={'failureStateFile','failureStateSHA256','geometryCheckpointFile','geometryCheckpointSHA256'};
failurePresent=isfield(options,failureFields);
assert(~any(failurePresent)||all(failurePresent),'ipm_cr:StepApiOptions', ...
    'The failure-state and geometry-checkpoint files and both SHA256 values must be supplied together.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'step_api.json')),'ipm_cr:StepApiOutput','Use a fresh output directory.');
timer=tic;source=ipm_cr.io.identity('source');
report=struct('kind','ipm_cr_narrow_step_api_v1','passed',false,'status','running', ...
    'baseline',struct('requested',isfield(options,'baselineFile'),'status','skipped', ...
    'file','','sha256',''),'sourceIdentity',source,'requestedChecksSatisfied',false, ...
    'failureRegression',struct('requested',all(failurePresent),'performed',false, ...
    'status','skipped','passed',false,'reason','No external real-failure fixture requested.'), ...
    'frozenCases',struct([]),'controls',struct(),'identities',struct(), ...
    'sourceUnchanged',false,'trajectoryAccuracyCertified',false, ...
    'scope','Normal prepared-context calls and direct time controls; historical bitwise evidence only when explicitly requested. No order or long-time claim.');
try
    assert(nargin('ipm_cr.time.step')==4&&nargin('ipm_cr.step')==4,'ipm_cr:StepApiSignature');
    baseline=struct('cases',struct([]));
    method=struct('geometry',struct('coreSpacing',.06, ...
        'grid',struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
        'solverOptions',struct('coreDegree',[48,48],'angularDegree',64,'polarDegree',16,'rOuter',1e10), ...
        'cfl',.3,'maxDt',.001,'timeIntegrator','rk4');
    if report.baseline.requested
        report.baseline.status='requested';report.baseline.file=char(options.baselineFile);
        assert(~isempty(report.baseline.file)&&isfile(report.baseline.file), ...
            'ipm_cr:StepApiBaseline','The explicitly requested baseline file is missing.');
        report.baseline.sha256=ipm_cr.io.sha(report.baseline.file);
        if isfield(options,'baselineSHA256')
            assert(strcmp(report.baseline.sha256,char(options.baselineSHA256)), ...
                'ipm_cr:StepApiBaseline','The explicitly bound baseline SHA changed.');
        end
        data=load(report.baseline.file,'baseline');baseline=data.baseline;
        assert(numel(baseline.cases)==2&&all(isfield(baseline.cases,{'initial','final','stepNumerics'})), ...
            'ipm_cr:StepApiBaseline','Use the complete two-case frozen trajectory.');
        cfg=ipm_cr.problem();
        assert(strcmp(ipm_cr.io.sha(cfg.authorityFile),baseline.authoritySHA),'ipm_cr:StepApiAuthority');
        method=baseline.cases(1).method;
    else
        assert(~isfield(options,'baselineSHA256'),'ipm_cr:StepApiOptions','baselineSHA256 needs baselineFile.');
    end
    for j=1:numel(baseline.cases)
        old=baseline.cases(j);[s,ctx]=ipm_cr.init(old.method);
        exact(snapshot(s),old.initial,'initial');
        [v,d]=ipm_cr.diag.compare(s,ctx.mesh,ctx.opts,baseline.queryCloud);
        exact(v,old.initialSample,'initialSample');exact(d,old.initialDiagnostic,'initialDiagnostic');
        steps=0;
        while s.scale.canonicalTime<baseline.tauEnd
            [s,info]=ipm_cr.step(s,ctx,baseline.tauEnd);
            steps=steps+1;
            assert(steps<=numel(old.stepNumerics),'ipm_cr:StepApiSteps','More steps than the frozen trajectory.');
            exact(rmfield(info,'seconds'),old.stepNumerics(steps),'stepNumerics');
        end
        assert(steps==numel(old.stepNumerics),'ipm_cr:StepApiSteps','The frozen step count changed.');
        exact(snapshot(s),old.final,'final');
        [v,d]=ipm_cr.diag.compare(s,ctx.mesh,ctx.opts,baseline.queryCloud);
        exact(v,old.finalSample,'finalSample');exact(d,old.finalDiagnostic,'finalDiagnostic');
        row=struct('caseIndex',j,'steps',steps,'tau',s.scale.canonicalTime, ...
            'initialAndFinalBitwiseIdentical',true,'allStepNumericsBitwiseIdentical',true, ...
            'commonCloudBitwiseIdentical',true);
        if isempty(report.frozenCases),report.frozenCases=row;else,report.frozenCases(end+1)=row;end %#ok<AGROW>
    end
    if report.baseline.requested,report.baseline.status='passed_exact_frozen_steps';end
    % Identity is inspected by this test, not enforced in the time kernel.
    % Distinct time policies must leave the mathematical input key unchanged.
    [a,ca]=ipm_cr.init(method);
    controls=struct('cfl',.2,'maxDt',.0005,'timeIntegrator','ssprk3');
    names=fieldnames(controls);
    for k=1:numel(names),method.(names{k})=controls.(names{k});end
    [b,cb]=ipm_cr.init(method);
    exact(ca.key,cb.key,'timeIndependentKey');exact(snapshot(a),snapshot(b),'timeIndependentInitialRhs');
    assert(~strcmp(ca.id,cb.id)&&strcmp(a.ctxId,ca.id)&&strcmp(b.ctxId,cb.id)&& ...
        isequaln(a.ctxKey,ca.key)&&isequaln(b.ctxKey,cb.key),'ipm_cr:StepApiInputIdentity');
    savedOptions=ca.opts;savedKey=ca.key;
    otherBefore=cb.stats();
    [nextA,ia]=ipm_cr.step(a,ca,.001,controls);
    exact(cb.stats(),otherBefore,'otherContextAfterOverrideStep');
    before=ca.stats();
    [nextB,ib]=ipm_cr.time.step(b,cb,.001);
    exact(ca.stats(),before,'otherContextAfterPreparedStep');
    exact(snapshot(nextA),snapshot(nextB),'controlledStep');
    exact(rmfield(ia,'seconds'),rmfield(ib,'seconds'),'controlledStepNumerics');
    exact(ca.opts,savedOptions,'contextOptions');exact(ca.key,savedKey,'contextKey');
    assert(strcmp(nextA.ctxId,ca.id)&&strcmp(nextB.ctxId,cb.id)&& ...
        isequaln(nextA.ctxKey,ca.key)&&isequaln(nextB.ctxKey,cb.key),'ipm_cr:StepApiOutputIdentity');
    expectedDt=min([controls.maxDt,.95*controls.cfl/a.rhs.rate,.001-a.scale.canonicalTime]);
    assert(ia.dt==expectedDt&&strcmp(ia.timeIntegrator,'ssprk3')&& ...
        ia.cflStepFormulaApplied&&isequal(ia.stageCfl,ia.dt*ia.stageRates)&& ...
        ia.realizedMaximumStageCfl==max(ia.stageCfl), ...
        'ipm_cr:StepApiControls','The independent policy or passive stage-CFL measurements changed.');
    report.controls=struct('keyIndependent',true,'initialRhsBitwiseIdentical',true, ...
        'overrideEqualsPreparedPolicyBitwise',true,'contextUnchanged',true, ...
        'initialCflFormulaExact',true,'stageCflMeasured',true, ...
        'dt',ia.dt,'maximumStageCfl',ia.realizedMaximumStageCfl,'timeIntegrator',ia.timeIntegrator);
    report.identities=struct('distinctContextIds',true,'inputKeysMatch',true,'outputKeysMatch',true, ...
        'mathematicalKeyIndependentOfTimePolicy',true,'contextsIsolated',true);
    if report.failureRegression.requested
        [report.failureRegression,failureProblem]=failure_regression(options);
        if ~isempty(failureProblem),rethrow(failureProblem);end
    end
    report.sourceUnchanged=isequaln(source,ipm_cr.io.identity('source'));
    if report.baseline.requested
        report.sourceUnchanged=report.sourceUnchanged&& ...
            strcmp(ipm_cr.io.sha(report.baseline.file),report.baseline.sha256);
    end
    assert(report.sourceUnchanged,'ipm_cr:StepApiSource','Source or the frozen baseline changed.');
    report.requestedChecksSatisfied=true;report.passed=true;report.status='passed_step_api';
catch problem
    if report.baseline.requested&&~strcmp(report.baseline.status,'passed_exact_frozen_steps')
        report.baseline.status='failed_requested';
    end
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'step_api.json'),report);rethrow(problem);
end
report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'step_api.json'),report);
fprintf('CR STEP API: independent controls, state identities and context isolation passed; baseline %s; failure regression %s (%.3fs)\n', ...
    report.baseline.status,report.failureRegression.status,report.seconds);
end

function [r,problem]=failure_regression(options)
timer=tic;problem=[];
r=struct('requested',true,'performed',false,'status','requested','passed',false, ...
    'scope','Real secant-zero fixture, frozen checkpoint spatial matrices and frozen initial RHS, candidate time kernel only; not a new active-only-D7 spatial trajectory.', ...
    'failureStateFile',char(options.failureStateFile),'failureStateSHA256',char(options.failureStateSHA256), ...
    'geometryCheckpointFile',char(options.geometryCheckpointFile), ...
    'geometryCheckpointSHA256',char(options.geometryCheckpointSHA256), ...
    'targetTau',12.5,'expectedFittedKnot',263,'completedSteps',0,'actualTau',NaN, ...
    'stepInformation',struct([]),'inputsUnchanged',false);
try
    assert(isfile(r.failureStateFile)&&isfile(r.geometryCheckpointFile), ...
        'ipm_cr:StepApiFailureFixture','An explicitly requested real-failure input is missing.');
    beforeFailure=ipm_cr.io.sha(r.failureStateFile);beforeGeometry=ipm_cr.io.sha(r.geometryCheckpointFile);
    assert(strcmp(beforeFailure,r.failureStateSHA256)&&strcmp(beforeGeometry,r.geometryCheckpointSHA256), ...
        'ipm_cr:StepApiFailureFixture','An explicitly bound real-failure input SHA differs.');
    failed=load(r.failureStateFile,'lastAcceptedInMemory');
    geometry=load(r.geometryCheckpointFile,'checkpoint');
    assert(isfield(failed,'lastAcceptedInMemory')&&isfield(geometry,'checkpoint'), ...
        'ipm_cr:StepApiFailureFixture','Expected lastAcceptedInMemory and checkpoint MAT variables.');
    initial=failed.lastAcceptedInMemory;cp=geometry.checkpoint;
    r.initialContextKeysMatch=isequaln(initial.ctxKey,cp.state.ctxKey);
    assert(r.initialContextKeysMatch,'ipm_cr:StepApiFailureFixture', ...
        'The in-memory state and geometry checkpoint must have identical frozen context keys.');
    assert(initial.scale.canonicalTime==12.479430743775856&&initial.step==3218, ...
        'ipm_cr:StepApiFailureFixture','Expected the actual tau12.479430743775856 step3218 failure state.');
    r.initialTau=initial.scale.canonicalTime;r.initialStep=initial.step;r.actualTau=r.initialTau;
    r.initialScale=initial.scale;r.frozenMethod=cp.explicitOptions;r.frozenSource=cp.sourceInventory;
    % Deliberately use the saved geometry matrices, not current mesh.build:
    % this isolates the time fix from the later spatial derivative change.
    ctx=ipm_cr.context(cp.geometry,cp.wallOps,cp.explicitOptions);
    r.candidateContextIdentity=ctx.identity;s=initial;s.ctxId=ctx.id;s.ctxKey=ctx.key;
    exact(rmfield(s,{'ctxId','ctxKey'}),rmfield(initial,{'ctxId','ctxKey'}),'failureFrozenInitialStateAndRhs');
    r.initialFrozenRhsPreserved=true;r.initialRefreshPerformed=false;r.performed=true;
    savedOptions=ctx.opts;savedKey=ctx.key;
    previous=s;[s,info]=ipm_cr.step(s,ctx,r.targetTau);
    r.firstInformation=info;r.firstScale=s.scale;r.completedSteps=1;r.actualTau=s.scale.canonicalTime;
    r.stepInformation=info;
    failure_step_checks(previous,s,ctx,info,savedOptions,savedKey);
    assert(info.timePartition.fitted&&info.timePartition.knotIndex==r.expectedFittedKnot, ...
        'ipm_cr:StepApiFailureKnot','The real failure fixture must fit original spline knot263.');
    while s.scale.canonicalTime<r.targetTau
        assert(r.completedSteps<16,'ipm_cr:StepApiFailureSteps', ...
            'The bounded real-failure regression exceeded sixteen completed short steps.');
        previous=s;[s,info]=ipm_cr.step(s,ctx,r.targetTau);
        r.completedSteps=r.completedSteps+1;r.actualTau=s.scale.canonicalTime;
        r.stepInformation(end+1)=info;
        failure_step_checks(previous,s,ctx,info,savedOptions,savedKey);
    end
    r.finalScale=s.scale;r.finalStep=s.step;r.maximumStageCfl=max([r.stepInformation.realizedMaximumStageCfl]);
    r.finalFieldInfinity=struct('F',max([abs(s.rhs.F.core(:));abs(s.rhs.F.polar(:))]), ...
        'Q',max([abs(s.rhs.Q.core(:));abs(s.rhs.Q.polar(:))]));
    r.failureStateSHA256After=ipm_cr.io.sha(r.failureStateFile);
    r.geometryCheckpointSHA256After=ipm_cr.io.sha(r.geometryCheckpointFile);
    r.inputsUnchanged=strcmp(beforeFailure,r.failureStateSHA256After)&&strcmp(beforeGeometry,r.geometryCheckpointSHA256After);
    assert(s.scale.canonicalTime==r.targetTau&&r.inputsUnchanged, ...
        'ipm_cr:StepApiFailureIdentity','The real-failure endpoint or bound inputs changed.');
    r.status='passed_real_secant_zero_frozen_space';r.passed=true;
catch failure
    r.status='failed_requested';r.errorIdentifier=failure.identifier;
    r.errorMessage=getReport(failure,'extended','hyperlinks','off');problem=failure;
end
r.seconds=toc(timer);
end

function failure_step_checks(previous,s,ctx,info,savedOptions,savedKey)
assert(isfinite(info.dt)&&info.dt>0&&info.dt<=info.initialProposedDt&& ...
    s.scale.canonicalTime>previous.scale.canonicalTime&&s.step==previous.step+1&& ...
    strcmp(info.timeIntegrator,'rk4')&&info.cflStepFormulaApplied&& ...
    isequal(info.stageCfl,info.dt*info.stageRates)&& ...
    all(isfinite(info.stageRates))&&isfinite(info.realizedMaximumStageCfl)&& ...
    info.realizedMaximumStageCfl<=.3,'ipm_cr:StepApiFailureFiniteStep');
for patch={'core','polar'}
    key=patch{1};
    assert(all(isfinite(s.rhs.F.(key)),'all')&&all(isfinite(s.rhs.Q.(key)),'all'), ...
        'ipm_cr:StepApiFailureFiniteFields','Complete F and Q must remain finite.');
end
exact(fieldnames(s.scale),fieldnames(previous.scale),'failureScaleApi');
names=fieldnames(s.scale);
for k=1:numel(names)
    value=s.scale.(names{k});
    if isnumeric(value),assert(all(isfinite(value),'all'),'ipm_cr:StepApiFailureFiniteScale');end
end
assert(info.physicalClockUpdate.wordCount==4&& ...
    isequal(clock_words(previous.scale),info.physicalClockUpdate.beforeWords)&& ...
    isequal(clock_words(s.scale),info.physicalClockUpdate.afterWords), ...
    'ipm_cr:StepApiFailureClockWords','Every physical-clock word must follow the public step API.');
exact(ctx.opts,savedOptions,'failureContextOptions');exact(ctx.key,savedKey,'failureContextKey');
assert(strcmp(s.ctxId,ctx.id)&&isequaln(s.ctxKey,ctx.key)&& ...
    isequaln(s.originalAuthority,previous.originalAuthority), ...
    'ipm_cr:StepApiFailureOutputIdentity');
end

function words=clock_words(scale)
tail=zeros(1,2);if isfield(scale,'physicalTimeTail'),tail=reshape(scale.physicalTimeTail,1,[]);end
words=[scale.physicalTime,scale.physicalTimeLow,tail];
end

function out=snapshot(s)
f=s.rhs;
out=struct('rho',s.W,'residual',s.eta,'scale',s.scale,'step',s.step, ...
    'flow',struct('q',f.q,'F',f.F,'Q',f.Q,'residualRhs',f.residualRhs, ...
    'physical',struct('core',f.physical.core,'polar',f.physical.polar), ...
    'rates',[f.c_l,f.c_omega,f.c_r,f.delta,f.Uref,f.transportOrigin]));
end

function exact(a,b,label)
path=mismatch(a,b,label);
assert(isempty(path),'ipm_cr:StepApiBitwise','Frozen numerical regression differs at %s.',path);
end

function path=mismatch(a,b,label)
path='';
if ~strcmp(class(a),class(b))||~isequal(size(a),size(b))||issparse(a)~=issparse(b)
    path=[label,' (class, shape or storage)'];return;
end
if isstruct(a)
    names=fieldnames(a);
    if ~isequal(names,fieldnames(b)),path=[label,' (fields)'];return;end
    for j=1:numel(a)
        for k=1:numel(names)
            name=names{k};path=mismatch(a(j).(name),b(j).(name),sprintf('%s(%d).%s',label,j,name));
            if ~isempty(path),return;end
        end
    end
elseif iscell(a)
    for j=1:numel(a)
        path=mismatch(a{j},b{j},sprintf('%s{%d}',label,j));if ~isempty(path),return;end
    end
elseif isnumeric(a)
    if issparse(a)
        [ar,ac,av]=find(a);[br,bc,bv]=find(b);
        if ~isequal(ar,br)||~isequal(ac,bc),path=[label,' (sparse indices)'];return;end
        path=mismatch(av,bv,[label,' (sparse values)']);
    elseif isreal(a)~=isreal(b)|| ...
            ~isequal(typecast(real(a(:)),'uint8'),typecast(real(b(:)),'uint8'))|| ...
            ~isequal(typecast(imag(a(:)),'uint8'),typecast(imag(b(:)),'uint8'))
        path=[label,' (numeric bytes)'];
    end
elseif ~isequaln(a,b)
    path=[label,' (value)'];
end
end
