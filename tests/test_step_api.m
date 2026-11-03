function report=test_step_api(outputDir,options)
%TEST_STEP_API Explicit context/time-control ownership, optional frozen oracle.
% OPTIONS.baselineFile explicitly requests a historical MAT comparison;
% optional baselineSHA256 binds that input. The default has no old-data need.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','step_api');end
if nargin<2,options=struct();end
assert(isstruct(options)&&isscalar(options)&& ...
    all(ismember(fieldnames(options),{'baselineFile','baselineSHA256'})), ...
    'ipm_cr:StepApiOptions','Use baselineFile and optional baselineSHA256.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'step_api.json')),'ipm_cr:StepApiOutput','Use a fresh output directory.');
timer=tic;source=ipm_cr.io.identity('source');
report=struct('kind','ipm_cr_narrow_step_api_v1','passed',false,'status','running', ...
    'baseline',struct('requested',isfield(options,'baselineFile'),'status','skipped', ...
    'file','','sha256',''),'sourceIdentity',source,'requestedChecksSatisfied',false, ...
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
fprintf('CR STEP API: independent controls, state identities and context isolation passed; baseline %s (%.3fs)\n', ...
    report.baseline.status,report.seconds);
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
