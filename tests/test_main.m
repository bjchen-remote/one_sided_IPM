function report=test_main(outputDir)
%TEST_MAIN Short expression-initialized entrypoint run from an unrelated cwd.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'src'),fullfile(root,'tools'),'-begin');
if nargin<1,outputDir=fullfile(root,'results','test_main');end
outputDir=char(java.io.File(outputDir).getCanonicalPath());
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'main.json')),'ipm_cr:MainTestFresh','Use a fresh output directory.');
before=ipm_cr.io.identity('source');mainSha=ipm_cr.io.sha(fullfile(root,'main.m'));
old=pwd;restore=onCleanup(@()cd(old)); %#ok<NASGU>
foreign=fullfile(outputDir,'foreign 工作目录');mkdir(foreign);cd(foreign);
timer=tic;report=struct('kind','ipm_cr_main_expression_test_v1','passed',false, ...
    'status','running','sourceUnchanged',false,'trajectoryAccuracyCertified',false);
try
    [c,example]=main('config');
    assert(isa(example,'function_handle')&&isstruct(c)&&strcmp(c.format,'ipm_cr_case_v1')&& ...
        strcmp(pwd,foreign),'ipm_cr:MainTestConfig');
    folder=fullfile(outputDir,'表达式 accepted run with spaces');
    assert(~isfolder(folder),'ipm_cr:MainTestFresh','The test must own a new run folder.');
    base=java.io.File(root);target=java.io.File(folder);
    c.run.outputDir=char(base.toPath().relativize(target.toPath()).toString());
    c.run.tauEnd=.002;c.run.wallSeconds=120;c.run.verbose=false;c.run.maxSteps=1;
    c.run.checkpointEveryTau=.001;c.run.checkpointWallSeconds=60;
    c.method=struct('geometry',struct('coreSpacing',.06, ...
        'grid',struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
        'solverOptions',struct('coreDegree',[32,32],'angularDegree',40, ...
        'polarDegree',12,'requireAccuracy',true),'maxDt',.001,'cfl',.3);
    c.stop.steady.enabled=false;c.save.label='函数初值 UTF-8';
    rho0=@(x,y)(1+tanh(x)).*exp(-y.^2);
    [a,s,ctx]=main(c,rho0);
    assert(strcmp(pwd,foreign)&&strcmp(a.status,'maximum_steps')&&s.step==1&& ...
        strcmp(char(java.io.File(a.outputDir).getCanonicalPath()),folder), ...
        'ipm_cr:MainTestFirst','Main must retain cwd and stop at its first accepted step.');
    assert(strcmp(ctx.opts.background.kind,'none')&&~isfield(s,'eta'), ...
        'ipm_cr:MainTestRepresentation','A generic expression must use full W without Gaussian eta.');
    [cp,pointer]=ipm_cr.io.checkpoint('load',folder);
    verify(a,s,ctx,cp,pointer);
    initial=ctx.opts.initial;identity=ipm_cr.io.identity('authority',root,initial);
    cfg=ipm_cr.problem(initial);loaded=load(cfg.authorityFile,'frame');
    native=loaded.frame.native;g=ctx.mesh;expected=rho0(g.core.X,g.core.Y);
    mask=g.core.active;
    assert(isequaln(native.geometryKey,g.geometryKey)&& ...
        isequaln(native.core(mask),expected(mask)),'ipm_cr:MainTestSamples');
    [s0,own]=ipm_cr.init(ctx.opts);
    reconstructed=s0.W.core+s0.scale.densityOrigin;
    scale=max(1,max(abs(expected(mask))));
    reconstructionError=max(abs(reconstructed(mask)-expected(mask)));
    assert(reconstructionError<=64*eps(scale)&&~isfield(s0,'eta')&& ...
        strcmp(own.opts.background.kind,'none'),'ipm_cr:MainTestInitial');
    [~,~,refs]=ipm_cr.model.reference(cfg);
    [~,~,oldRefs]=ipm_cr.model.reference(ipm_cr.problem('data/initial.json'));
    assert(isequaln(refs,s.originalAuthority)&&refs.I0>0&&refs.P0>0, ...
        'ipm_cr:MainTestReference','I/P/C must be rebuilt from the actual expression.');
    assert(~isequaln(refs,oldRefs),'ipm_cr:MainTestReference','Do not reuse the original Gaussian I/P/C.');
    c.run.maxSteps=20;[b,s2,ctx2]=main(c,rho0);
    assert(strcmp(pwd,foreign)&&b.resumed&&strcmp(b.status,'target_reached')&& ...
        s2.scale.canonicalTime==.002&&s2.step>s.step&& ...
        strcmp(ctx2.opts.initial,initial),'ipm_cr:MainTestResume');
    assert(isequaln(identity,ipm_cr.io.identity('authority',root,initial)), ...
        'ipm_cr:MainTestReuse','The same materialized function and mesh must reuse immutable input.');
    [cp2,pointer2]=ipm_cr.io.checkpoint('load',folder);
    verify(b,s2,ctx2,cp2,pointer2);
    raw=fileread(fullfile(folder,'summary.json'));saved=jsondecode(raw);
    assert(strcmp(saved.status,b.status)&& ...
        strcmp(char(java.io.File(saved.outputDir).getCanonicalPath()),folder)&& ...
        strcmp(saved.runControls.save.label,c.save.label)&&contains(raw,'函数初值'), ...
        'ipm_cr:MainTestUtf8','Summary JSON must preserve Unicode paths and labels.');
    report.sourceUnchanged=isequaln(before,ipm_cr.io.identity('source'))&& ...
        strcmp(mainSha,ipm_cr.io.sha(fullfile(root,'main.m')));
    assert(report.sourceUnchanged,'ipm_cr:MainTestSource','Maintained input sources changed.');
    report.passed=true;report.status='passed_short_expression_entrypoint';
    report.foreignCwdRetained=true;report.relativeOutputWithSpacesAndUnicode=true;
    report.materializedInitial=initial;report.initialIdentity=identity;
    report.activeCoreInitialReconstructionInfinity=reconstructionError;
    report.reference=refs;report.finalTau=s2.scale.canonicalTime;report.finalSteps=s2.step;
    report.maximumStageCfl=b.maximumAcceptedStageCfl;
    report.expressionManifestReused=true;report.packageAndCheckpointVerified=true;
    report.utf8SummaryReadback=true;
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'main.json'),report);rethrow(problem);
end
report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'main.json'),report);
end

function verify(summary,s,ctx,cp,pointer)
assert(strcmp(summary.kind,'ipm_cr_run_v2')&&isa(ctx,'ipm_cr.context')&& ...
    pointer.sequence==s.step&&pointer.tau==s.scale.canonicalTime&&cp.acceptedOnly&& ...
    isequaln(cp.state.W,s.W)&&isequaln(cp.state.scale,s.scale)&& ...
    isequaln(cp.explicitOptions,ctx.opts),'ipm_cr:MainTestCheckpoint');
assert(ctx.opts.solverOptions.requireAccuracy&&s.rhs.poisson.accepted&& ...
    summary.maximumAcceptedStageCfl<=.3*(1+64*eps)&&ctx.opts.cfl==.3, ...
    'ipm_cr:MainTestCfl','Main must retain original Poisson acceptance and all-stage CFL.');
if isfield(s.rhs.poisson,'linearSuperposition')
    assert(~s.rhs.poisson.linearSuperposition.accuracyGatesRelaxed,'ipm_cr:MainTestGate');
end
package=summary.dataPackage;
file=fullfile(package.folder,package.fieldFile);
assert(strcmp(ipm_cr.io.sha(file),package.fieldSHA256)&& ...
    strcmp(ipm_cr.io.sha(fullfile(package.folder,package.metadataFile)),package.metadataSHA256), ...
    'ipm_cr:MainTestPackage');
loaded=load(file,'data');data=loaded.data;
assert(data.step==s.step&&isequaln(data.scale,s.scale)&&isequaln(data.method,ctx.opts), ...
    'ipm_cr:MainTestPackage');
for patch={'core','polar'}
    key=patch{1};
    assert(isequaln(data.fields.(key).W,s.W.(key))&& ...
        isequaln(data.fields.(key).F,s.rhs.F.(key))&& ...
        isequaln(data.fields.(key).Q,s.rhs.Q.(key))&& ...
        isequaln(cp.state.rhs.F.(key),s.rhs.F.(key))&& ...
        isequaln(cp.state.rhs.Q.(key),s.rhs.Q.(key)), ...
        'ipm_cr:MainTestPackage','Native fields must agree with the returned and accepted state.');
end
end
