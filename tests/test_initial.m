function report=test_initial(outputDir)
%TEST_INITIAL New Gaussian authorities, provenance and initialization contract.
% Generated cases are retained under data/cases with unique names. This is
% an initialization/identity test, not a trajectory accuracy qualification.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'),fullfile(root,'tools'));
if nargin<1||isempty(outputDir)
    outputDir=fullfile(root,'results',['initial_',char(java.util.UUID.randomUUID())]);
end
outputDir=char(outputDir);
assert(startsWith(outputDir,filesep),'ipm_cr:InitialTestOutput','Use an absolute output directory.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'initial.json')),'ipm_cr:InitialTestFresh','Choose a fresh test output directory.');
timer=tic;source=ipm_cr.io.identity('source');base=ipm_cr.problem();
defaultIdentity=ipm_cr.io.identity('authority');
assert(isequaln(defaultIdentity,ipm_cr.io.identity('authority',root,'data/initial.json')), ...
    'ipm_cr:InitialTestDefault','Default and explicit original manifests must select the same identity.');
name=['test_',char(java.util.UUID.randomUUID())];
parameters=struct('amplitude',.63,'centerX',.1,'widthX',.85,'widthY',.9,'skewWeight',.25);
[initial,generated]=make_initial(name,parameters);cfg=ipm_cr.problem(initial);
customIdentity=ipm_cr.io.identity('authority',root,initial);
assert(strcmp(cfg.manifestFile,initial)&&strcmp(cfg.manifestSha,customIdentity.manifestSha256)&& ...
    strcmp(customIdentity.manifestSha256,ipm_cr.io.sha(fullfile(root,initial)))&& ...
    strcmp(cfg.authoritySha,customIdentity.sha256)&&~isequaln(defaultIdentity,customIdentity), ...
    'ipm_cr:InitialTestIdentity','A new Gaussian manifest must have a distinct bound identity.');
assert(isequal([cfg.physics.gaussianAmplitude,cfg.physics.gaussianCenterX, ...
    cfg.physics.gaussianWidthX,cfg.physics.gaussianWidthY,cfg.physics.gaussianSkewWeight], ...
    [.63,.1,.85,.9,.25]),'ipm_cr:InitialTestParameters');
old=load(base.authorityFile,'frame');data=load(cfg.authorityFile,'frame');frame=data.frame;
assert(isequaln(frame.x,old.frame.x)&&isequaln(frame.y,old.frame.y)&& ...
    isequal([frame.step,frame.tau,frame.t,frame.Xshift],[0,0,0,0])&& ...
    isequal([frame.Cx,frame.Cy,frame.Comega],[1,1,1]), ...
    'ipm_cr:InitialTestFrame','New authorities must retain the fixed axes and exact t0 frame.');
[wall,scale,refs,guess]=ipm_cr.model.reference(cfg);[~,~,oldRefs]=ipm_cr.model.reference(base);
assert(~isequaln(refs,oldRefs)&&isequaln(refs,generated.references)&& ...
    all(isfinite([refs.I0,refs.P0,refs.C0,guess]))&&refs.I0>0&&refs.P0>0&& ...
    scale.canonicalTime==0&&scale.physicalTime==0&&scale.physicalTimeLow==0, ...
    'ipm_cr:InitialTestReference','New I/P/C must be rebuilt from the selected original measurement axes.');
manifest=jsondecode(fileread(fullfile(root,initial)));caseFolder=fileparts(fullfile(root,initial));
manifestSHA=ipm_cr.io.sha(fullfile(root,initial));matSHA=ipm_cr.io.sha(cfg.authorityFile);
expect(@()make_initial(name,parameters),'ipm_cr:InitialExists');
assert(strcmp(manifestSHA,ipm_cr.io.sha(fullfile(root,initial)))&& ...
    strcmp(matSHA,ipm_cr.io.sha(cfg.authorityFile)), ...
    'ipm_cr:InitialTestOverwrite','Refusing an existing name must leave its files unchanged.');
expect(@()make_initial([name,'x'],struct('eta',zeros(2))),'ipm_cr:InitialParameters');
expect(@()make_initial([name,'x'],struct('widthX',.8,'gaussianWidthX',.8)), ...
    'ipm_cr:InitialParameters');
expect(@()ipm_cr.problem('../data/initial.json'),'ipm_cr:ProblemPath');
expect(@()ipm_cr.problem(fullfile(root,'data','initial.json')),'ipm_cr:ProblemPath');
expect(@()ipm_cr.problem('data\initial.json'),'ipm_cr:ProblemPath');
expect(@()ipm_cr.problem('data//initial.json'),'ipm_cr:ProblemPath');
expect(@()ipm_cr.io.identity('authority',root,'../data/initial.json'),'ipm_cr:AuthorityPath');

% Corrupt only new disposable manifests, never the generated valid input.
bad=manifest;bad.sha256=repmat('0',1,64);
badInitial=fixture('bad_sha.json',bad);
expect(@()ipm_cr.problem(badInitial),'ipm_cr:ProblemIdentity');
expect(@()ipm_cr.io.identity('authority',root,badInitial),'ipm_cr:AuthoritySha');
bad=manifest;bad.file='../escape.mat';badInitial=fixture('bad_path.json',bad);
expect(@()ipm_cr.problem(badInitial),'ipm_cr:ProblemPath');
expect(@()ipm_cr.io.identity('authority',root,badInitial),'ipm_cr:AuthorityPath');
bad=manifest;bad.physics.gaussianWidthX=.97*bad.physics.gaussianWidthX;
badCfg=ipm_cr.problem(fixture('bad_parameters.json',bad));
expect(@()ipm_cr.model.reference(badCfg),'polar:T0AnalyticProvenance');
late=frame;late.tau=1;frame=late; %#ok<NASGU>
latePath=fullfile(caseFolder,'late.mat');save(latePath,'frame','-v7.3');
bad=manifest;bad.file=['data/cases/',name,'/late.mat'];bad.sha256=ipm_cr.io.sha(latePath);
badCfg=ipm_cr.problem(fixture('late.json',bad));
expect(@()ipm_cr.model.reference(badCfg),'polar:T0AuthorityClock');
bad=manifest;bad.note='disposable metadata identity test';metadataInitial=fixture('metadata.json',bad);
first=ipm_cr.io.identity('authority',root,metadataInitial);
bad.note='changed disposable metadata';fixture('metadata.json',bad);
second=ipm_cr.io.identity('authority',root,metadataInitial);
assert(strcmp(first.file,second.file)&&strcmp(first.sha256,second.sha256)&& ...
    ~strcmp(first.manifestSha256,second.manifestSha256)&&~isequaln(first,second), ...
    'ipm_cr:InitialTestManifest','Manifest metadata must be bound even when MAT bytes are unchanged.');

clone=fullfile(outputDir,'project_copy');mkdir(clone);
for path={initial,manifest.file}
    target=fullfile(clone,path{1});parent=fileparts(target);if ~isfolder(parent),mkdir(parent);end
    [ok,message]=copyfile(fullfile(root,path{1}),target);assert(ok,'ipm_cr:InitialTestCopy','%s',message);
end
assert(isequaln(customIdentity,ipm_cr.io.identity('authority',clone,initial)), ...
    'ipm_cr:InitialTestPortable','Machine paths must not enter a Gaussian authority identity.');
method=ipm_cr.config(struct('initial',initial));
expected=struct('kind','gaussian_primitive','amplitude',parameters.amplitude, ...
    'centerX',parameters.centerX,'widthX',parameters.widthX,'widthY',parameters.widthY, ...
    'skewWeight',parameters.skewWeight,'tailCorrection',true);
assert(strcmp(method.initial,initial)&&isequaln(method.background,expected), ...
    'ipm_cr:InitialTestBackground','The default background must use the selected Gaussian parameters.');
rejectFolder=fullfile(outputDir,'wrong_authority');mkdir(rejectFolder);
payload=struct('sourceInventory',source,'authority',defaultIdentity, ...
    'explicitOptions',struct('initial',initial));
expect(@()ipm_cr.io.checkpoint('save',rejectFolder,payload,[]),'ipm_cr:CheckpointAuthorityChanged');
assert(~isfile(fullfile(rejectFolder,'latest.json')),'ipm_cr:InitialTestCheckpoint', ...
    'An old authority must not be published under a new initial manifest.');
assert(isequaln(source,ipm_cr.io.identity('source'))&& ...
    isequaln(defaultIdentity,ipm_cr.io.identity('authority')), ...
    'ipm_cr:InitialTestOriginal','Generating a new case must not change source or the original authority.');
report=struct('kind','Gaussian_initial_authority_release_contract_v1','passed',true, ...
    'initial',initial,'authority',customIdentity,'originalAuthority',defaultIdentity, ...
    'references',refs,'fixedMeasurementAxes',true,'newProblemIdentity',true, ...
    'manifestAndMatSHABound',true,'identityPortable',true,'overwriteRejected',true, ...
    'escapedPathsRejected',true,'lateFrameRejected',true,'parameterMismatchRejected',true, ...
    'oldCheckpointAuthorityRejected',true,'defaultBackgroundMatchesSelectedGaussian',true, ...
    'initializationIntegrationNotRun',true,'nativeRequired',false, ...
    'arbitraryEtaAndCallbacksExposed',false,'timeAdvanced',0, ...
    'trajectoryAccuracyCertified',false,'seconds',toc(timer));
ipm_cr.io.json(fullfile(outputDir,'initial.json'),report);
fprintf('CR INITIAL: new Gaussian authority, fixed I/P/C axes and identity contracts passed (%.3fs)\n',report.seconds);

    function relative=fixture(file,value)
        ipm_cr.io.json(fullfile(caseFolder,file),value);
        relative=['data/cases/',name,'/',file];
    end
end

function expect(callback,identifier)
try
    callback();error('ipm_cr:InitialTestExpected','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:InitialTestWrongFailure', ...
        'Expected %s, got %s: %s',identifier,problem.identifier,problem.message);
end
end
