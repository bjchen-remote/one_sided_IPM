function report=test_runtime(outputDir)
%TEST_RUNTIME Short restart, portable identity and accepted-pointer checks.
% This suite advances only two short original-datum intervals; it is not a
% resolution study or a long-time trajectory qualification.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','runtime');end
outputDir=char(outputDir);
assert(java.io.File(outputDir).isAbsolute(),'ipm_cr:RuntimeTestOutput','Use an absolute output directory.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'runtime.json')),'ipm_cr:RuntimeTestFresh','Use a fresh output directory.');
timer=tic;source=ipm_cr.io.identity('source');
method=struct('geometry',struct('coreSpacing',.06, ...
    'grid',struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
    'solverOptions',struct('coreDegree',[48,48],'angularDegree',64,'polarDegree',16,'rOuter',1e10), ...
    'cfl',.3,'maxDt',.2,'timeIntegrator','rk4');
first=fullfile(outputDir,'original');moved=fullfile(outputDir,'moved');
options=struct('outputDir',first,'tauEnd',.003,'wallSeconds',120,'verbose',false,'explicitOptions',method);
[initial,state,ctx]=ipm_cr.run(options);
assert(strcmp(initial.status,'target_reached'),'ipm_cr:RuntimeTestInitial');
[cp,pointer]=ipm_cr.io.checkpoint('load',first);
assert(strcmp(pointer.path,'accepted_1.mat')||strcmp(pointer.path,'accepted_0.mat'),'ipm_cr:RuntimeTestRelative');
assert(~java.io.File(pointer.path).isAbsolute()&&cp.state.step==state.step,'ipm_cr:RuntimeTestPayload');

% Match the split endpoint and cold refresh exactly in the direct reference.
reference=ipm_cr.refresh(state,ctx);
while reference.scale.canonicalTime<.006
    [reference,info]=ipm_cr.step(reference,ctx,.006);
    assert(info.realizedMaximumStageCfl<=method.cfl*(1+64*eps),'ipm_cr:RuntimeTestCfl');
end
[ok,message]=copyfile(first,moved);assert(ok,'ipm_cr:RuntimeTestCopy','%s',message);
options.tauEnd=.006;
[resumed,a]=ipm_cr.run(options);
options.outputDir=moved;[relocated,b]=ipm_cr.run(options);
assert(resumed.resumed&&relocated.resumed&&strcmp(resumed.status,'target_reached')&& ...
    strcmp(relocated.status,'target_reached'),'ipm_cr:RuntimeTestResume');
restart=compare_state(reference,a);migration=compare_state(a,b);
assert(restart.maximumRelativeDifference<=2e-11&&migration.maximumRelativeDifference<=2e-11, ...
    'ipm_cr:RuntimeTestEquivalence','Restart or directory migration changed the short numerical trajectory.');

% A requested STOP must not replace an already-published accepted state.
cursor=fullfile(moved,'latest.json');before=fileread(cursor);
stop=fopen(fullfile(moved,'STOP'),'w');assert(stop>=0,'ipm_cr:RuntimeTestStop');fclose(stop);
options.tauEnd=.009;stopped=ipm_cr.run(options);
assert(strcmp(stopped.status,'stopped')&&strcmp(before,fileread(cursor)), ...
    'ipm_cr:RuntimeTestStop','STOP changed the durable accepted pointer.');

% Corrupt a disposable pointer, never the valid accepted-state fixture.
corrupt=fullfile(outputDir,'corrupt');[ok,message]=copyfile(first,corrupt);
assert(ok,'ipm_cr:RuntimeTestCopy','%s',message);
bad=jsondecode(fileread(fullfile(corrupt,'latest.json')));bad.sha256=repmat('0',1,64);
ipm_cr.io.json(fullfile(corrupt,'latest.json'),bad);
expect_failure(@()ipm_cr.io.checkpoint('load',corrupt),'ipm_cr:CheckpointIdentity');

% Reject an obsolete writer cursor before either currently accepted slot is touched.
stale=fullfile(outputDir,'stale_cursor');[ok,message]=copyfile(first,stale);
assert(ok,'ipm_cr:RuntimeTestCopy','%s',message);
[payload,oldPointer]=ipm_cr.io.checkpoint('load',stale);
ipm_cr.io.checkpoint('save',stale,payload,oldPointer);
before=fileread(fullfile(stale,'latest.json'));
expect_failure(@()ipm_cr.io.checkpoint('save',stale,payload,oldPointer),'ipm_cr:CheckpointStale');
assert(strcmp(before,fileread(fullfile(stale,'latest.json'))),'ipm_cr:RuntimeTestStale');
ipm_cr.io.checkpoint('load',stale);

% Relocating unchanged project inputs must preserve content identities.
clone=fullfile(outputDir,'project_copy');mkdir(clone);
[ok,message]=copyfile(fullfile(root,'src'),fullfile(clone,'src'));assert(ok,'ipm_cr:RuntimeTestCopy','%s',message);
[ok,message]=copyfile(fullfile(root,'data'),fullfile(clone,'data'));assert(ok,'ipm_cr:RuntimeTestCopy','%s',message);
assert(isequaln(source,ipm_cr.io.identity('source',clone))&& ...
    isequaln(ipm_cr.io.identity('authority'),ipm_cr.io.identity('authority',clone)), ...
    'ipm_cr:RuntimeTestProjectMove','Machine root leaked into a content identity.');

% Historical imports remain explicit and read-only, never a v2 continuation.
legacy=fullfile(outputDir,'legacy_fixture');mkdir(legacy);
checkpoint=cp;checkpoint.kind='ipm_cr_explicit_accepted_checkpoint_v1';
% Construct the historical schema only in this disposable v1 fixture.
checkpoint.state.rho=checkpoint.state.W;checkpoint.state=rmfield(checkpoint.state,'W');
if isfield(checkpoint.state,'eta')
    checkpoint.state.residual=checkpoint.state.eta;checkpoint.state=rmfield(checkpoint.state,'eta');
end
checkpoint.state.flow=checkpoint.state.rhs;checkpoint.state=rmfield(checkpoint.state,'rhs');
checkpoint.rho=checkpoint.state.rho;
path=fullfile(legacy,'accepted_0.mat');save(path,'checkpoint','-v7.3');
p=struct('kind','ipm_cr_explicit_accepted_pointer_v1','path',path,'generation',0, ...
    'sequence',checkpoint.state.step,'tau',checkpoint.state.scale.canonicalTime,'sha256',ipm_cr.io.sha(path));
ipm_cr.io.json(fullfile(legacy,'latest.json'),p);legacyBefore=fileread(fullfile(legacy,'latest.json'));
expect_failure(@()ipm_cr.io.checkpoint('load',legacy),'ipm_cr:CheckpointLegacy');
imported=ipm_cr.io.checkpoint('legacy',legacy);
assert(isequaln(imported,checkpoint)&&strcmp(legacyBefore,fileread(fullfile(legacy,'latest.json'))), ...
    'ipm_cr:RuntimeTestLegacy','Read-only legacy import changed historical data.');
normalized=ipm_cr.io.checkpoint('import',legacy);
assert(all(isfield(normalized.state,{'W','eta','rhs'}))&& ...
    ~any(isfield(normalized.state,{'rho','residual','flow'}))&& ...
    isequaln(normalized.state.W,checkpoint.state.rho)&& ...
    isequaln(normalized.state.eta,checkpoint.state.residual)&& ...
    isequaln(normalized.state.rhs,checkpoint.state.flow)&& ...
    isequaln(rmfield(normalized,'state'),rmfield(checkpoint,'state'))&& ...
    strcmp(legacyBefore,fileread(fullfile(legacy,'latest.json'))), ...
    'ipm_cr:RuntimeTestImport','Normalized read-only import changed numerical data or historical files.');
assert(isequaln(source,ipm_cr.io.identity('source')),'ipm_cr:RuntimeTestSources','Source changed during the runtime test.');
report=struct('kind','ipm_cr_runtime_transactions_v2','passed',true, ...
    'restart',restart,'outputMigration',migration,'relativePointer',true, ...
    'stopPreservedPointer',true,'corruptShaRejected',true,'staleWriterRejected',true,'projectIdentityPortable',true, ...
    'legacyImportExplicitReadOnly',true,'sourceUnchanged',true, ...
    'trajectoryAccuracyCertified',false,'tauEnd',.006,'seconds',toc(timer));
ipm_cr.io.json(fullfile(outputDir,'runtime.json'),report);
end

function metrics=compare_state(a,b)
fields={'W','scale'};if isfield(a,'eta'),fields{end+1}='eta';end
metrics=struct('maximumRelativeDifference',0,'sameStep',a.step==b.step);
assert(metrics.sameStep,'ipm_cr:RuntimeTestSteps','Compared trajectories took different accepted steps.');
for k=1:numel(fields)
    value=relative_difference(a.(fields{k}),b.(fields{k}));
    metrics.maximumRelativeDifference=max(metrics.maximumRelativeDifference,value);
end
end

function value=relative_difference(a,b)
assert(strcmp(class(a),class(b))&&isequal(size(a),size(b)),'ipm_cr:RuntimeTestShape');
value=0;
if isstruct(a)
    names=fieldnames(a);assert(isequal(names,fieldnames(b)),'ipm_cr:RuntimeTestFields');
    for k=1:numel(names),value=max(value,relative_difference(a.(names{k}),b.(names{k})));end
elseif isnumeric(a)
    assert(all(isfinite(a),'all')&&all(isfinite(b),'all'),'ipm_cr:RuntimeTestFinite');
    value=max(abs(a-b),[],'all')/max(1,max(abs(a),[],'all'));
else
    assert(isequaln(a,b),'ipm_cr:RuntimeTestValues');
end
end

function expect_failure(callback,identifier)
try
    callback();error('ipm_cr:RuntimeTestExpectedFailure','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:RuntimeTestWrongFailure', ...
        'Expected %s, received %s.',identifier,problem.identifier);
end
end
