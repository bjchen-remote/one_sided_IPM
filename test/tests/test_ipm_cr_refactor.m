function report=test_ipm_cr_refactor(outputDir,baselineMat)
%TEST_IPM_CR_REFACTOR Exact regression against an independently frozen run.
% BASELINEMAT contains baseline.cases (redistribution off/on), queryCloud,
% tauEnd, baselineRuntime and authoritySHA. Each case stores the complete
% initialization method, geometry key, wall operators, initial/final numeric
% snapshots, and the old query samples/diagnostics. Generate this baseline
% in a separate MATLAB process using its original frozen runtime.
% Only numerical state and observations are compared; solver cache objects
% and construction/timing metadata are deliberately outside the snapshot.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'time_evolution'));
assert(nargin==2,'ipm_cr:RefactorArguments','Supply a fresh absolute output directory and baseline MAT file.');
outputDir=char(outputDir);baselineMat=char(baselineMat);
assert(startsWith(outputDir,filesep)&&startsWith(baselineMat,filesep)&&isfile(baselineMat), ...
    'ipm_cr:RefactorPaths','Use absolute paths and an existing baseline file.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'refactor.json'))&&~isfile(fullfile(outputDir,'refactor.mat')), ...
    'ipm_cr:RefactorPreserve','Use a fresh output directory.');
timer=tic;
report=struct('kind','exact_ipm_cr_maintenance_regression_v2','status','initializing', ...
    'baselineMat',baselineMat,'baselineSHA',ipm_c.output.sha256(baselineMat), ...
    'sourceInventory',source_inventory(root),'sourceHashesUnchanged',false, ...
    'cases',struct([]),'passed',false, ...
    'comparisonRule','Numeric values have identical class, dimensions and bytes; sparse operators compare FIND indices and value bytes.', ...
    'scope','Two short original-datum RK4 trajectories, complete common-cloud readout, explicit-stage CFL and constant annihilation.', ...
    'withdrawnCacheExperiment',struct('status','withdrawn_from_maintained_source', ...
        'reason','The measured derivative microbenchmark did not show a speed benefit, so its additional cache and invalidation code were removed.', ...
        'candidateEvidence',fullfile(root,'test','results','cleanup_20261022','regression','refactor.json'), ...
        'candidateEvidenceCertifiesCurrentSource',false), ...
    'limitations',{{'This maintenance regression preserves the frozen algorithm; it is not a new spatial-order or long-time accuracy certificate.', ...
        'No speed optimization is certified. The withdrawn cache experiment is historical evidence only.'}}, ...
    'matlabVersion',version,'computer',computer,'elapsedSeconds',0);
persist(outputDir,report);
try
    data=load(baselineMat,'baseline');
    assert(isfield(data,'baseline'),'ipm_cr:RefactorBaseline','The MAT file must contain baseline.');
    baseline=data.baseline;
    assert(isstruct(baseline)&&isscalar(baseline)&& ...
        all(isfield(baseline,{'cases','queryCloud','tauEnd','baselineRuntime','authoritySHA'}))&& ...
        numel(baseline.cases)==2,'ipm_cr:RefactorBaseline','Use the complete two-case frozen baseline.');
    validateattributes(baseline.tauEnd,{'numeric'},{'scalar','real','finite','positive'});
    report.baselineRuntime=baseline.baselineRuntime;report.tauEnd=baseline.tauEnd;
    cfg=ipm_c.config();authoritySHA=ipm_c.output.sha256(cfg.authorityFile);
    require_bitwise(authoritySHA,baseline.authoritySHA,'authoritySHA');
    report.authoritySHA=authoritySHA;report.status='running';
    for j=1:numel(baseline.cases)
        caseTimer=tic;reference=baseline.cases(j);
        label=sprintf('case(%d)',j);
        [state,g,wall,method]=ipm_cr.explicit.initialize(reference.method);
        row=struct('caseIndex',j,'geometryKeyIdentical',false,'wallOpsIdentical',false, ...
            'methodIdentical',false,'authorityIdentical',false,'initialSnapshotIdentical',false, ...
            'initialSampleIdentical',false,'initialDiagnosticIdentical',false, ...
            'finalSnapshotIdentical',false,'finalSampleIdentical',false, ...
            'finalDiagnosticIdentical',false,'constantDerivative',struct(), ...
            'stepsAdvanced',0,'maximumStageCfl',0,'retries',0,'seconds',0,'passed',false);
        require_bitwise(method,reference.method,[label,'.method']);row.methodIdentical=true;
        require_bitwise(g.geometryKey,reference.geometryKey,[label,'.geometryKey']);row.geometryKeyIdentical=true;
        require_bitwise(wall,reference.wallOps,[label,'.wallOps']);row.wallOpsIdentical=true;
        require_bitwise(authoritySHA,reference.authoritySHA,[label,'.authoritySHA']);row.authorityIdentical=true;
        require_bitwise(snapshot(state),reference.initial,[label,'.initial']);row.initialSnapshotIdentical=true;
        [sample,diagnostic]=ipm_cr.explicit.sampleState(state,g,method,baseline.queryCloud);
        require_bitwise(sample,reference.initialSample,[label,'.initialSample']);row.initialSampleIdentical=true;
        require_bitwise(diagnostic,reference.initialDiagnostic,[label,'.initialDiagnostic']);row.initialDiagnosticIdentical=true;
        row.constantDerivative=constant_check(g);
        while state.scale.canonicalTime<baseline.tauEnd
            [state,step]=ipm_cr.explicit.advance(state,g,wall,method,baseline.tauEnd);
            assert(strcmp(step.timeIntegrator,'rk4')&& ...
                step.realizedMaximumStageCfl<=method.cfl*(1+128*eps), ...
                'ipm_cr:RefactorCfl','The regression must retain complete-stage explicit RK4 CFL enforcement.');
            row.stepsAdvanced=row.stepsAdvanced+1;
            row.maximumStageCfl=max(row.maximumStageCfl,step.realizedMaximumStageCfl);
            row.retries=row.retries+step.retries;
        end
        assert(state.scale.canonicalTime==baseline.tauEnd,'ipm_cr:RefactorEndpoint');
        require_bitwise(snapshot(state),reference.final,[label,'.final']);row.finalSnapshotIdentical=true;
        [sample,diagnostic]=ipm_cr.explicit.sampleState(state,g,method,baseline.queryCloud);
        require_bitwise(sample,reference.finalSample,[label,'.finalSample']);row.finalSampleIdentical=true;
        require_bitwise(diagnostic,reference.finalDiagnostic,[label,'.finalDiagnostic']);row.finalDiagnosticIdentical=true;
        row.seconds=toc(caseTimer);row.passed=true;report.cases=append(report.cases,row);
        report.elapsedSeconds=toc(timer);persist(outputDir,report);
        fprintf('CR REFACTOR case %d: numerical state and query readout identical, %d steps, CFL %.6g, %.3fs\n', ...
            j,row.stepsAdvanced,row.maximumStageCfl,row.seconds);
    end
    report.sourceHashesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
    assert(report.sourceHashesUnchanged,'ipm_cr:RefactorSources','Recorded source files changed during the regression.');
    report.passed=all([report.cases.passed]);report.status='passed_exact_maintenance_regression';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function out=snapshot(state)
f=state.flow;
out=struct('rho',state.rho,'residual',state.residual,'scale',state.scale,'step',state.step, ...
    'flow',struct('q',f.q,'F',f.F,'Q',f.Q,'residualRhs',f.residualRhs, ...
    'physical',struct('core',f.physical.core,'polar',f.physical.polar), ...
    'rates',[f.c_l,f.c_omega,f.c_r,f.delta,f.Uref,f.transportOrigin]));
end

function info=constant_check(g)
% An exact mathematical invariant, without mirroring the implementation.
value=3.125;field=struct('core',value*ones(size(g.core.X)), ...
    'polar',value*ones(size(g.polar.X)));
d=ipm_cr.explicit.differentiate(field,g);
assert(all(d.core.x==0,'all')&&all(d.core.y==0,'all')&& ...
    all(d.polar.x==0,'all')&&all(d.polar.y==0,'all'), ...
    'ipm_cr:RefactorConstant','The derivative of a constant must be exactly zero.');
info=struct('inputConstant',value,'coreXExactlyZero',true,'coreYExactlyZero',true, ...
    'polarXExactlyZero',true,'polarYExactlyZero',true,'passed',true);
end

function require_bitwise(actual,expected,label)
path=first_mismatch(actual,expected,label);
assert(isempty(path),'ipm_cr:RefactorBitwise','Frozen numerical regression differs at %s.',path);
end

function path=first_mismatch(a,b,label)
path='';
if ~strcmp(class(a),class(b))||~isequal(size(a),size(b))||issparse(a)~=issparse(b)
    path=[label,' (class, shape or storage)'];return;
end
if isstruct(a)
    names=fieldnames(a);
    if ~isequal(names,fieldnames(b)),path=[label,' (fields)'];return;end
    for j=1:numel(a)
        for k=1:numel(names)
            name=names{k};path=first_mismatch(a(j).(name),b(j).(name),sprintf('%s(%d).%s',label,j,name));
            if ~isempty(path),return;end
        end
    end
elseif iscell(a)
    for j=1:numel(a)
        path=first_mismatch(a{j},b{j},sprintf('%s{%d}',label,j));
        if ~isempty(path),return;end
    end
elseif isnumeric(a)
    if issparse(a)
        [ar,ac,av]=find(a);[br,bc,bv]=find(b);
        if ~isequal(ar,br)||~isequal(ac,bc),path=[label,' (sparse indices)'];return;end
        path=first_mismatch(av,bv,[label,' (sparse values)']);
    elseif isreal(a)~=isreal(b)|| ...
            ~isequal(typecast(real(a(:)),'uint8'),typecast(real(b(:)),'uint8'))|| ...
            ~isequal(typecast(imag(a(:)),'uint8'),typecast(imag(b(:)),'uint8'))
        path=[label,' (numeric bytes)'];
    end
elseif ~isequaln(a,b)
    path=[label,' (value)'];
end
end

function inventory=source_inventory(root)
paths={mfilename('fullpath')};if ~endsWith(paths{1},'.m'),paths{1}=[paths{1},'.m'];end
for folder={'test/+ipm_cr','time_evolution/+ipm','time_evolution/+ipm_c'}
    files=dir(fullfile(root,folder{1},'**','*'));
    for j=1:numel(files)
        [~,~,extension]=fileparts(files(j).name);
        if ~files(j).isdir&&ismember(extension,{'.m','.c','.cpp',['.',mexext]})
            paths{end+1}=fullfile(files(j).folder,files(j).name); %#ok<AGROW>
        end
    end
end
cfg=ipm_c.config();paths{end+1}=cfg.authorityFile;
paths=sort(unique(paths));inventory=struct([]);
for j=1:numel(paths)
    item=dir(paths{j});inventory=append(inventory,struct('path',paths{j},'bytes',item.bytes, ...
        'sha256',ipm_c.output.sha256(paths{j})));
end
end

function persist(folder,report)
temporary=fullfile(folder,'refactor.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'refactor.mat'),'f');assert(ok,'ipm_cr:RefactorSave','%s',message);
temporary=fullfile(folder,'refactor.json.writing');fid=fopen(temporary,'w');assert(fid>=0,'ipm_cr:RefactorSave');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'refactor.json'),'f');assert(ok,'ipm_cr:RefactorSave','%s',message);
end

function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
