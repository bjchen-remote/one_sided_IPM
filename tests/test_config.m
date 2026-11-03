function report=test_config(outputDir)
%TEST_CONFIG Public case roundtrips and mathematical input contracts.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'result','test_config');end
if ~isfolder(outputDir),mkdir(outputDir);end
c=ipm_cr.config('case');[r,m,requested]=ipm_cr.config(c,'run');
assert(strcmp(r.outputDir,fullfile(pwd,'result'))&&r.tauEnd==200&&r.wallSeconds==600);
assert(c.domain.rmax==1e10&&m.geometry.grid.rmax==1e10&&m.solverOptions.rOuter==1e10);
assert(strcmp(m.timeIntegrator,'rk4')&&m.cfl==.3&& ...
    isequal(fieldnames(c.stop),{'physicalTime'})&& ...
    isequal(fieldnames(r.stop),{'physicalTime'}));
[~,minimal]=ipm_cr.config(struct('format','ipm_cr_case_v1'),'run');
[~,flat]=ipm_cr.config(struct(),'run');defaults=ipm_cr.config();
exampleFile=fullfile(root,'examples','case.json');
[~,example]=ipm_cr.config(jsondecode(fileread(exampleFile)),'run');
for value={m,minimal,flat,defaults,example}
    current=value{1};
    assert(strcmp(current.geometry.preset,'medium')&&strcmp(current.poissonBackend,'fem')&& ...
        strcmp(current.transport.reconstruction,'linear5')&&current.remesh.enabled, ...
        'ipm_cr:ConfigMainlineDefault','All default entrypoints must select the same production method.');
    assert(current.remesh.level==.25&&current.remesh.threshold==.2&& ...
        current.remesh.checkEveryTau==.5&&current.remesh.minCheckEveryTau==.005&& ...
        current.remesh.minCheckEveryTau<=current.remesh.checkEveryTau&& ...
        numel(fieldnames(current.remesh))==9&&~isfield(current.remesh,'intervalTau'), ...
        'ipm_cr:ConfigRemeshTrigger', ...
        'Remeshing must use source-region changes with the original adaptive check schedule, not periodic interpolation.');
    assert(~isfield(current,'minDt')&& ...
        ~isfield(current.transport,'outerInflowPolicy')&& ...
        ~isfield(current.transport,'outerInflowTolerance'), ...
        'ipm_cr:ConfigNumericalGateFields','Production methods must not contain numerical rejection switches.');
end
partial=ipm_cr.config(struct('remesh',struct('enabled',true), ...
    'transport',struct('centralDissipation',.25)));
assert(strcmp(partial.transport.reconstruction,'linear5')&& ...
    partial.transport.centralDissipation==.25&& ...
    isequaln(partial,ipm_cr.config(partial)),'ipm_cr:ConfigPartialRoundtrip');
[r2,m2,p2]=ipm_cr.config(jsondecode(jsonencode(c)),'run');
assert(isequaln(r,r2)&&isequaln(m,m2)&&isequaln(requested,p2));
spaceFolder=char(java.io.File(fullfile(outputDir,'path with spaces')).getAbsolutePath());
if ~isfolder(spaceFolder),mkdir(spaceFolder);end
paths={spaceFolder,fullfile('result','relative path with spaces')};
expected={spaceFolder,fullfile(pwd,paths{2})};
if ispc
    paths=[paths,{'C:\ipm cr case\output','\\server\share\ipm cr case'}];
    expected=[expected,paths(3:4)];
end
for k=1:numel(paths)
    pathCase=c;pathCase.run.outputDir=paths{k};resolved=ipm_cr.config(pathCase,'run');
    assert(strcmp(resolved.outputDir,expected{k}),'ipm_cr:ConfigOutputPath', ...
        'Absolute paths must be retained, and relative paths must resolve from the current directory.');
end
% A custom radius uses the explicit generator, not the fixed-axis preset.
c.method.geometry=struct();c.domain.rmax=1e8;c.run.outputDir='result/config_case';c.run.tauEnd=12;
c.stop.physicalTime=.5;c.save.label='test';
[r,m]=ipm_cr.config(c,'run');
assert(m.geometry.grid.rmax==1e8&&m.solverOptions.rOuter==1e8&&r.tauEnd==12);
assert(strcmp(r.outputDir,fullfile(pwd,'result/config_case'))&&r.stop.physicalTime==.5);
assert(isequal(fieldnames(r.stop),{'physicalTime'})&&strcmp(r.save.label,'test'));
[m,cfg]=ipm_cr.config();assert(isequaln(m.background,cfg.background));
a=ipm_cr.config(struct('mesh',struct(),'poisson',struct('rOuter',1e10)));
b=ipm_cr.config(struct('geometry',struct(),'solverOptions',struct('rOuter',1e10)));
assert(isequaln(a,b));
expect(@()ipm_cr.config(struct('mesh',struct(),'geometry',struct())),'ipm_cr:Config');
expect(@()ipm_cr.config(struct('typo',1)),'ipm_cr:Config');
expect(@()ipm_cr.config(struct('remesh',struct('intervalTau',2))),'ipm_cr:RemeshOptions');
expect(@()ipm_cr.config(struct('remesh',struct('checkEveryTau',.1, ...
    'minCheckEveryTau',.2))),'ipm_cr:RemeshOptions');
bad=c;bad.method.solverOptions.rOuter=1e10;
expect(@()ipm_cr.config(bad,'run'),'ipm_cr:DomainRadius');
bad=c;bad.run.typo=1;expect(@()ipm_cr.config(bad,'run'),'ipm_cr:Config');
bad=c;bad.stop.typo=1;expect(@()ipm_cr.config(bad,'run'),'ipm_cr:Config');
bad=c;bad.save.typo=true;expect(@()ipm_cr.config(bad,'run'),'ipm_cr:Config');
bad=c;bad.run.outputDir=123;expect(@()ipm_cr.config(bad,'run'),'ipm_cr:Config');
bad=c;bad.save.label=123;expect(@()ipm_cr.config(bad,'run'),'ipm_cr:Config');
expect(@()ipm_cr.config(struct('geometry',struct('grid',struct('originFitDegree',6)))), ...
    'ipm_cr:RetiredOption');
expect(@()ipm_cr.config(struct('explicitOptions',struct('cfl',.3),'cfl',.2),'run'), ...
    'ipm_cr:Config');
report=struct('kind','ipm_cr_case_configuration_v1','passed',true,'pdeSolved',false, ...
    'caseRoundtrip',true,'mainlineDefaultsAcrossEntrypoints',true, ...
    'partialMethodRoundtrip',true,'radiusCoupled',true,'aliasesAtOneBoundary',true, ...
    'shapeChangeRemeshDefaults',true,'originalAdaptiveCheckSchedule',true, ...
    'periodicInterpolationOptionRemoved',true, ...
    'physicalTimeOnlyStop',true,'productionNumericalGateFieldsAbsent',true, ...
    'nativeAbsolutePathAndSpaces',true,'relativeOutputUsesCurrentDirectory',true, ...
    'windowsDriveAndUNCChecksRun',ispc, ...
    'duplicateAndUnknownRejected',true,'backgroundFromSelectedDatum',true);
ipm_cr.io.json(fullfile(outputDir,'config.json'),report);
end

function expect(fn,id)
try,fn();error('ipm_cr:ExpectedFailure','Expected %s.',id);
catch problem,assert(strcmp(problem.identifier,id),'ipm_cr:WrongFailure', ...
        'Expected %s, got %s.',id,problem.identifier);end
end
