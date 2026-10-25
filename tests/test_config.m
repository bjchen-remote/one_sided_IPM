function report=test_config(outputDir)
%TEST_CONFIG Public case resolution and rejection checks, without a PDE solve.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'result','test_config');end
if ~isfolder(outputDir),mkdir(outputDir);end
c=ipm_cr.config('case');[r,m,requested]=ipm_cr.config(c,'run');
assert(strcmp(r.outputDir,fullfile(pwd,'result'))&&r.tauEnd==200&&r.wallSeconds==600);
assert(c.domain.rmax==1e10&&m.geometry.grid.rmax==1e10&&m.solverOptions.rOuter==1e10);
assert(strcmp(m.timeIntegrator,'rk4')&&m.cfl==.3&&~r.stop.steady.enabled);
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
c.domain.rmax=1e8;c.run.outputDir='result/config_case';c.run.tauEnd=12;
c.stop.physicalTime=.5;c.stop.steady.enabled=true;c.save.label='test';
[r,m]=ipm_cr.config(c,'run');
assert(m.geometry.grid.rmax==1e8&&m.solverOptions.rOuter==1e8&&r.tauEnd==12);
assert(strcmp(r.outputDir,fullfile(pwd,'result/config_case'))&&r.stop.physicalTime==.5);
assert(r.stop.steady.enabled&&strcmp(r.save.label,'test'));
[m,cfg]=ipm_cr.config();assert(isequaln(m.background,cfg.background));
a=ipm_cr.config(struct('mesh',struct(),'poisson',struct('rOuter',1e10)));
b=ipm_cr.config(struct('geometry',struct(),'solverOptions',struct('rOuter',1e10)));
assert(isequaln(a,b));
expect(@()ipm_cr.config(struct('mesh',struct(),'geometry',struct())),'ipm_cr:Config');
expect(@()ipm_cr.config(struct('typo',1)),'ipm_cr:Config');
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
    'caseRoundtrip',true,'radiusCoupled',true,'aliasesAtOneBoundary',true, ...
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
