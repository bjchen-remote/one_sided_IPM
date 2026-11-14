function report=test_config(output)
%TEST_CONFIG Portable case defaults, schema and method identity.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,output=tempname;end
c=ipm_cr.config('case');[r,m]=ipm_cr.config(c,'run');
assert(c.domain.rmax==1e10&&m.geometry.grid.rmax==1e10&&m.cfl==.4);
assert(strcmp(m.geometry.preset,'medium')&&m.geometry.staircase&& ...
    m.geometry.interpolationWidth==7&&strcmp(m.initial,'data/initial_small_source.json'));
assert(m.remesh.enabled&&isequal(m.remesh.level,[.25,.8])&&m.remesh.threshold==.25);
assert(isequaln(m,ipm_cr.config(m)));
assert(isequaln(r,ipm_cr.config(jsondecode(jsonencode(c)),'run')));
example=fullfile(root,'examples','bspline4.json');ipm_cr.config(jsondecode(fileread(example)),'run');
coarse=ipm_cr.config(struct('geometry',struct('preset','coarse')));
assert(coarse.remesh.stopTau==50);
coarse=ipm_cr.config(struct('geometry',struct('preset','coarse'),'remesh',struct('stopTau',[])));
assert(isempty(coarse.remesh.stopTau));
c.domain.rmax=1e8;c.run.tauEnd=12;[r,m]=ipm_cr.config(c,'run');
assert(m.geometry.grid.rmax==1e8&&r.tauEnd==12);
bad=c;bad.method.geometry.grid.rmax=1e10;
expect(@()ipm_cr.config(bad,'run'),'ipm_cr:DomainRadius');
expect(@()ipm_cr.config(struct('poissonBackend','fem')),'ipm_cr:Config');
expect(@()ipm_cr.config(struct('timeIntegrator','rk4')),'ipm_cr:Config');
expect(@()ipm_cr.config(struct('remesh',struct('intervalTau',2))),'ipm_cr:RemeshOptions');
expect(@()ipm_cr.config(struct(),'run'),'ipm_cr:CaseFormat');
bad=c;bad.run.outputDir=123;
expect(@()ipm_cr.config(bad,'run'),'ipm_cr:Config');
if ~isfolder(output),mkdir(output);end
report=struct('passed',true,'format',c.format,'cfl',m.cfl,'ghostDonorWidth',m.geometry.interpolationWidth);
ipm_cr.io.identity('json',fullfile(output,'config.json'),report);
end

function expect(callback,id)
try
    callback();error('ipm_cr:ExpectedConfigError','Expected %s.',id);
catch problem
    assert(strcmp(problem.identifier,id),'ipm_cr:UnexpectedConfigError','Expected %s, got %s.',id,problem.identifier);
end
end
