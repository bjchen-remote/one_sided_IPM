function report=test_native_sampling(output)
%TEST_NATIVE_SAMPLING Compare native readout with its sparse MATLAB formula.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,output=tempname;end
g=ipm_cr.mesh.build(struct('preset','coarse','staircase',true,'interpolationWidth',7));
x=[0,.13,-.21,12,-19];y=[0,.37,.05,4,7];
f=struct('core',sin(g.core.X).*cos(g.core.Y), ...
    'polar',sin(g.polar.X).*cos(g.polar.Y));
[actual,map]=ipm_cr.mesh.interpolate.sample(f,g,x,y);
assert(map.allValid&&isequaln(actual,ipm_cr.mesh.interpolate.sample(f,map)));
map.referenceOnly=true;
reference=ipm_cr.mesh.interpolate.sample(f,map);
assert(max(abs(actual-reference),[],'all')<5e-14,'ipm_cr:NativeSampling');
constant=struct('core',ones(size(f.core)),'polar',ones(size(f.polar)));
assert(all(ipm_cr.mesh.interpolate.sample(constant,map)==1));
report=struct('passed',true,'samples',actual,'constantPreserved',true);
if ~isfolder(output),mkdir(output);end
ipm_cr.io.identity('json',fullfile(output,'native_sampling.json'),report);
end
