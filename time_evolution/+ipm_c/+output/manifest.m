function rows=manifest()
%MANIFEST Bind the complete shared ipm dependency as well as the new solver.
parent=fileparts(fileparts(fileparts(mfilename('fullpath'))));
files=[dir(fullfile(parent,'+ipm_c','**','*'));dir(fullfile(parent,'+ipm','**','*'))];
files=files(~[files.isdir]);names=string({files.name});
files=files(endsWith(names,{'.m','.cpp','.c'})|contains(names,'.mex'));
for name={'ipm.mesh.fdMatrix','ipm.field.weno5FluxDerivative', ...
        'ipm.evolve.globalSplineMaximum','ipm.evolve.wallGradientIntegral','ipm.evolve.wallGradientCentroid'}
    actual=which(name{1});
    assert(startsWith(actual,[parent,filesep,'+ipm',filesep]),'ipm_c:ShadowedDependency','Validation failed.');
end
[~,order]=sort(string({files.folder})+"/"+string({files.name}));files=files(order);
rows=repmat(struct('path','','sha256',''),numel(files),1);
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    rows(i)=struct('path',path,'sha256',ipm_c.output.sha256(path));
end
end
