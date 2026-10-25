function initial=set_initial(c,rho0)
%SET_INITIAL Materialize a density expression as an immutable input.
% rho0=@(x,y)... must return matching finite real arrays (a scalar is allowed).
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
[~,method]=ipm_cr.config(c,'run');g=ipm_cr.mesh.build(method.geometry);
base=ipm_cr.problem();old=load(base.authorityFile,'frame');
assert(isa(rho0,'function_handle'),'ipm_cr:InitialFunction','Write rho0 as @(x,y) your_expression.');
evaluate=rho0;description=func2str(rho0);
[X,Y]=meshgrid(old.frame.x,old.frame.y);rho=values(X,Y);
native=struct('geometryKey',g.geometryKey,'core',values(g.core.X,g.core.Y), ...
    'polar',values(g.polar.X,g.polar.Y),'origin',values(0,0));
frame=struct('x',old.frame.x,'y',old.frame.y,'rho',rho, ...
    'step',0,'tau',0,'t',0,'Cx',1,'Cy',1,'Comega',1,'Xshift',0,'native',native); %#ok<NASGU>
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native([base.manifestSha,jsonencode(g.geometryKey)],'UTF-8'));
for v={rho,native.core,native.polar,native.origin}
    digest.update(typecast(double(v{1}(:)),'uint8'));
end
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
name=['rho_',hash(1:24)];initial=['data/cases/',name,'/initial.json'];
parent=fullfile(root,'data','cases');target=fullfile(parent,name);
if isfolder(target)
    cfg=ipm_cr.problem(initial);saved=load(cfg.authorityFile,'frame');
    assert(isequaln(saved.frame,frame),'ipm_cr:InitialIdentity','Existing initial data differ; preserve and inspect %s.',initial);
    return
end
if ~isfolder(parent),mkdir(parent);end
stage=tempname(parent);mkdir(stage);cleanup=onCleanup(@()remove_stage(stage)); %#ok<NASGU>
file=fullfile(stage,'t0.mat');save(file,'frame','-v7.3');
manifest=jsondecode(fileread(fullfile(root,'data','initial.json')));
manifest.physics=struct('initialCondition','sampled_density','symmetryMode','half_plane');
manifest.file=['data/cases/',name,'/t0.mat'];manifest.sha256=ipm_cr.io.sha(file);
manifest.initialDescription=description;
manifest.representation='Exact materialized native samples; transport mesh identity is pinned. No function handle required on restart.';
for key={'referenceWallIntegral','referenceWallGradientCentroid','referenceWallPeak'}
    if isfield(manifest.scaling,key{1}),manifest.scaling=rmfield(manifest.scaling,key{1});end
end
cfg=base;cfg.physics=manifest.physics;cfg.scaling=manifest.scaling;
cfg.authorityFile=file;cfg.authoritySha=manifest.sha256;
[~,~,refs]=ipm_cr.model.reference(cfg);
assert(all(isfinite([refs.I0,refs.P0,refs.C0]))&&refs.I0>0&&refs.P0>0, ...
    'ipm_cr:InitialGauge','This initial density does not admit the configured I/P/C normalization.');
manifest.scaling.referenceWallIntegral=refs.I0;
manifest.scaling.referenceWallGradientCentroid=refs.C0;
ipm_cr.io.json(fullfile(stage,'initial.json'),manifest);
java.nio.file.Files.move(java.io.File(stage).toPath(),java.io.File(target).toPath(), ...
    javaArray('java.nio.file.CopyOption',0));
fprintf('Initial density saved: %s\n',initial);

    function R=values(X,Y)
        R=evaluate(X,Y);
        if isnumeric(R)&&isscalar(R),R=repmat(R,size(X));end
        assert(isnumeric(R)&&isreal(R)&&isequal(size(R),size(X))&&all(isfinite(R),'all'), ...
            'ipm_cr:InitialFunction','rho0(x,y) must return matching finite real values.');
        R=double(R);
    end

end

function remove_stage(folder)
if isfolder(folder),rmdir(folder,'s');end
end
