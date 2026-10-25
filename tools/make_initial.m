function [initial,report]=make_initial(name,overrides)
%MAKE_INITIAL Create a new Gaussian t0 authority, never overwrite a case.
% [INITIAL,REPORT]=MAKE_INITIAL(NAME,PARAMETERS) returns a project-relative
% manifest accepted by method.initial. PARAMETERS may contain amplitude,
% centerX, widthX, widthY, skewWeight (or their gaussian... physics names).
% The original authority supplies fixed measurement axes, not its I/P/C.
% This creates a new IVP identity; previous trajectory qualifications do not
% transfer. It performs no Poisson solve or time advancement.
if nargin<2,overrides=struct();end
assert(nargin>=1&&((ischar(name)&&isrow(name))||(isstring(name)&&isscalar(name))), ...
    'ipm_cr:InitialName','Supply a new case name.');
name=char(name);
assert(~isempty(regexp(name,'^[A-Za-z][A-Za-z0-9_-]{0,63}$','once')), ...
    'ipm_cr:InitialName','Use a letter followed by at most 63 letters, digits, underscores or hyphens.');
assert(isempty(regexpi(name,'^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])$','once')), ...
    'ipm_cr:InitialName','Choose a portable name that is not a reserved Windows device name.');
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:InitialParameters', ...
    'Use a scalar Gaussian parameter structure, without eta arrays or callbacks.');
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
base=ipm_cr.problem();p=base.physics;
short={'amplitude','centerX','widthX','widthY','skewWeight'};
full={'gaussianAmplitude','gaussianCenterX','gaussianWidthX','gaussianWidthY','gaussianSkewWeight'};
assert(all(ismember(fieldnames(overrides),[short,full])), ...
    'ipm_cr:InitialParameters','Only the five Gaussian parameters may be changed.');
for k=1:numel(short)
    assert(~(isfield(overrides,short{k})&&isfield(overrides,full{k})), ...
        'ipm_cr:InitialParameters','Duplicate Gaussian parameter %s.',short{k});
    if isfield(overrides,short{k}),value=overrides.(short{k});
    elseif isfield(overrides,full{k}),value=overrides.(full{k});
    else,value=p.(full{k});end
    validateattributes(value,{'numeric'},{'scalar','real','finite'});
    p.(full{k})=double(value);
end
assert(p.gaussianAmplitude>0&&p.gaussianWidthX>0&&p.gaussianWidthY>0&& ...
    p.gaussianSkewWeight>=0&&p.gaussianSkewWeight<=1, ...
    'ipm_cr:InitialParameters','Positive amplitude/widths and skewWeight in [0,1] are required.');
cases=fullfile(root,'data','cases');target=fullfile(cases,name);
canonicalRoot=char(java.io.File(root).getCanonicalPath());
canonicalTarget=char(java.io.File(target).getCanonicalPath());
assert(startsWith(canonicalTarget,[canonicalRoot,filesep]), ...
    'ipm_cr:InitialPath','The case directory resolves outside the project.');
assert(~isfolder(target)&&~isfile(target),'ipm_cr:InitialExists','Case already exists: %s',target);
if ~isfolder(cases),[ok,message]=mkdir(cases);assert(ok,'ipm_cr:InitialWrite','%s',message);end
stage=tempname(cases);[ok,message]=mkdir(stage);assert(ok,'ipm_cr:InitialWrite','%s',message);
cleanup=onCleanup(@()remove_stage(stage)); %#ok<NASGU>
loaded=load(base.authorityFile,'frame');old=loaded.frame;
[X,Y]=meshgrid(old.x,old.y);
rho=ipm_cr.model.initial(struct('X',X,'Y',Y,'ny',size(X,1),'nx',size(X,2)),p);
Dx=ipm_cr.mesh.fd(old.x,1,7);
frame=struct('x',old.x,'y',old.y,'rho',rho,'w',rho*Dx', ...
    'step',0,'tau',0,'t',0,'Cx',1,'Cy',1,'Comega',1,'Xshift',0); %#ok<NASGU>
matFile=fullfile(stage,'t0.mat');save(matFile,'frame','-v7.3');
manifest=jsondecode(fileread(fullfile(root,base.manifestFile)));
manifest.physics=p;
if isfield(manifest.grid,'originFitDegree'),manifest.grid=rmfield(manifest.grid,'originFitDegree');end
for field={'referenceWallIntegral','referenceWallGradientCentroid','referenceWallPeak'}
    if isfield(manifest.scaling,field{1}),manifest.scaling=rmfield(manifest.scaling,field{1});end
end
candidate=base;candidate.physics=p;candidate.scaling=manifest.scaling;
candidate.authorityFile=matFile;candidate.authoritySha=ipm_cr.io.sha(matFile);
[wall,scale,refs,guess]=ipm_cr.model.reference(candidate);
assert(all(isfinite([refs.I0,refs.P0,refs.C0,guess]))&&refs.I0>0&&refs.P0>0, ...
    'ipm_cr:InitialGauge','The fixed measurement axes must resolve positive, finite I/P/C gauges.');
manifest.scaling.referenceWallIntegral=refs.I0;
manifest.scaling.referenceWallGradientCentroid=refs.C0;
candidate.scaling=manifest.scaling;
[verifiedWall,verifiedScale,verifiedRefs,verifiedGuess]=ipm_cr.model.reference(candidate);
assert(isequaln(wall,verifiedWall)&&isequaln(scale,verifiedScale)&& ...
    isequaln(refs,verifiedRefs)&&isequaln(guess,verifiedGuess), ...
    'ipm_cr:InitialReference','Persisted reference verification changed the Gaussian authority.');
initial=['data/cases/',name,'/initial.json'];
manifest.file=['data/cases/',name,'/t0.mat'];manifest.sha256=candidate.authoritySha;
ipm_cr.io.json(fullfile(stage,'initial.json'),manifest);
assert(strcmp(ipm_cr.io.sha(base.authorityFile),base.authoritySha)&& ...
    strcmp(ipm_cr.io.sha(fullfile(root,base.manifestFile)),base.manifestSha), ...
    'ipm_cr:InitialBaseChanged','The fixed measurement authority changed during generation.');
assert(~isfolder(target)&&~isfile(target),'ipm_cr:InitialExists','Case appeared during generation.');
% No REPLACE_EXISTING option: a concurrent case creator cannot be overwritten.
java.nio.file.Files.move(java.io.File(stage).toPath(),java.io.File(target).toPath(), ...
    javaArray('java.nio.file.CopyOption',0));
cfg=ipm_cr.problem(initial);[~,~,finalRefs]=ipm_cr.model.reference(cfg);
assert(isequaln(finalRefs,refs),'ipm_cr:InitialReference','Published reference differs from staged verification.');
report=struct('kind','generated_Gaussian_t0_authority_v1','initial',initial, ...
    'authority',ipm_cr.io.identity('authority',root,initial),'physics',p,'references',refs, ...
    'measurementAxesFrom',base.manifestFile,'measurementAuthoritySHA256',base.authoritySha, ...
    'newProblemIdentity',true,'timeAdvanced',0,'poissonCalls',0, ...
    'previousTrajectoryQualificationInherited',false);
fprintf('Created %s; new I0=%.17g P0=%.17g C0=%.17g\n',initial,refs.I0,refs.P0,refs.C0);
end

function remove_stage(folder)
if isfolder(folder),rmdir(folder,'s');end
end
