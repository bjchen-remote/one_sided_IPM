function data=data(source,selector,mode)
%DATA Load a polar trajectory and select one native field slice.
% d=ipm_c.data('/absolute/run/data.mat');       % latest spatial slice
% d=ipm_c.data(d,3.5);                         % nearest saved tau
% d=ipm_c.data(d,2,'index');                    % exact frame index
% d=ipm_c.data(d,3.5,'exact');                  % require a saved tau
% d.tau, d.t and scalar diagnostics remain whole-trajectory vectors.
% d.selected has the chosen clocks, scale and rates. Selected field arrays
% are on d.grid; physical aliases are computed from that frame's own scale.
if ischar(source)||(isstring(source)&&isscalar(source))
    file=char(source);if isfolder(file),file=fullfile(file,'data.mat');end
    assert(isfile(file),'ipm_c:DatasetMissing','Dataset does not exist: %s',file);
    saved=load(file,'data');assert(isfield(saved,'data'),'ipm_c:DatasetSchema');data=saved.data;
    [ok,attributes]=fileattrib(file);assert(ok,'ipm_c:DatasetPath');
    if ~isfield(data,'meta'),data.meta=struct();end
    data.meta.datasetFile=attributes.Name;data.meta.baseFolder=fileparts(attributes.Name);
    data=upgrade_legacy(data);
else
    assert(isstruct(source)&&isscalar(source),'ipm_c:DatasetInput');data=source;
end
assert(isfield(data,'schemaVersion')&&ismember(data.schemaVersion,[3,4])&& ...
    isfield(data,'frames')&&~isempty(data.frames),'ipm_c:DatasetSchema', ...
    'Expected a nonempty polar trajectory loaded by ipm_c.data.');
if nargin<2||isempty(selector)
    selector=numel(data.frames);mode='index';
elseif nargin<3||isempty(mode)
    mode='nearest';
end
mode=validatestring(char(mode),{'nearest','exact','index'});
validateattributes(selector,{'numeric'},{'scalar','real','finite'});
times=reshape([data.frames.tau],[],1);
switch mode
    case 'index'
        validateattributes(selector,{'numeric'},{'integer','>=',1,'<=',numel(times)});index=selector;
    case 'exact'
        index=find(times==selector,1,'first');
        assert(~isempty(index),'ipm_c:DatasetTimeNotSaved', ...
            'No spatial slice at tau=%.17g. Inspect data.fieldTau or use nearest.',selector);
    case 'nearest'
        [~,index]=min(abs(times-selector));
end
entry=data.frames(index);file=resolve_path(data.meta.baseFolder,entry.file);
assert(isfile(file),'ipm_c:DatasetFrameMissing','Saved frame is missing: %s',file);
loaded=load(file,'frame');f=loaded.frame;
if isfield(f,'schemaVersion')&&f.schemaVersion>=3
 frameKinds={'ipm_c_total_rho_polar_frame_v3','ipm_c_centered_density_polar_frame_v4'};
 mixed=data.schemaVersion==4&&isfield(data.meta,'mixedFrameSchemas')&&any(data.meta.mixedFrameSchemas(:));
 assert(ismember(f.schemaVersion,[3,4])&&(f.schemaVersion==data.schemaVersion|| ...
     (mixed&&f.schemaVersion==3))&&strcmp(f.kind,frameKinds{f.schemaVersion-2}), ...
     'ipm_c:DatasetFrameSchema','Mixed frame schemas require explicit migration metadata.');
end
if isfield(f,'schemaVersion')&&ismember(f.schemaVersion,[3,4])
    meshFile=resolve_path(data.meta.baseFolder,f.meshFile);saved=load(meshFile,'grid');g=saved.grid;
    assert(strcmp(g.meshId,f.meshId)&&strcmp(mesh_key(g.r,g.theta),f.meshId), ...
        'ipm_c:DatasetMeshMismatch','The saved coordinates do not match this frame mesh identity.');
    f.r=g.r;f.theta=g.theta;[f.X,f.Y]=cartesian_nodes(g.r,g.theta);
else
    g=struct('r',f.r(:),'theta',f.theta(:).','nr',numel(f.r), ...
        'ntheta',numel(f.theta),'rmax',f.r(end),'domain','upper_half_disk', ...
        'originIsSinglePhysicalPoint',true,'meshId','legacy_per_frame_mesh');
end
assert(f.step==entry.step&&f.tau==entry.tau,'ipm_c:DatasetFrameMismatch');
data.grid=g;data.r=f.r(:);data.theta=f.theta(:).';data.X=f.X;data.Y=f.Y;data.scale=f.scale;
names={'rho','rho_X','rho_Y','rho_tau','Q','psi','u1','u2','vx','vy'};
for k=1:numel(names)
    name=names{k};if isfield(f,name),data.(name)=f.(name);else,data.(name)=NaN(size(f.rho));end
end
for name={'rho_centered','rho_origin','rho_tau_centered','rho_origin_tau','transportOrigin'}
    if isfield(f,name{1}),data.(name{1})=f.(name{1});
    elseif isfield(data,name{1}),data=rmfield(data,name{1});end
end
data.rho_x=data.rho_X;data.rho_y=data.rho_Y;data.F=data.rho_tau;
data.rho_xt=data.Q;data.transport_u1=data.vx;data.transport_u2=data.vy;
data.vmag=hypot(data.u1,data.u2);data.transport_speed=hypot(data.vx,data.vy);
s=f.scale;Cl=exp(s.logC_l);Co=exp(s.logC_omega);
data.physical_X=(data.X-s.X_shift)/Cl;data.physical_Y=data.Y/Cl;
if isfield(f,'rho_centered')
    data.rho_physical=data.rho_centered/Co+data.rho_origin/Co;
else,data.rho_physical=data.rho/Co;end
data.rho_x_physical=(Cl/Co)*data.rho_X;
data.rho_y_physical=(Cl/Co)*data.rho_Y;
data.u1_physical=data.u1/Co;data.u2_physical=data.u2/Co;data.psi_physical=data.psi/(Cl*Co);
data.vmag_physical=data.vmag/Co;
frameSchema=data.schemaVersion;if isfield(f,'schemaVersion'),frameSchema=f.schemaVersion;end
data.selected=struct('schemaVersion',frameSchema,'frameIndex',index,'step',f.step,'tau',f.tau,'t',f.t, ...
    'scale',s,'C_l',Cl,'C_omega',Co,'c_l',f.c_l,'c_omega',f.c_omega,'c_r',f.c_r, ...
    'physicalTimeLow',scale_field(s,'physicalTimeLow',0), ...
    'densityOrigin',scale_field(s,'densityOrigin',0), ...
    'meshId',g.meshId,'mode',mode,'requested',selector,'timeInterpolationUsed',false, ...
    'file',entry.file,'source','native_accepted_state');
if isfield(f,'diagnostics'),data.selected.diagnostics=f.diagnostics;end
if isfield(f,'remeshCount'),data.selected.remeshCount=f.remeshCount;end
data.selectedIndex=index;data.selectedTau=f.tau;data.selectedStep=f.step;
data.meta.selection=data.selected;
end

function [X,Y]=cartesian_nodes(r,theta)
c=cos(theta(:).');s=sin(theta(:).');c([1,end])=[1,-1];s([1,end])=0;
if mod(numel(theta),2)==1
    middle=(numel(theta)+1)/2;
    if theta(middle)==pi/2,c(middle)=0;s(middle)=1;end
end
X=r(:)*c;Y=r(:)*s;
end

function path=resolve_path(base,file)
if startsWith(file,filesep),path=file;else,path=fullfile(base,file);end
end

function id=mesh_key(r,theta)
d=java.security.MessageDigest.getInstance('SHA-256');
d.update(typecast(uint64([numel(r),numel(theta)]),'uint8'));
d.update(typecast(r(:),'uint8'));d.update(typecast(theta(:),'uint8'));
id=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2)',1,[]));
end

function data=upgrade_legacy(data)
if isfield(data,'schemaVersion')&&ismember(data.schemaVersion,[3,4]),return;end
assert(isfield(data,'kind')&&strcmp(data.kind,'ipm_c_polar_native_data_v2'), ...
    'ipm_c:DatasetSchema','Expected an ipm_c polar trajectory.');
base=data.meta.baseFolder;listing=dir(fullfile(base,'data_frames','frame_*.mat'));
entries=struct('file',{},'step',{},'tau',{},'t',{},'meshId',{},'meshFile',{});
for k=1:numel(listing)
    relative=fullfile('data_frames',listing(k).name);loaded=load(fullfile(base,relative),'frame');f=loaded.frame;
    if f.step>data.history(end,1),continue;end
    entries(end+1)=struct('file',relative,'step',f.step,'tau',f.tau,'t',f.t, ...
        'meshId','legacy_per_frame_mesh','meshFile',''); %#ok<AGROW>
end
assert(~isempty(entries),'ipm_c:DatasetEmpty');[~,order]=sort([entries.step]);entries=entries(order);
data.frames=entries;data.schemaVersion=3;data.meta.loadedLegacySchema=2;
names=data.historyColumns;
for k=1:numel(names),data.(names{k})=data.history(:,k);end
data.dt=data.H;data.c_l=data.cl;data.c_omega=data.cw;data.c_r=data.cr;
data.c=[data.cl,data.cw,data.cr];data.residual_inf=data.Qinf;
data.fieldTau=reshape([entries.tau],[],1);data.frameStep=reshape([entries.step],[],1);
end

function value=scale_field(scale,name,fallback)
if isfield(scale,name),value=scale.(name);else,value=fallback;end
end
