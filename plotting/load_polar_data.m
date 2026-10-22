function data=load_polar_data(source,tau,box,resolution)
%LOAD_POLAR_DATA Ordinary plotting struct: one Cartesian slice + all scalars.
% data=load_polar_data()                         % latest published run07 slice
% data=load_polar_data(source,50,[-2 2 0 4],301)  % nearest saved tau to 50
% source may be a run folder, its native folder, or native/data.mat.
% data.x/y are Cartesian axis vectors; data.rho is Ny-by-Nx.
% data.polar retains the selected native Nr-by-Ntheta arrays, without interpolation.
% data.t and data.tau are COMPLETE physical/canonical time trajectories;
% data.c_r has the same length. Selected slice scalars are in data.selected.
% rho is the complete rescaled density R. W=R-R(0) is stored separately.
% This reader performs saved-field interpolation only, with no PDE solve.
project=fileparts(fileparts(mfilename('fullpath')));
if nargin<1||isempty(source)
    control=fullfile(project,'output','polar_q_oscillation_repair_20261021','control.json');
    if isfile(control)
        current=jsondecode(fileread(control));source=current.activeRun;
    else
        source=fullfile(project,'runs','ipm_c_tau200_originfix_20261021_07');
    end
end
if nargin<2,tau=[];end
if nargin<3||isempty(box),box=[-2 2 0 4];end
if nargin<4||isempty(resolution),resolution=201;end
assert(ischar(source)||(isstring(source)&&isscalar(source)), ...
    'polar:PlotSource','source must name a saved run, native folder, or data.mat.');
source=char(source);
if isfolder(source)&&~isfile(fullfile(source,'data.mat'))&& ...
        isfile(fullfile(source,'native','data.mat'))
    source=fullfile(source,'native');
end
addpath(fullfile(project,'time_evolution'),'-begin');
if isempty(tau)
    native=ipm_c.data(source);
else
    native=ipm_c.data(source,tau,'nearest');
end

% Every public field below is reconstructed onto the SAME Cartesian grid.
fields={'rho','rho_centered','rho_X','rho_Y','rho_tau', ...
    'rho_tau_centered','psi','u1','u2','vx','vy','Q'};
fields=fields(cellfun(@(name)isfield(native,name)&& ...
    all(isfinite(native.(name)),'all'),fields));
options=struct('fields',{fields},'coordinates','native','method','cubic');
if isfield(native,'rho_centered')
    options.gradients.rho_centered=struct('X',native.rho_X,'Y',native.rho_Y);
end
data=ipm_c.cartesian(native,box,resolution,options);
data.kind='polar_plot_data_v1';data.schemaVersion=1;
data.rho_x=data.rho_X;data.rho_y=data.rho_Y;
data.F=data.rho_tau;
if isfield(native,'rho_origin'),data.rho_origin=native.rho_origin;
else,data.rho_origin=native.rho(1,1);end
if ~isfield(data,'rho_centered'),data.rho_centered=data.rho-data.rho_origin;end
data.W=data.rho_centered;
if isfield(native,'rho_origin_tau'),data.rho_origin_tau=native.rho_origin_tau;end
if isfield(native,'transportOrigin'),data.transportOrigin=native.transportOrigin;end
data.speed=hypot(data.u1,data.u2);
data.transport_speed=hypot(data.vx,data.vy);
if isfield(native.selected,'diagnostics')&& ...
        isfield(native.selected.diagnostics,'referenceVelocity')
    data.reference_velocity=native.selected.diagnostics.referenceVelocity;
    data.u1_relative=data.u1-data.reference_velocity;
    data.u2_relative=data.u2;
end

% Selecting one spatial slice does not shorten any time trajectory.
for k=1:numel(native.historyColumns)
    name=native.historyColumns{k};data.(name)=native.history(:,k);
end
data.c_l=native.c_l;data.c_omega=native.c_omega;data.c_r=native.c_r;
data.c=[data.c_l,data.c_omega,data.c_r];
data.dtau=native.dt;data.Q_inf=native.Qinf;
data.C_l=native.C_l;data.C_omega=native.C_omega;
if isfield(native,'C_s'),data.C_s=native.C_s;end
data.history=native.history;data.historyColumns=native.historyColumns;
data.fieldTau=native.fieldTau;data.frames=native.frames;
data.polar_grid=native.grid;
% Native polar arrays keep every original node of this frame's actual mesh.
p=struct('r',native.r,'theta',native.theta,'X',native.X,'Y',native.Y, ...
    'grid',native.grid,'selected',native.selected,'rho_origin',data.rho_origin, ...
    'fieldOrientation','Nr-by-Ntheta');
for name={'rho','rho_X','rho_Y','rho_tau','rho_centered', ...
        'rho_tau_centered','Q','psi','u1','u2','vx','vy'}
    if isfield(native,name{1}),p.(name{1})=native.(name{1});end
end
p.rho_x=p.rho_X;p.rho_y=p.rho_Y;
if isfield(p,'rho_centered'),p.W=p.rho_centered;
else,p.W=p.rho-p.rho_origin;end
p.rho_r=p.rho_X.*cos(p.theta)+p.rho_Y.*sin(p.theta);
p.rho_theta=p.r.*(-p.rho_X.*sin(p.theta)+p.rho_Y.*cos(p.theta));
p.ur=p.u1.*cos(p.theta)+p.u2.*sin(p.theta);
p.utheta=-p.u1.*sin(p.theta)+p.u2.*cos(p.theta);
p.vr=p.vx.*cos(p.theta)+p.vy.*sin(p.theta);
p.vtheta=-p.vx.*sin(p.theta)+p.vy.*cos(p.theta);
data.polar=p;
data.grid=struct('x',data.x,'y',data.y,'box',data.box, ...
    'resolution',data.resolution,'coordinates','rescaled Cartesian X,Y');
data.meta=struct('datasetFile',native.meta.datasetFile, ...
    'selectedFrameFile',fullfile(native.meta.baseFolder,native.selected.file), ...
    'requestedTau',tau,'actualTau',native.selected.tau, ...
    'latestHistoryTau',native.tau(end),'nativeSchemaVersion',native.schemaVersion, ...
    'reader',mfilename('fullpath'), ...
    'definitions',struct('x_y','rescaled Cartesian axes X,Y; rho arrays are Ny-by-Nx', ...
    'rho','complete rescaled density R(X,Y,tau)', ...
    'W','rho_centered=R-R(0,tau); rho_origin is the selected scalar R(0,tau)', ...
    'rho_x_y','derivatives dR/dX and dR/dY, not physical-coordinate gradients', ...
    't','complete recorded physical-time vector; late high-part values may repeat', ...
    'tau','complete recorded canonical/rescaled-time vector', ...
    'dtau','accepted canonical time steps; initial entry is NaN', ...
    'c_r','complete recorded rescaled translation coefficient, aligned with t and tau', ...
    'Q','canonical D_X(R_tau) on the selected slice; Q_inf is the native nodal norm trajectory', ...
    'u1_u2','rescaled Darcy velocity', ...
    'u1_relative_u2_relative','Darcy velocity relative to the saved wall-centroid reference velocity', ...
    'vx_vy','complete rescaled transport velocity including dilation/translation'), ...
    'interpolation','existing local P3 finite elements on this frame actual polar mesh; no extrapolation', ...
    'polar','original Nr-by-Ntheta arrays on selected frame mesh; no interpolation', ...
    'polar_rho_theta','dR/dtheta = r*(-sin(theta)*R_X+cos(theta)*R_Y)', ...
    'polar_utheta_vtheta','tangential velocity components, not angular speeds');
assert(isequal(size(data.rho),[numel(data.y),numel(data.x)],size(data.mask)), ...
    'polar:PlotFieldShape','Cartesian axes and fields must match.');
assert(numel(data.t)==numel(data.tau)&&numel(data.t)==numel(data.c_r), ...
    'polar:PlotTimeShape','t, tau and c_r must describe the same full trajectory.');
end
