function data = selectTime(data,selector,mode)
%IPM.OUTPUT.SELECTTIME Select a saved native field without interpolation.
%   data=ipm.output.selectTime(data,tau); mesh(data.x,data.y,data.rho_xt)
%   mode='exact' requires a saved time; mode='index' selects a frame number.
if nargin<3 || isempty(mode), mode='nearest'; end
validateattributes(selector,{'numeric'},{'scalar','real','finite'});
mode = validatestring(char(mode),{'nearest','exact','index'});
assert(isfield(data,'schemaVersion') && data.schemaVersion==2, ...
    'ipm:DatasetSchema','Expected a schema-2 native plot dataset.');
data=ipm.output.withResidualNorm(data);
assert(~isempty(data.frames),'ipm:DatasetEmpty','No native spatial frame is available.');
tau = reshape([data.frames.tau],[],1);
switch mode
    case 'index'
        validateattributes(selector,{'numeric'},{'integer','>=',1,'<=',numel(tau)});
        k=selector;
    case 'exact'
        k=find(tau==selector,1,'first');
        if isempty(k)
            error('ipm:DatasetTimeNotSaved','No spatial frame at tau=%.17g. Choose data.fieldTau or use nearest.',selector);
        end
    case 'nearest'
        [~,k]=min(abs(tau-selector));
        if tau(k)~=selector
            warning('ipm:DatasetNearestTime', ...
                'Requested tau=%.17g; selected saved tau=%.17g. No time interpolation.',selector,tau(k));
        end
end
f = ipm.output.loadFrame(data,k);
data.x=f.x(:).'; data.y=f.y(:);
names={'rho','w','rho_t','rho_xt','rho_y','psi','psiLow','u1','u2', ...
    'transport_u1','transport_u2'};
for j=1:numel(names), data.(names{j})=f.(names{j}); end
data.vx=data.transport_u1;
data.vy=data.transport_u2;
data.vmag=hypot(data.vx,data.vy);
data.meta.definitions.compositeVelocity= ...
    ['Saved composite transport velocity in computational X,Y and canonical tau: ' ...
     'vx=U1+cx*X+cr; vy=U2+cy*Y. Aliases of transport_u1/transport_u2; ' ...
     'vmag=hypot(vx,vy).'];
data.rho_x=f.w;
data.psi_low=f.psiLow; data=rmfield(data,'psiLow');
data.physical_x=(data.x-f.Xshift)/f.Cx;
data.physical_y=data.y/f.Cy;
data.rho_physical=f.rho/f.Comega;
data.rho_x_physical=(f.Cx/f.Comega)*f.w;
data.rho_y_physical=(f.Cy/f.Comega)*f.rho_y;
data.u1_physical=f.u1/f.Comega;
data.u2_physical=(f.Cx/f.Cy)*f.u2/f.Comega;
data.psi_physical=f.psi/(f.Cy*f.Comega);
data.psi_low_physical=f.psiLow/(f.Cy*f.Comega);
data.selectedIndex=k; data.selectedTau=f.tau; data.selectedStep=f.step;
data.meta.selection=struct('mode',mode,'requested',selector, ...
    'frameIndex',k,'tau',f.tau,'t',f.t,'step',f.step, ...
    'interpolationUsed',false,'source','native_accepted_state');
end
