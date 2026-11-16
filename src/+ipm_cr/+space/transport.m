function [prepared,probe]=transport(rho,g,a,sites)
%PREPARE Freeze spatial reconstruction for one RHS only.
% RESULT=PREPARED(VX,VY,REQUEST) returns named B/rate/position/partials.
% REQUEST defaults to cfl=true, partials=false; no output-count dispatch.
c=g.core;p=g.polar;
check_array(rho.core,size(c.X),'rho.core');
check_array(rho.polar,size(p.X),'rho.polar');
assert(isequal(size(c.active),size(c.X))&&any(c.active,'all')&& ...
    isequal(size(p.active),size(p.X))&&any(p.active,'all'), ...
    'ipm_cr:ExplicitTransportGeometry','Both patches require actual active masks.');
validateattributes(c.dx,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c.dy,{'numeric'},{'scalar','real','finite','positive'});
nr=numel(p.r);nt=numel(p.theta);
assert(isequal(size(p.X),[nr,nt])&&all(p.r>0)&& ...
    numel(p.Jr)==nr&&all(isfinite(p.Jr))&&all(p.Jr>0), ...
    'ipm_cr:ExplicitTransportGeometry','Use a positive-radius mapped annulus.');
angularSpacing=pi/(nt-1);
co=cos(reshape(p.theta,1,[]));si=sin(reshape(p.theta,1,[]));
co([1,end])=[1,-1];si([1,end])=0;
% Center both patches before the shared spline reconstruction.
origin=interp1(c.x,rho.core(1,:),0,'linear');
qc=rho.core-origin;qp=rho.polar-origin;
coreX=spline_pair(a.x,qc');coreY=spline_pair(a.y,qc);
polarRadialOut=spline_pair(a.t,qp);polarTheta=spline_pair(a.theta,qp');
dtMetric=p.r(:);
xSpacing=axis_spacing(c.x).';ySpacing=axis_spacing(c.y);
radialSpacing=p.r(:).*axis_spacing(a.t.nodes);
prepared=@apply;

probe=[];
if nargin>=4
    ix=sites.core;ir=sites.polar;
    coreNormal=mod(ix-1,c.ny)==0;
    polarColumn=floor((ir-1)/nr)+1;
    polarWall=polarColumn==1|polarColumn==nt;
    probeCo=reshape(co(polarColumn),[],1);probeSi=reshape(si(polarColumn),[],1);
    probeR=p.r(mod(ir-1,nr)+1);
    probeX=pair_at(coreX,ix);
    probeY=pair_at(coreY,transpose_index(ix,c.ny,c.nx));
    probeRadial=pair_at(polarRadialOut,transpose_index(ir,nr,nt));
    probeTheta=pair_at(polarTheta,ir);
    probe=@wall_probe;
end

    function [B,partials]=wall_probe(vx,vy)
        uy=vy.core;uy(coreNormal)=0;
        py=vy.polar;py(polarWall)=0;
        vr=vx.polar.*probeCo+py.*probeSi;
        angular=(-vx.polar.*probeSi+py.*probeCo)./probeR;
        angular(polarWall)=0;
        B=struct('core',advection(vx.core,probeX)+advection(uy,probeY), ...
            'polar',advection(vr,probeRadial)./probeR+advection(angular,probeTheta));
        dx=selected(vx.core,probeX);dy=selected(uy,probeY);
        dr=selected(vr,probeRadial)./probeR;da=selected(angular,probeTheta);
        partials=struct('c_l',struct('core',-c.X(ix).*dx-c.Y(ix).*dy,'polar',-probeR.*dr), ...
            'delta',struct('core',-dx,'polar',-probeCo.*dr+probeSi.*da./probeR));
    end

    function result=apply(vx,vy,request)
        if nargin<3,request=struct();end
        request=resolve_request(request);
        rate=[];position=[];partials=[];
        check_array(vx.core,size(c.X),'vx.core');
        check_array(vy.core,size(c.X),'vy.core');
        check_array(vx.polar,size(p.X),'vx.polar');
        check_array(vy.polar,size(p.X),'vy.polar');
        ux=vx.core;uy=vy.core;uy(1,:)=0;
        px=vx.polar;py=vy.polar;py(:,[1,end])=0;
        vr=px.*co+py.*si;
        vt=-px.*si+py.*co;
        angularVelocity=vt./p.r(:);angularVelocity(:,[1,end])=0;
        radial=polarRadialOut;incoming=vr(end,:)'<0;
        outer=outer_boundary(vr(end,:),p.r(end),incoming,request.cfl);
        B=struct('core',advection(ux,coreX)+advection(uy',coreY)', ...
            'polar',advection(vr',radial)'./dtMetric+ ...
                advection(angularVelocity,polarTheta));
        if request.partials
            % Speed derivatives of the frozen sign-selected reconstruction.
            dx=selected(ux,coreX);dy=selected(uy',coreY)';
            dr=selected(vr',radial)'./dtMetric;
            da=selected(angularVelocity,polarTheta);
            partials=struct('c_l',struct('core',-c.X.*dx-c.Y.*dy, ...
                'polar',-p.r(:).*dr), ...
                'delta',struct('core',-dx, ...
                'polar',-co.*dr+si.*da./p.r(:)));
        end
        if request.cfl
            cx=abs(ux)./xSpacing;cy=abs(uy)./ySpacing;
            pr=abs(vr)./radialSpacing;pa=abs(angularVelocity)/angularSpacing;
            pa(:,[1,end])=0;
            coreRate=cx+cy;polarRate=pr+pa;
            coreRate(~c.active)=-Inf;polarRate(~p.active)=-Inf;
            [maximumCore,ic]=max(coreRate(:));
            [maximumPolar,ip]=max(polarRate(:));
            rate=max(maximumCore,maximumPolar);
            if maximumCore>=maximumPolar
                [row,column]=ind2sub(size(c.X),ic);
                position=struct('patch','core','row',row,'column',column, ...
                    'X',c.X(ic),'Y',c.Y(ic),'r',hypot(c.X(ic),c.Y(ic)), ...
                    'theta',atan2(c.Y(ic),c.X(ic)), ...
                    'firstAxisRate',cx(ic),'secondAxisRate',cy(ic), ...
                    'firstAxisSpacing',c.dx,'secondAxisSpacing',c.dy, ...
                    'u1',ux(ic),'u2',uy(ic));
                position.firstAxisSpacing=xSpacing(column);
                position.secondAxisSpacing=ySpacing(row);
            else
                [row,column]=ind2sub(size(p.X),ip);
                position=struct('patch','polar','row',row,'column',column, ...
                    'X',p.X(ip),'Y',p.Y(ip),'r',p.r(row), ...
                    'theta',p.theta(column), ...
                    'firstAxisRate',pr(ip),'secondAxisRate',pa(ip), ...
                    'firstAxisSpacing',radialSpacing(row), ...
                    'secondAxisSpacing',p.r(row)*angularSpacing, ...
                    'u1',px(ip),'u2',py(ip));
            end
            position.maximumCoreRate=maximumCore;
            position.maximumPolarRate=maximumPolar;
            position.reactionIncluded=false;
            position.reconstruction='bspline4';
            position.outerBoundary=outer;
            position.allActiveOverlapIncluded=true;
        end
        result=struct('B',B,'rate',rate,'position',position,'partials',partials, ...
            'outerBoundary',outer);
    end
end

function pair=pair_at(pair,index)
pair=struct('plus',pair.plus(index),'minus',pair.minus(index));
end

function index=transpose_index(index,rows,columns)
index=floor((index-1)/rows)+1+columns*mod(index-1,rows);
end

function pair=spline_pair(axis,field)
[D,Q]=ipm_cr.space.basis.differentiate(axis,field,true);
pair=struct('plus',(D+Q).','minus',(D-Q).');
end
function h=axis_spacing(sites)
step=diff(sites(:));h=min([step;Inf],[Inf;step]);
end

function info=outer_boundary(speed,radius,incoming,complete)
closure='open_quintic_spline_one_sided_volume_derivative';
unresolved=any(incoming);
if unresolved,status='unprovided_physical_inflow';
else,status='no_physical_inflow';end
info=struct('radius',radius,'minimumRadialVelocity',min(speed), ...
    'maximumRadialVelocity',max(speed),'incomingCount',nnz(incoming), ...
    'maximumInwardSpeed',max(0,-min(speed)), ...
    'closure',closure,'physicalInflowDataSupplied',false, ...
    'unresolvedInflow',unresolved,'exteriorTailCertified',false, ...
    'completeEvaluation',complete, ...
    'status',status);
end

function result=advection(velocity,pair)
result=-(max(velocity,0).*pair.plus+min(velocity,0).*pair.minus);
end

function result=selected(velocity,pair)
result=pair.minus;positive=velocity>=0;result(positive)=pair.plus(positive);
end

function check_array(value,shape,name)
assert(isnumeric(value)&&isreal(value)&&isequal(size(value),shape), ...
    'ipm_cr:ExplicitTransportArray', ...
    '%s must be a real native array of the correct shape.',name);
end

function request=resolve_request(overrides)
assert(isstruct(overrides)&&isscalar(overrides)&& ...
    all(ismember(fieldnames(overrides),{'cfl','partials'})), ...
    'ipm_cr:TransportRequest','Use only named logical cfl/partials requests.');
request=struct('cfl',true,'partials',false);
names=fieldnames(overrides);
for k=1:numel(names),request.(names{k})=overrides.(names{k});end
validateattributes(request.cfl,{'logical'},{'scalar'});
validateattributes(request.partials,{'logical'},{'scalar'});
end
