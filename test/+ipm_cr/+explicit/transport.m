function [result,rate,position,partials]=transport(rho,varargin)
%TRANSPORT Explicit fifth-order advection on both overset transport patches.
%   PREPARED=TRANSPORT(RHO,G,OPTIONS) caches the unit-speed left/right
%   derivatives of this density only. Then [B,RATE,POSITION]=PREPARED(VX,VY)
%   applies the actual net Cartesian velocities, stored in core/polar arrays.
%   PREPARED(VX,VY,false) skips only the discarded CFL reduction in gauge
%   probes and returns empty RATE/POSITION. Real stages use the true default.
%   [B,RATE,POSITION]=TRANSPORT(RHO,VX,VY,G,OPTIONS) is the direct equivalent.
%   OPTIONS may be the transport options themselves or contain .transport.
%
%   Both patch interiors evolve, including their overlap. Only algebraic
%   ghost collars marked inactive by G are excluded from the CFL maximum.
%   With smoothSpeedCoefficient=nu>0, use alpha=hypot(v,nu*h_core)
%   in a local LF split and in the full mapped CFL. Physical-wall normal
%   transport is exactly zero. The reaction rate is added by the caller.
%   reconstruction: legacy_weno5, smooth_weno5 or linear5. The latter two
%   require build_ipm_cr_native and use uniform degree-five outer ghosts.
%   No small-velocity cutoff or CFL exemption is applied.

if numel(varargin)==1||numel(varargin)==2
    g=varargin{1};options=struct();
    if numel(varargin)==2,options=varargin{2};end
    result=prepare(rho,g,options);rate=[];position=[];partials=[];
elseif numel(varargin)==3||numel(varargin)==4
    vx=varargin{1};vy=varargin{2};g=varargin{3};options=struct();
    if numel(varargin)==4,options=varargin{4};end
    prepared=prepare(rho,g,options);
    if nargout>3,[result,rate,position,partials]=prepared(vx,vy);
    else,[result,rate,position]=prepared(vx,vy);partials=[];end
else
    error('ipm_cr:ExplicitTransportArguments', ...
        'Use transport(rho,g,options) or transport(rho,vx,vy,g,options).');
end
end

function evaluate=prepare(rho,g,options)
c=g.core;p=g.polar;
tr=options;if isfield(options,'transport'),tr=options.transport;end
assert(isstruct(tr)&&isscalar(tr),'ipm_cr:ExplicitTransportOptions');
nu=.1;if isfield(tr,'smoothSpeedCoefficient'),nu=tr.smoothSpeedCoefficient;end
validateattributes(nu,{'numeric'},{'scalar','real','finite','nonnegative'});
% One physical speed floor across patches: it shrinks with core refinement
% even in the logarithmic far field. The complete CFL includes this floor.
sigma=nu*min(c.dx,c.dy);
reconstruction='linear5';
if isfield(tr,'reconstruction'),reconstruction=char(tr.reconstruction);end
highOrderOuter=sigma>0||~strcmp(reconstruction,'legacy_weno5');
% The selected mode keeps one outer extrapolation operator for both speed
% signs. This is a numerical closure at finite R, not exact inflow data;
% accuracy still requires the separate exterior-tail/R audit.
assert(isstruct(options)&&isscalar(options),'ipm_cr:ExplicitTransportOptions', ...
    'Use a scalar transport options struct.');
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
radialMetric=p.Jr(:);radialSpacing=radialMetric/(nr-1);
angularSpacing=pi/(nt-1);
co=cos(reshape(p.theta,1,[]));si=sin(reshape(p.theta,1,[]));
co([1,end])=[1,-1];si([1,end])=0;
% The same origin value is subtracted in both patches before the nonlinear
% WENO weights are formed. No negative-radius extension is needed because
% the annulus's inner ghost collar already contains Cartesian donor data.
origin=interp1(c.x,rho.core(1,:),0,'linear');
assert(isfinite(origin),'ipm_cr:ExplicitTransportOrigin', ...
    'The Cartesian transport patch must contain the physical origin.');
qc=rho.core-origin;qp=rho.polar-origin;
coreX=unit_pair(qc,c.dx,options,false(size(qc,1),1),false(size(qc,1),1));
coreY=unit_pair(qc',c.dy,options,false(size(qc,2),1),false(size(qc,2),1));
polarTheta=unit_pair(qp,angularSpacing,options,false(nr,1),false(nr,1));
[polarRadialOut,radialNative]=unit_pair(qp',1/(nr-1),options,false(nt,1),false(nt,1));
polarRadialMixedOut=polarRadialOut;
if ~highOrderOuter&&~radialNative
    % The maintained shared fallback disables its local MEX for a whole
    % batch if any line uses low-degree inflow data. A zero extra line
    % forces that same dispatch while preserving each real outflow line.
    dummyField=[qp';zeros(1,nr,'like',qp)];
    upperMask=false(nt+1,1);upperMask(end)=true;
    polarRadialMixedOut=unit_pair(dummyField,1/(nr-1),options,false(nt+1,1),upperMask);
    polarRadialMixedOut.plus(end,:)=[];
    polarRadialMixedOut.minus(end,:)=[];
end
polarRadialIn=[];
if ~highOrderOuter
    polarRadialIn=unit_pair(qp',1/(nr-1),options,false(nt,1),true(nt,1));
end
sy=sigma*ones(size(qc'));sy(:,1)=0;
sa=repmat(sigma./p.r(:),1,nt);sa(:,[1,end])=0;
evaluate=@apply;

    function [B,rate,position,partials]=apply(vx,vy,computeCfl)
        if nargin<3,computeCfl=true;end
        validateattributes(computeCfl,{'logical'},{'scalar'});
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
        if ~highOrderOuter&&any(incoming)
            radial=polarRadialMixedOut;
            radial.plus(incoming,:)=polarRadialIn.plus(incoming,:);
            radial.minus(incoming,:)=polarRadialIn.minus(incoming,:);
        end
        alphaX=[];alphaY=[];alphaR=[];alphaTheta=[];
        if sigma>0
            % Reuse one speed magnitude in the LF split, its exact speed
            % derivative, and the complete physical/mapped CFL reduction.
            alphaX=hypot(ux,sigma);alphaY=hypot(uy',sy);
            alphaR=hypot(vr',sigma);alphaTheta=hypot(angularVelocity,sa);
            alphaY(:,1)=0;alphaTheta(:,[1,end])=0;
        end
        B=struct('core',advection(ux,coreX,alphaX)+advection(uy',coreY,alphaY)', ...
            'polar',advection(vr',radial,alphaR)'./radialMetric+ ...
                advection(angularVelocity,polarTheta,alphaTheta));
        assert(all(isfinite(B.core),'all')&&all(isfinite(B.polar),'all'), ...
            'ipm_cr:ExplicitTransportNonfinite','Both complete advection arrays must remain finite.');
        if nargout>3
            % Exact speed derivatives of the frozen reconstruction. For
            % smooth LF, D_eff=.5*(Dp+Dm+(v/alpha)*(Dp-Dm)); the legacy
            % zero-floor path retains its selected upwind derivative.
            % The gauge line search and full RHS audit remain binding.
            dx=selected(ux,coreX,alphaX);dy=selected(uy',coreY,alphaY)';
            dr=selected(vr',radial,alphaR)'./radialMetric;
            da=selected(angularVelocity,polarTheta,alphaTheta);
            partials=struct('c_l',struct('core',-c.X.*dx-c.Y.*dy, ...
                'polar',-p.r(:).*dr), ...
                'delta',struct('core',-dx, ...
                'polar',-co.*dr+si.*da./p.r(:)));
        end
        if nargout>1&&computeCfl
            if sigma>0
                cx=alphaX/c.dx;cy=alphaY'/c.dy;
                pr=alphaR'./radialSpacing;pa=alphaTheta/angularSpacing;
            else
                cx=abs(ux)/c.dx;cy=abs(uy)/c.dy;
                pr=abs(vr)./radialSpacing;pa=abs(angularVelocity)/angularSpacing;
            end
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
            position.smoothPhysicalSpeedFloor=sigma;
            position.reconstruction=reconstruction;
            position.outerDegreeSwitching=~highOrderOuter;
            position.allActiveOverlapIncluded=true;
            assert(isfinite(rate)&&rate>=0,'ipm_cr:ExplicitTransportCfl', ...
                'The complete active physical CFL rate must be finite and nonnegative.');
        end
    end
end

function [pair,useNative]=unit_pair(field,spacing,options,lowerIncoming,upperIncoming)
opts=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate', ...
    'extrapolationDegree',5,'inflowDensityDegree',1,'epsilon',1e-12, ...
    'alphaMode','local','physicalLowerInflow',lowerIncoming, ...
    'physicalUpperInflow',upperIncoming);
tr=options;if isfield(options,'transport'),tr=options.transport;end
assert(isstruct(tr)&&isscalar(tr),'ipm_cr:ExplicitTransportOptions', ...
    'The transport settings must be a scalar struct.');
for name={'wenoEpsilon','wenoAlphaMode','wenoInflowDensityDegree'}
    if isfield(tr,name{1})
        switch name{1}
            case 'wenoEpsilon',opts.epsilon=tr.(name{1});
            case 'wenoAlphaMode',opts.alphaMode=tr.(name{1});
            case 'wenoInflowDensityDegree',opts.inflowDensityDegree=tr.(name{1});
        end
    end
end
disableNative=false;
if isfield(options,'disableUnitPair'),disableNative=options.disableUnitPair;end
assert(isscalar(disableNative)&&(islogical(disableNative)||isnumeric(disableNative)) ...
    &&isreal(disableNative)&&isfinite(disableNative)&&any(disableNative==[0,1]), ...
    'ipm_cr:ExplicitTransportOptions','disableUnitPair must be a scalar logical flag.');
reconstruction='linear5';
if isfield(tr,'reconstruction'),reconstruction=char(tr.reconstruction);end
assert(any(strcmp(reconstruction,{'legacy_weno5','smooth_weno5','linear5'})), ...
    'ipm_cr:ExplicitReconstruction','Use legacy_weno5, smooth_weno5, or linear5.');
if ~strcmp(reconstruction,'legacy_weno5')
    nativeFile=which('ipm_cr.explicit.wenoPairSmoothMex');
    assert(~disableNative&&strcmp(opts.alphaMode,'local')&&isa(field,'double')&& ...
        ~isempty(nativeFile)&&exist(nativeFile,'file')==3, ...
        'ipm_cr:ExplicitNativeRequired','Run build_ipm_cr_native before using the smooth or linear reconstruction.');
    mode=double(strcmp(reconstruction,'smooth_weno5'));
    [pair.plus,pair.minus]=ipm_cr.explicit.wenoPairSmoothMex( ...
        field,spacing,opts.epsilon,opts.inflowDensityDegree, ...
        opts.physicalLowerInflow,opts.physicalUpperInflow,mode);
    useNative=true;return
end
nativeFile=which('ipm_c.field.weno5UnitPairMex');
useNative=~disableNative&&strcmp(opts.alphaMode,'local')&&isa(field,'double') ...
    &&~isempty(nativeFile)&&exist(nativeFile,'file')==3;
if useNative
    [pair.plus,pair.minus]=ipm_c.field.weno5UnitPairMex( ...
        field,spacing,opts.epsilon,opts.inflowDensityDegree, ...
        opts.physicalLowerInflow,opts.physicalUpperInflow);
else
    unit=ones(size(field),'like',field);
    pair.plus=ipm.field.weno5FluxDerivative(field,unit,1,spacing,opts);
    pair.minus=-ipm.field.weno5FluxDerivative(field,-unit,1,spacing,opts);
end
end

function result=advection(velocity,pair,alpha)
if isempty(alpha)
    result=-(max(velocity,0).*pair.plus+min(velocity,0).*pair.minus);
else
    result=-.5*(velocity.*(pair.plus+pair.minus)+ ...
        alpha.*(pair.plus-pair.minus));
end
end

function result=selected(velocity,pair,alpha)
if isempty(alpha)
    result=pair.minus;positive=velocity>=0;result(positive)=pair.plus(positive);
else
    ratio=zeros(size(velocity));
    nonzero=alpha>0;ratio(nonzero)=velocity(nonzero)./alpha(nonzero);
    result=.5*((pair.plus+pair.minus)+ratio.*(pair.plus-pair.minus));
end
end

function check_array(value,shape,name)
assert(isnumeric(value)&&isreal(value)&&isequal(size(value),shape) ...
    &&all(isfinite(value),'all'),'ipm_cr:ExplicitTransportArray', ...
    '%s must be a finite real native array of the correct shape.',name);
end
