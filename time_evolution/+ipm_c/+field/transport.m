function [B,info]=transport(rho,vx,vy,g,options)
%POLAR_TRANSPORT Complete density-centered mapped WENO advective transport.
% options.disableUnitPair=true forces the unchanged shared reference path.
% B=-V dot grad(rho). No reaction, conservative correction, or PDE mask.
% The joined signed wall supplies the single Cartesian origin value.
% transport(rho,[],[],g,options) prepares a velocity-to-B function handle.
% Its WENO derivatives belong to this rho only, never to a later RHS state.
% Radial WENO receives three Cartesian-consistent
% negative-radius samples; angular and joined-wall operators are unchanged.
if nargin<5,options=struct();end
assert(isequal(size(rho),[g.nr,g.ntheta]),'polar:TransportShape', ...
    'Density must match the complete polar grid.');
assert(all(rho(1,:)==rho(1,1)),'polar:OriginValueMismatch', ...
    'The complete origin row must represent one shared scalar.');
assert(~isfield(options,'poleRegularization')||isequal(options.poleRegularization,false), ...
    'polar:RetiredPoleRegularization','The experimental pole partition has been retired.');
q=rho-rho(1,1);
info=struct('scheme','density_centered_mapped_weno5','singleCartesianOrigin',true, ...
    'sourceOrPdeMask',false);
if isempty(vx) && isempty(vy)
    B=prepare_transport(q,g,options);
    return;
end
assert(isequal(size(vx),size(vy),[g.nr,g.ntheta]),'polar:TransportShape', ...
    'Both Cartesian velocity arrays must match the complete polar grid.');
radialGhost=radial_origin_ghost(q,g);
vr=vx.*g.cosTheta+vy.*g.sinTheta;
angularRate=-vx.*g.sinTheta+vy.*g.cosTheta;
angularRate(2:end,:)=angularRate(2:end,:)./g.r(2:end);
angularRate(1,:)=0;angularRate(:,[1,end])=0;
B=direction(q',vr',g.radialMetric',1/(g.nr-1),'radial')'+ ...
    direction(q,angularRate,1,pi/(g.ntheta-1),'angular');
wallQ=[flipud(q(2:end,end));q(:,1)]';
wallV=[flipud(vx(2:end,end));vx(:,1)]';
wallB=direction(wallQ,wallV,g.signedWallMetric,1/(g.nr-1),'wall');
B(:,1)=wallB(g.nr:end)';B(:,end)=flipud(wallB(1:g.nr)');
B(1,:)=wallB(g.nr);
assert(all(isfinite(B),'all'),'polar:NonfiniteRhs','Complete polar transport must be finite.');

    function result=direction(field,velocity,J,spacing,kind)
        opts=derivative_options(size(field,1),options);
        if strcmp(kind,'radial'),opts.physicalUpperInflow=velocity(:,end)<0;end
        if strcmp(kind,'wall')
            opts.physicalLowerInflow=velocity(:,1)>0;opts.physicalUpperInflow=velocity(:,end)<0;
        end
        if strcmp(kind,'radial')
            pair=direction_derivatives(field,spacing,opts,radialGhost);
        else
            pair=direction_derivatives(field,spacing,opts);
        end
        result=apply_direction(velocity,J,pair);
    end
end

function evaluate=prepare_transport(q,g,options)
% Boundary continuation is the only velocity dependence of these WENO
% derivatives. Preserve the original batch-wide MEX/MATLAB dispatch as well
% as each row's actual inflow continuation.
spacing=1/(g.nr-1);
radialField=q';opts=derivative_options(size(radialField,1),options);
radialGhost=radial_origin_ghost(q,g);
radialOut=direction_derivatives(radialField,spacing,opts,radialGhost);
% Any low-degree inflow disables MEX for the WHOLE original radial batch.
% A zero dummy line forces that same fallback for cached outflow rows;
% weno5FluxDerivative differentiates rows independently. Shared +ipm is
% unchanged. With degree-five continuation MEX remains valid in both cases.
if unit_native_available(opts)&&isa(radialField,'double')
    % The dedicated unit kernel handles each line's boundary degree. No
    % whole-batch MEX/MATLAB dispatch difference or dummy outflow pass exists.
    radialOutMatlab=radialOut;
else
    dummyField=[radialField;zeros(1,size(radialField,2),'like',radialField)];
    dummyOpts=derivative_options(size(dummyField,1),options);
    dummyOpts.physicalUpperInflow(end)=true;
    dummyGhost=[radialGhost;zeros(1,3,'like',radialGhost)];
    radialOutMatlab=direction_derivatives(dummyField,spacing,dummyOpts,dummyGhost);
    radialOutMatlab.plus(end,:)=[];radialOutMatlab.minus(end,:)=[];
end
opts.physicalUpperInflow(:)=true;
radialIn=direction_derivatives(radialField,spacing,opts,radialGhost);
angular=direction_derivatives(q,pi/(g.ntheta-1), ...
    derivative_options(size(q,1),options));
wallField=[flipud(q(2:end,end));q(:,1)]';wall=cell(1,4);
for upper=0:1
    for lower=0:1
        wallOpts=derivative_options(1,options);
        wallOpts.physicalLowerInflow=logical(lower);
        wallOpts.physicalUpperInflow=logical(upper);
        wall{1+lower+2*upper}=direction_derivatives(wallField,spacing,wallOpts);
    end
end
evaluate=@apply;

    function B=apply(vx,vy)
        assert(isequal(size(vx),size(vy),[g.nr,g.ntheta]),'polar:TransportShape', ...
            'Both Cartesian velocity arrays must match the prepared polar grid.');
        vr=vx.*g.cosTheta+vy.*g.sinTheta;
        angularRate=-vx.*g.sinTheta+vy.*g.cosTheta;
        angularRate(2:end,:)=angularRate(2:end,:)./g.r(2:end);
        angularRate(1,:)=0;angularRate(:,[1,end])=0;
        incoming=vr(end,:)'<0;
        if any(incoming)
            radial=radialOutMatlab;
            radial.plus(incoming,:)=radialIn.plus(incoming,:);
            radial.minus(incoming,:)=radialIn.minus(incoming,:);
        else
            radial=radialOut;
        end
        B=apply_direction(vr',g.radialMetric',radial)'+ ...
            apply_direction(angularRate,1,angular);
        wallV=[flipud(vx(2:end,end));vx(:,1)]';
        wallIndex=1+double(wallV(1)>0)+2*double(wallV(end)<0);
        wallB=apply_direction(wallV,g.signedWallMetric,wall{wallIndex});
        B(:,1)=wallB(g.nr:end)';B(:,end)=flipud(wallB(1:g.nr)');
        B(1,:)=wallB(g.nr);
        assert(all(isfinite(B),'all'),'polar:NonfiniteRhs','Complete polar transport must be finite.');
    end
end

function opts=derivative_options(numberOfLines,options)
opts=struct('lowerBoundary','extrapolate','upperBoundary','extrapolate', ...
    'extrapolationDegree',5,'epsilon',1e-12,'alphaMode','local', ...
    'inflowDensityDegree',5,'physicalLowerInflow',false(numberOfLines,1), ...
    'physicalUpperInflow',false(numberOfLines,1));
if isfield(options,'transport')
    tr=options.transport;
    for key={'wenoEpsilon','wenoAlphaMode','wenoInflowDensityDegree'}
        if isfield(tr,key{1})
            switch key{1}
                case 'wenoEpsilon',opts.epsilon=tr.(key{1});
                case 'wenoAlphaMode',opts.alphaMode=tr.(key{1});
                case 'wenoInflowDensityDegree',opts.inflowDensityDegree=tr.(key{1});
            end
        end
    end
end
opts.disableUnitPair=false;
if isfield(options,'disableUnitPair')
    value=options.disableUnitPair;
    assert(isscalar(value)&&(islogical(value)||isnumeric(value))&& ...
        isreal(value)&&isfinite(value)&&any(value==[0,1]), ...
        'ipm_c:UnitWenoOption','disableUnitPair must be a scalar logical.');
    opts.disableUnitPair=logical(value);
end
end

function pair=direction_derivatives(field,spacing,opts,lowerGhost)
% Prepending makes the physical origin an interior WENO node. Three
% supplied lower samples cover every retained derivative's full stencil;
% the kernel's extrapolated samples below these are discarded entirely.
ghostColumns=0;
if nargin>=4
    assert(isequal(size(lowerGhost),[size(field,1),3])&& ...
        all(isfinite(lowerGhost),'all'),'polar:OriginGhostShape', ...
        'Radial origin continuation needs three finite samples per ray.');
    field=[lowerGhost,field];ghostColumns=3;
end
if unit_native_available(opts)&&isa(field,'double')
    [pair.plus,pair.minus]=ipm_c.field.weno5UnitPairMex( ...
        field,spacing,opts.epsilon,opts.inflowDensityDegree, ...
        opts.physicalLowerInflow,opts.physicalUpperInflow);
else
    unit=ones(size(field));
    pair.plus=ipm.field.weno5FluxDerivative(field,unit,1,spacing,opts);
    pair.minus=-ipm.field.weno5FluxDerivative(field,-unit,1,spacing,opts);
end
if ghostColumns
    pair.plus=pair.plus(:,ghostColumns+1:end);
    pair.minus=pair.minus(:,ghostColumns+1:end);
end
end

function ghost=radial_origin_ghost(q,g)
% A negative signed radius on theta has Cartesian coordinates
% (-r*cos(theta),-r*sin(theta)). Read its actual upper-half reflected sample
% (-r*cos(theta),+r*sin(theta)), then subtract twice the odd-in-Y part of
% the existing origin Cartesian jet. Thus the wall rays use the exact
% opposite wall data and every fitted Cartesian polynomial extends with
% its correct signed-ray values. No density or RHS is filtered.
fit=g.originFit;
assert(isfield(fit,'coefficientPowers')&&ismember(fit.degree,[4,6]), ...
    'polar:OriginGhostFit','Use the existing degree-four or degree-six Cartesian origin jet.');
fitValues=q(fit.rows,:)-q(1,1);
coefficients=fit.map*fitValues(:);powers=fit.coefficientPowers;
assert(size(powers,1)==numel(coefficients)&&size(powers,2)==2&& ...
    all(isfinite(coefficients)),'polar:OriginGhostFit', ...
    'The origin jet coefficients must match their Cartesian powers.');
positiveRadii=g.r(4:-1:2)/fit.scale;
x=-positiveRadii*g.cosTheta;y=positiveRadii*g.sinTheta;
oddPart=zeros(size(x),'like',q);
for index=find(mod(powers(:,2),2)==1)'
    oddPart=oddPart+coefficients(index)* ...
        (x.^powers(index,1)).*(y.^powers(index,2));
end
ghost=(q(4:-1:2,end:-1:1)-2*oddPart)';
assert(all(isfinite(ghost),'all'),'polar:OriginGhostFinite', ...
    'Cartesian signed-ray continuation must remain finite.');
end

function available=unit_native_available(opts)
% Look up each time: compiling the optional kernel needs no session restart.
% The original shared path remains selectable and is the missing-MEX fallback.
available=false;
if opts.disableUnitPair||~strcmp(opts.alphaMode,'local'),return;end
nativeFile=which('ipm_c.field.weno5UnitPairMex');
available=~isempty(nativeFile)&&exist(nativeFile,'file')==3;
end

function result=apply_direction(velocity,J,pair)
result=-(max(velocity,0).*pair.plus+min(velocity,0).*pair.minus)./J;
end
