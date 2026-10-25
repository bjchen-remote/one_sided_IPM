function [wall,scale,refs,guess]=reference(cfg)
%REFERENCE Original wall authority, independent of transport evolution.
% Read the audited t0 frame once. Rebuild its original Dx and I/P/C without
% constructing Cartesian/polar state carriers or evaluating a PDE RHS.
assert(isstruct(cfg)&&isscalar(cfg)&&all(isfield(cfg, ...
    {'physics','transport','scaling','authorityFile','authoritySha'})), ...
    'polar:T0Configuration','Pass the resolved original-datum configuration.');
p=cfg.physics;
sampled=isstruct(p)&&isfield(p,'initialCondition')&&strcmp(p.initialCondition,'sampled_density');
if ~sampled
    names={'gaussianAmplitude','gaussianCenterX','gaussianWidthX','gaussianWidthY','gaussianSkewWeight'};
    assert(isstruct(p)&&isfield(p,'initialCondition')&& ...
        strcmp(p.initialCondition,'localized_gaussian_primitive')&&all(isfield(p,names)), ...
        'polar:T0GaussianConfiguration','The original Gaussian parameters are required.');
    for j=1:numel(names)
        validateattributes(p.(names{j}),{'numeric'},{'real','finite','scalar'});
    end
    assert(p.gaussianAmplitude>0&&p.gaussianWidthX>0&&p.gaussianWidthY>0&& ...
        p.gaussianSkewWeight>=0&&p.gaussianSkewWeight<=1, ...
        'polar:T0GaussianConfiguration','Invalid Gaussian datum parameters.');
end
assert(strcmp(cfg.transport.spatialDiscretization,'high_order')&&strcmp(p.symmetryMode,'half_plane'), ...
    'polar:T0OriginalSpatialMode','Keep the original high_order half-plane datum.');
fixed={'scalingContract','gaugeTransportAssembly','wallIntegralRadius', ...
    'peakPositionRestoringRate','cOmegaGauge','cOmegaGaugeFloor','pointStrainConditionFloor'};
assert(all(isfield(cfg.scaling,fixed))&& ...
    strcmp(cfg.scaling.scalingContract,'wall_integral_centroid_practical_v1')&& ...
    strcmp(cfg.scaling.gaugeTransportAssembly,'net_velocity_density_centered_advective'), ...
    'polar:T0GaugeSettings','Keep the original complete I/P/C and advective settings.');
assert(isfile(cfg.authorityFile)&&strcmp(ipm_cr.io.sha(cfg.authorityFile),cfg.authoritySha), ...
    'ipm_cr:ReferenceIdentity','The original datum SHA differs from the resolved configuration.');
loaded=load(cfg.authorityFile,'frame');
assert(isfield(loaded,'frame'),'polar:T0AuthorityFrame','The authority file must contain frame.');
f=loaded.frame;
required={'x','y','rho','step','tau','t','Cx','Cy','Comega','Xshift'};
assert(isstruct(f)&&isscalar(f)&&all(isfield(f,required)), ...
    'polar:T0AuthorityFrame','The authority must be a complete original t0 frame.');
assert(isequal([f.step,f.tau,f.t,f.Xshift],[0,0,0,0])&& ...
    isequal([f.Cx,f.Cy,f.Comega],[1,1,1]), ...
    'polar:T0AuthorityClock','The authority must have zero clocks/shift and identity scales.');
x=f.x(:)';y=f.y(:);
assert(numel(x)>=7&&numel(y)>=7&&all(isfinite(x))&&all(diff(x)>0)&& ...
    all(isfinite(y))&&all(diff(y)>0)&&y(1)==0, ...
    'polar:T0OriginalAxes','Original axes must be finite and increasing from the wall.');
[X,Y]=meshgrid(x,y);
if sampled
    assert(isnumeric(f.rho)&&isreal(f.rho)&&isequal(size(f.rho),size(X))&& ...
        all(isfinite(f.rho),'all'),'polar:T0FrameDensity', ...
        'The sampled original density must be a finite real y-by-x array on its measurement axes.');
    analytic=f.rho;
else
    analytic=ipm_cr.model.initial(struct('X',X,'Y',Y,'ny',size(X,1),'nx',size(X,2)),p);
    assert(isequal(size(f.rho),size(analytic))&&all(isfinite(f.rho),'all'), ...
        'polar:T0FrameDensity','The original density has an incompatible shape.');
    assert(max(abs(analytic-f.rho),[],'all')<=256*eps(max(1,max(abs(f.rho),[],'all'))), ...
        'polar:T0AnalyticProvenance','Gaussian parameters do not reproduce original t0 rho.');
end
Dx=ipm_cr.mesh.fd(x,1,7);q=analytic(1,:)*Dx';
if isfield(f,'w')
    assert(max(abs(q-f.w(1,:)))<=256*eps(max(1,max(abs(q)))), ...
        'polar:T0DerivativeProvenance','Rebuilt original wall Dx differs from the t0 source.');
end
settings=cfg.scaling;
runtime={'pinIndex','pinX','originIndex','strainTarget','trackedPeakAnchor', ...
    'gaugeRateGuess','peakEvaluationTime','referenceWallIntegral', ...
    'referenceWallGradientCentroid','referenceWallPeak'};
present=runtime(isfield(settings,runtime));
if ~isempty(present),settings=rmfield(settings,present);end
settings.enabled=true;wall=struct('x',x,'Dx',Dx,'rescaling',settings);
integral=ipm_cr.gauge.integral(q,wall,settings);
peak=ipm_cr.gauge.peak(q,wall);centroid=ipm_cr.gauge.centroid(q,wall);
settings.referenceWallIntegral=integral.value;
settings.referenceWallGradientCentroid=centroid.value;settings.pinX=peak.x;
[~,settings.pinIndex]=min(abs(x-peak.x));[~,settings.originIndex]=min(abs(x));
assert(x(settings.originIndex)==0,'polar:T0Origin','The original wall must contain the exact origin.');
wall.rescaling=settings;
for name={'referenceWallIntegral','referenceWallGradientCentroid'}
    key=name{1};
    if isfield(cfg.scaling,key)
        expected=cfg.scaling.(key);value=settings.(key);
        assert(isnumeric(expected)&&isscalar(expected)&&isfinite(expected)&& ...
            abs(value-expected)<=64*eps(max([1,abs(expected),abs(value)])), ...
            'polar:T0ReferenceMismatch','The original reference differs from persisted %s.',key);
    end
end
scale=struct('logC_l',0,'logC_omega',0,'physicalTime',0,'physicalTimeLow',0, ...
    'X_shift',0,'canonicalTime',0);
refs=struct('I0',integral.value,'P0',peak.value,'C0',centroid.value, ...
    'peakX0',peak.x,'wallIntegralRadius',settings.wallIntegralRadius);
guess=[0,0];
if all(isfield(f,{'transport_u1','u1'}))
    difference=f.transport_u1(1,:)-f.u1(1,:);
    guess(1)=(difference(end)-difference(1))/(x(end)-x(1));
    guess(2)=difference(settings.originIndex);
end
end
