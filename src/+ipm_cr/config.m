function opts=config(overrides)
%CONFIG Resolve the complete numerical method once; PCONFIG is elliptic only.
if nargin<1,overrides=struct();end
cfg=ipm_cr.problem();p=cfg.physics;
background=struct('kind','gaussian_primitive','amplitude',p.gaussianAmplitude, ...
    'centerX',p.gaussianCenterX,'widthX',p.gaussianWidthX,'widthY',p.gaussianWidthY, ...
    'skewWeight',p.gaussianSkewWeight,'tailCorrection',true);
transport=cfg.transport;transport.reconstruction='linear5';transport.smoothSpeedCoefficient=.1;
opts=struct('geometry',struct(),'solverOptions',struct('coreDegree',[48,48], ...
    'angularDegree',64,'polarDegree',16,'rOuter',1e10), ...
    'transport',transport,'cfl',.3,'maxDt',.2,'minDt',1e-10, ...
    'timeIntegrator','rk4','background',background);
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:Config','Use scalar method options.');
names=fieldnames(overrides);
for k=1:numel(names)
    name=names{k};
    % Short input aliases do not introduce a second set of numerical defaults.
    if strcmp(name,'mesh'),target='geometry';elseif strcmp(name,'poisson'),target='solverOptions';else,target=name;end
    assert(isfield(opts,target),'ipm_cr:Config','Unknown numerical option %s.',name);
    if ~strcmp(target,name),assert(~isfield(overrides,target),'ipm_cr:Config','Duplicate option alias %s.',target);end
    opts.(target)=overrides.(name);
end
for name={'cfl','maxDt','minDt'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
assert(opts.minDt<=opts.maxDt,'ipm_cr:Config','minDt must not exceed maxDt.');
opts.timeIntegrator=char(validatestring(opts.timeIntegrator,{'rk4','ssprk3'}));
assert(isstruct(opts.geometry)&&isscalar(opts.geometry)&&isstruct(opts.solverOptions)&& ...
    isscalar(opts.solverOptions)&&isstruct(opts.transport)&&isscalar(opts.transport)&& ...
    isstruct(opts.background)&&isscalar(opts.background),'ipm_cr:Config','Use scalar module options.');
end
