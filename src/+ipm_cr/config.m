function [opts,extra,requested]=config(overrides,mode)
%CONFIG One resolver for method options and portable experiment cases.
% config() retains the method API. config('case') returns an editable case.
% config(case,'run') resolves its run controls and method, without a PDE solve.
if nargin<1,overrides=struct();end
if nargin<2,mode='method';end
if (ischar(overrides)||isstring(overrides))&&strcmp(overrides,'case')
    opts=case_defaults();extra=[];requested=[];return
end
assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:Config','Use a scalar option structure.');
if strcmp(mode,'run')
    [opts,requested]=run_options(overrides);
    [extra,~]=method_options(requested);
else
    assert(strcmp(mode,'method'),'ipm_cr:Config','Unknown configuration mode.');
    [opts,extra]=method_options(overrides);requested=[];
end
end

function [opts,cfg]=method_options(overrides)
initial='data/initial.json';
if isfield(overrides,'initial'),initial=overrides.initial;end
cfg=ipm_cr.problem(initial);p=cfg.physics;
sampled=isfield(p,'initialCondition')&&strcmp(p.initialCondition,'sampled_density');
if sampled
    cfg.background=struct('kind','none');
else
    cfg.background=struct('kind','gaussian_primitive','amplitude',p.gaussianAmplitude, ...
        'centerX',p.gaussianCenterX,'widthX',p.gaussianWidthX,'widthY',p.gaussianWidthY, ...
        'skewWeight',p.gaussianSkewWeight,'tailCorrection',true);
end
transport=cfg.transport;transport.reconstruction='linear5';transport.smoothSpeedCoefficient=.1;
transport.outerInflowPolicy='error';transport.outerInflowTolerance=0;
opts=struct('geometry',struct(),'solverOptions',struct('coreDegree',[48,48], ...
    'angularDegree',64,'polarDegree',16,'rOuter',1e10), ...
    'transport',transport,'cfl',.3,'maxDt',.2,'minDt',1e-10, ...
    'timeIntegrator','rk4','background',cfg.background,'initial',cfg.manifestFile);
overrides=aliases(overrides);
names=fieldnames(overrides);
for k=1:numel(names)
    name=names{k};assert(isfield(opts,name),'ipm_cr:Config','Unknown numerical option %s.',name);
    opts.(name)=overrides.(name);
end
for name={'cfl','maxDt','minDt'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
assert(opts.minDt<=opts.maxDt,'ipm_cr:Config','minDt must not exceed maxDt.');
opts.timeIntegrator=char(validatestring(opts.timeIntegrator,{'rk4','ssprk3'}));
for name={'geometry','solverOptions','transport','background'}
    assert(isstruct(opts.(name{1}))&&isscalar(opts.(name{1})),'ipm_cr:Config','Use scalar module options.');
end
% A partial transport override must not silently disable the production gate.
if ~isfield(opts.transport,'outerInflowPolicy'),opts.transport.outerInflowPolicy='error';end
if ~isfield(opts.transport,'outerInflowTolerance'),opts.transport.outerInflowTolerance=0;end
if sampled
    assert(isfield(opts.background,'kind')&&strcmp(opts.background.kind,'none'), ...
        'ipm_cr:SampledBackground','A sampled density requires background.kind=''none''; no analytic background is defined.');
end
if isfield(opts.geometry,'grid'),retire_origin(opts.geometry.grid);end
if isfield(opts.geometry,'grid')&&isfield(opts.geometry.grid,'rmax')&&isfield(opts.solverOptions,'rOuter')
    assert(opts.geometry.grid.rmax==opts.solverOptions.rOuter,'ipm_cr:DomainRadius', ...
        'Transport rmax and Poisson rOuter must agree. Set case.domain.rmax once.');
end
opts.initial=cfg.manifestFile;
end

function c=case_defaults()
c=struct('format','ipm_cr_case_v1','initial','data/initial.json', ...
    'domain',struct('rmax',1e10),'method',struct(),'run',run_defaults(), ...
    'stop',stop_defaults(),'save',struct('data',true,'history',true,'label',''));
end

function r=run_defaults()
r=struct('outputDir',fullfile(pwd,'result'),'tauEnd',200,'wallSeconds',600, ...
    'checkpointEveryTau',1,'checkpointWallSeconds',60,'logEveryTau',1, ...
    'maxSteps',2000000,'resume',true,'verbose',true,'rethrowOnFailure',true);
end

function s=stop_defaults()
s=struct('physicalTime',[],'steady',struct('enabled',false,'residualTol',1e-6, ...
    'changeTol',1e-5,'windowTau',10,'windows',3,'sampleEveryTau',1, ...
    'minTau',0,'nearRadius',3,'bodyRadius',200));
end

function [opts,requested]=run_options(input)
opts=run_defaults();requested=struct();stop=stop_defaults();
save=struct('data',true,'history',true,'label','');
isCase=isfield(input,'format');
if isCase
    assert(strcmp(input.format,'ipm_cr_case_v1'),'ipm_cr:CaseFormat','Unknown experiment case format.');
    allowed(input,{'format','initial','domain','method','run','stop','save'});
    if isfield(input,'run'),opts=merge(opts,input.run);end
    if isfield(input,'method'),requested=aliases(input.method);end
    if isfield(input,'initial')
        assert(~isfield(requested,'initial'),'ipm_cr:Config','Put initial in case.initial only.');
        requested.initial=input.initial;
    end
    radius=1e10;
    if isfield(input,'domain'),allowed(input.domain,{'rmax'});radius=input.domain.rmax;end
    validateattributes(radius,{'numeric'},{'scalar','real','finite','positive'});
    if ~isfield(requested,'geometry'),requested.geometry=struct();end
    if ~isfield(requested.geometry,'grid'),requested.geometry.grid=struct();end
    if ~isfield(requested,'solverOptions'),requested.solverOptions=struct();end
    if isfield(requested.geometry.grid,'rmax')
        assert(requested.geometry.grid.rmax==radius,'ipm_cr:DomainRadius','Conflicting domain radii.');
    end
    if isfield(requested.solverOptions,'rOuter')
        assert(requested.solverOptions.rOuter==radius,'ipm_cr:DomainRadius','Conflicting domain radii.');
    end
    requested.geometry.grid.rmax=radius;requested.solverOptions.rOuter=radius;
else
    % Compatibility exists at this one boundary, never in run/init loops.
    names=fieldnames(input);method={'geometry','solverOptions','transport','cfl', ...
        'maxDt','minDt','timeIntegrator','background','initial','mesh','poisson'};
    if isfield(input,'explicitOptions'),requested=aliases(input.explicitOptions);end
    direct=struct();
    for k=1:numel(names)
        name=names{k};
        if strcmp(name,'explicitOptions')||ismember(name,{'stop','save'}),continue;end
        if ismember(name,method),direct.(name)=input.(name);
        else
            assert(isfield(opts,name),'ipm_cr:Config','Unknown run option %s.',name);
            opts.(name)=input.(name);
        end
    end
    direct=aliases(direct);names=fieldnames(direct);
    for k=1:numel(names)
        name=names{k};assert(~isfield(requested,name),'ipm_cr:Config','Duplicate numerical option %s.',name);
        requested.(name)=direct.(name);
    end
end
if isfield(input,'stop')
    value=input.stop;
    if isfield(value,'steady'),stop.steady=merge(stop.steady,value.steady);value=rmfield(value,'steady');end
    stop=merge(stop,value);
end
if isfield(input,'save'),save=merge(save,input.save);end
assert((ischar(opts.outputDir)&&isrow(opts.outputDir))|| ...
    (isstring(opts.outputDir)&&isscalar(opts.outputDir)),'ipm_cr:Config','outputDir must be text.');
opts.outputDir=char(opts.outputDir);
assert(~isempty(opts.outputDir),'ipm_cr:RunOutput','Choose a nonempty output directory.');
if ~java.io.File(opts.outputDir).isAbsolute(),opts.outputDir=fullfile(pwd,opts.outputDir);end
for name={'tauEnd','wallSeconds','checkpointEveryTau','checkpointWallSeconds','logEveryTau'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(opts.maxSteps,{'numeric'},{'scalar','real','finite','integer','positive'});
for name={'resume','verbose','rethrowOnFailure'},flag(opts.(name{1}));end
if ~isempty(stop.physicalTime)
    validateattributes(stop.physicalTime,{'numeric'},{'scalar','real','finite','positive'});
end
st=stop.steady;flag(st.enabled);
for name={'residualTol','changeTol','windowTau','sampleEveryTau','nearRadius','bodyRadius'}
    validateattributes(st.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(st.windows,{'numeric'},{'scalar','real','finite','integer','positive'});
validateattributes(st.minTau,{'numeric'},{'scalar','real','finite','nonnegative'});
assert(st.bodyRadius>=st.nearRadius,'ipm_cr:StopRegions','bodyRadius must cover nearRadius.');
flag(save.data);flag(save.history);
assert(ischar(save.label)||(isstring(save.label)&&isscalar(save.label)), ...
    'ipm_cr:Config','save.label must be text.');
save.label=char(save.label);
opts.stop=stop;opts.save=save;
end

function s=aliases(s)
assert(isstruct(s)&&isscalar(s),'ipm_cr:Config','Use scalar numerical options.');
for pair={{'mesh','geometry'},{'poisson','solverOptions'}}
    p=pair{1};
    if isfield(s,p{1})
        assert(~isfield(s,p{2}),'ipm_cr:Config','Duplicate option alias %s.',p{2});
        s.(p{2})=s.(p{1});s=rmfield(s,p{1});
    end
end
end

function out=merge(out,input)
allowed(input,fieldnames(out));names=fieldnames(input);
for k=1:numel(names),out.(names{k})=input.(names{k});end
end

function allowed(s,names)
assert(isstruct(s)&&isscalar(s),'ipm_cr:Config','Use a scalar configuration block.');
unknown=setdiff(fieldnames(s),names);
assert(isempty(unknown),'ipm_cr:Config','Unknown configuration fields: %s.',strjoin(unknown,', '));
end

function flag(value)
validateattributes(value,{'logical'},{'scalar'});
end

function retire_origin(g)
for name={'originRings','originFitDegree'}
    if isfield(g,name{1})
        error('ipm_cr:RetiredOption','%s no longer affects the origin-free annulus; remove it.',name{1});
    end
end
end
