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
initial='data/initial_small_source.json';
if isfield(overrides,'initial'),initial=overrides.initial;end
cfg=ipm_cr.model.initial.problem(initial);p=cfg.physics;
sampled=isfield(p,'initialCondition')&&strcmp(p.initialCondition,'sampled_density');
if sampled
    cfg.background=struct('kind','none');
else
    cfg.background=struct('kind','gaussian_primitive','amplitude',p.gaussianAmplitude, ...
        'centerX',p.gaussianCenterX,'widthX',p.gaussianWidthX,'widthY',p.gaussianWidthY, ...
        'skewWeight',p.gaussianSkewWeight,'tailCorrection',true);
    if isfield(p,'gaussianSourceCutoff')
        cfg.background.kind='compact_gaussian_source';
        cfg.background.cutoff=p.gaussianSourceCutoff;
        cfg.background.tailCorrection=false;
    end
end
opts=struct('geometry',default_geometry(),'cfl',.4,'maxDt',.2, ...
    'background',cfg.background,'initial',cfg.manifestFile,'remesh',struct('enabled',true));
opts=merge(opts,overrides);
for name={'cfl','maxDt'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for name={'geometry','background'}
    assert(isstruct(opts.(name{1}))&&isscalar(opts.(name{1})),'ipm_cr:Config','Use scalar module options.');
end
g=opts.geometry;
if ~isfield(g,'interpolationWidth'),g.interpolationWidth=7;end
assert(g.interpolationWidth>=7,'ipm_cr:SplineExchangeOrder', ...
    'The spline method uses at least seven ghost donor points per axis.');
if isfield(g,'staircase')&&g.staircase&&isfield(g,'preset')&& ...
        ismember(string(g.preset),["coarse","medium"])&&~isfield(g,'interfaceAllocation')&& ...
        (~isfield(g,'staircaseFine')||~g.staircaseFine)
    g.interfaceAllocation=[-3,-2,.96,.60,.66,2];
end
opts.geometry=g;
if isstruct(opts.remesh)&&isscalar(opts.remesh)&&~isfield(opts.remesh,'stopTau')&& ...
        isfield(opts.geometry,'preset')&&strcmp(opts.geometry.preset,'coarse')
    opts.remesh.stopTau=50;
end
opts.remesh=ipm_cr.mesh.region.policy(opts.remesh);
if sampled
    assert(isfield(opts.background,'kind')&&strcmp(opts.background.kind,'none'), ...
        'ipm_cr:SampledBackground','A sampled density requires background.kind=''none''; no analytic background is defined.');
end
opts.initial=cfg.manifestFile;
end

function c=case_defaults()
c=struct('format','ipm_cr_case_v1','initial','data/initial_small_source.json', ...
    'domain',struct('rmax',1e10),'method',struct('geometry',default_geometry()),'run',run_defaults(), ...
    'stop',stop_defaults(),'save',struct('data',true,'history',true,'label',''));
end

function g=default_geometry()
g=struct('preset','medium','staircase',true);
end

function r=run_defaults()
r=struct('outputDir',fullfile(pwd,'result'),'tauEnd',200,'wallSeconds',600, ...
    'checkpointEveryTau',1,'checkpointWallSeconds',60,'logEveryTau',1, ...
    'maxSteps',2000000,'resume',true,'verbose',true,'rethrowOnFailure',true,'seed',[]);
end

function s=stop_defaults()
s=struct('physicalTime',[]);
end

function [opts,requested]=run_options(input)
opts=run_defaults();requested=struct();stop=stop_defaults();
save=struct('data',true,'history',true,'label','');
assert(isfield(input,'format')&&strcmp(input.format,'ipm_cr_case_v1'), ...
    'ipm_cr:CaseFormat','Use an ipm_cr_case_v1 case.');
allowed(input,{'format','initial','domain','method','run','stop','save'});
if isfield(input,'run'),opts=merge(opts,input.run);end
if isfield(input,'method'),requested=input.method;end
if isfield(input,'initial')
    assert(~isfield(requested,'initial'),'ipm_cr:Config','Put initial in case.initial only.');
    requested.initial=input.initial;
end
radius=1e10;
if isfield(input,'domain'),allowed(input.domain,{'rmax'});radius=input.domain.rmax;end
validateattributes(radius,{'numeric'},{'scalar','real','finite','positive'});
if ~isfield(requested,'geometry'),requested.geometry=default_geometry();end
if ~isfield(requested.geometry,'grid'),requested.geometry.grid=struct();end
if isfield(requested.geometry.grid,'rmax')
    assert(requested.geometry.grid.rmax==radius,'ipm_cr:DomainRadius','Conflicting domain radii.');
end
requested.geometry.grid.rmax=radius;
if isfield(input,'stop')
    value=input.stop;
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
if ~isempty(opts.seed)
    allowed(opts.seed,{'file','sha256','tau'});
    assert(all(isfield(opts.seed,{'file','sha256','tau'}))&& ...
        (ischar(opts.seed.file)||(isstring(opts.seed.file)&&isscalar(opts.seed.file)))&& ...
        ischar(opts.seed.sha256)&&~isempty(regexp(opts.seed.sha256,'^[0-9a-f]{64}$','once')), ...
        'ipm_cr:Seed','Supply an exact seed file, SHA256 and actual tau.');
    validateattributes(opts.seed.tau,{'numeric'},{'scalar','real','finite','nonnegative'});
    opts.seed.file=char(java.io.File(char(opts.seed.file)).getCanonicalPath());
end
if ~isempty(stop.physicalTime)
    validateattributes(stop.physicalTime,{'numeric'},{'scalar','real','finite','positive'});
end
flag(save.data);flag(save.history);
assert(ischar(save.label)||(isstring(save.label)&&isscalar(save.label)), ...
    'ipm_cr:Config','save.label must be text.');
save.label=char(save.label);
opts.stop=stop;opts.save=save;
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
