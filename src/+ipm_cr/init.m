function [s,ctx]=init(opts)
%INIT Original datum and its own prepared simulation context.
if nargin<1,opts=struct();end
[options,cfg]=ipm_cr.config(opts);
g=ipm_cr.mesh.build(options.geometry);
[wallOps,scale,refs,guess]=ipm_cr.model.reference(cfg);
sampled=strcmp(cfg.physics.initialCondition,'sampled_density');
if sampled
    loaded=load(cfg.authorityFile,'frame');
    assert(isfield(loaded,'frame')&&isfield(loaded.frame,'native')&& ...
        isstruct(loaded.frame.native)&&isscalar(loaded.frame.native)&& ...
        all(isfield(loaded.frame.native,{'geometryKey','core','polar','origin'})), ...
        'ipm_cr:SampledNative','A sampled authority must contain its complete native density and geometry identity.');
    native=loaded.frame.native;
    assert(isequaln(native.geometryKey,g.geometryKey),'ipm_cr:SampledGeometry', ...
        'The sampled authority belongs to another transport grid; materialize the initial density for this grid.');
    for patch={'core','polar'}
        key=patch{1};value=native.(key);
        assert(isnumeric(value)&&isreal(value)&&isequal(size(value),size(g.(key).X))&& ...
            all(isfinite(value),'all'),'ipm_cr:SampledNative', ...
            'Sampled native densities must be finite real arrays matching the complete transport patches.');
    end
    a=native.origin;validateattributes(a,{'numeric'},{'scalar','real','finite'});
end
ctx=ipm_cr.context(g,wallOps,options);
assert(strcmp(cfg.manifestSha,ctx.identity.authority.manifestSha256)&& ...
    strcmp(cfg.authoritySha,ctx.identity.authority.sha256), ...
    'ipm_cr:InitialIdentity','Initial datum changed during initialization.');
if sampled
    W=struct('core',native.core-a,'polar',native.polar-a);
else
    density=@(X,Y)ipm_cr.model.initial(struct('X',X,'Y',Y,'ny',size(X,1),'nx',size(X,2)),cfg.physics);
    a=density(0,0);
    W=struct('core',density(g.core.X,g.core.Y)-a,'polar',density(g.polar.X,g.polar.Y)-a);
end
W=ipm_cr.mesh.fill(W,g);W.core(1,g.core.originIndex)=0;
scale.physicalTimeLow=0;scale.densityOrigin=a;
% Keep the public state field order while constructing only one state form.
s=struct('scale',scale,'step',0,'referencePeak',refs.P0, ...
    'originalAuthority',refs,'ctxId',ctx.id,'ctxKey',ctx.key,'W',W);
if ~strcmp(options.background.kind,'none')
    if isequaln(options.background,cfg.background)
        s.eta=struct('core',zeros(size(W.core)),'polar',zeros(size(W.polar)));
    else
        bc=ipm_cr.model.background(options.background,scale,g.core.X,g.core.Y);
        bp=ipm_cr.model.background(options.background,scale,g.polar.X,g.polar.Y);
        s.eta=struct('core',W.core-bc.H,'polar',W.polar-bp.H);
    end
end
s.rhs=struct('c_l',guess(1),'c_r',guess(2));
[s.rhs,s.W]=ipm_cr.rhs(s,ctx);
end
