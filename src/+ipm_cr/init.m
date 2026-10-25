function [s,ctx]=init(opts)
%INIT Original datum and its own prepared simulation context.
if nargin<1,opts=struct();end
cfg=ipm_cr.problem();
physics=cfg.physics;
if strcmp(physics.initialCondition,'localized_gaussian_primitive')
    background=struct('kind','gaussian_primitive','amplitude',physics.gaussianAmplitude, ...
        'centerX',physics.gaussianCenterX,'widthX',physics.gaussianWidthX, ...
        'widthY',physics.gaussianWidthY,'skewWeight',physics.gaussianSkewWeight, ...
        'tailCorrection',true);
else,background=struct('kind','none');end
options=ipm_cr.config(opts);
g=ipm_cr.mesh.build(options.geometry);
[wallOps,scale,refs,guess]=ipm_cr.model.reference(cfg);
ctx=ipm_cr.context(g,wallOps,options);
density=@(X,Y)ipm_cr.model.initial(struct('X',X,'Y',Y,'ny',size(X,1),'nx',size(X,2)),cfg.physics);
a=density(0,0);
W=struct('core',density(g.core.X,g.core.Y)-a,'polar',density(g.polar.X,g.polar.Y)-a);
W=ipm_cr.mesh.fill(W,g);W.core(1,g.core.originIndex)=0;
scale.physicalTimeLow=0;scale.densityOrigin=a;
% Keep the public state field order while constructing only one state form.
s=struct('scale',scale,'step',0,'referencePeak',refs.P0, ...
    'originalAuthority',refs,'ctxId',ctx.id,'ctxKey',ctx.key,'W',W);
if ~strcmp(options.background.kind,'none')
    if isequaln(options.background,background)
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
