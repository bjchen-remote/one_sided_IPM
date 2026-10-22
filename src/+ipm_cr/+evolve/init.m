function [state,g,wallOps,options,ctx]=init(overrides)
%INITIALIZE Analytic original datum on the actual explicit hybrid grids.
if nargin<1,overrides=struct();end
cfg=ipm_cr.problem();
physics=cfg.physics;
if strcmp(physics.initialCondition,'localized_gaussian_primitive')
    background=struct('kind','gaussian_primitive','amplitude',physics.gaussianAmplitude, ...
        'centerX',physics.gaussianCenterX,'widthX',physics.gaussianWidthX, ...
        'widthY',physics.gaussianWidthY,'skewWeight',physics.gaussianSkewWeight, ...
        'tailCorrection',true);
else,background=struct('kind','none');end
options=ipm_cr.config(overrides);
g=ipm_cr.mesh.build(options.geometry);
frame=load(cfg.authorityFile,'frame');native=cfg;
native.grid.customX=frame.frame.x;native.grid.customY=frame.frame.y;
authority=struct('rescaling',cfg.scaling,'t0AuthorityFile',cfg.authorityFile);
referenceGrid=ipm_cr.mesh.polar(cfg.grid);
[cart,~,shared]=ipm_cr.model.reference(native,authority,referenceGrid);
wallOps=cart.originalWallOps;
ctx=ipm_cr.context(g,wallOps,options);
density=@(X,Y)ipm_cr.model.initial(struct('X',X,'Y',Y,'ny',size(X,1),'nx',size(X,2)),cfg.physics);
a=density(0,0);
rho=struct('core',density(g.core.X,g.core.Y)-a,'polar',density(g.polar.X,g.polar.Y)-a);
rho=ipm_cr.mesh.fill(rho,g);rho.core(1,g.core.originIndex)=0;
scale=cart.scale;scale.physicalTimeLow=0;scale.densityOrigin=a;
if strcmp(options.background.kind,'none')
    flow=ipm_cr.evolve.rhs(rho,g,wallOps,options,shared.initialGaugeGuess,[],true,ctx);
else
    if isequaln(options.background,background)
        residual=struct('core',zeros(size(rho.core)),'polar',zeros(size(rho.polar)));
    else
        bc=ipm_cr.model.background(options.background,scale,g.core.X,g.core.Y);
        bp=ipm_cr.model.background(options.background,scale,g.polar.X,g.polar.Y);
        residual=struct('core',rho.core-bc.H,'polar',rho.polar-bp.H);
    end
    [flow,rho]=ipm_cr.evolve.residual(residual,scale,g,wallOps,options,shared.initialGaugeGuess,true,ctx);
end
state=struct('rho',rho,'scale',scale,'flow',flow,'step',0, ...
    'referencePeak',shared.references.P0,'originalAuthority',shared.references);
if exist('residual','var'),state.residual=residual;end
state.ctxId=ctx.id;state.ctxKey=ctx.key;
end
