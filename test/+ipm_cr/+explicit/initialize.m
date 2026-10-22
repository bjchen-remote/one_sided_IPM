function [state,g,wallOps,options]=initialize(overrides)
%INITIALIZE Analytic original datum on the actual explicit hybrid grids.
if nargin<1,overrides=struct();end
cfg=ipm_c.config();
physics=cfg.physics;
if strcmp(physics.initialCondition,'localized_gaussian_primitive')
    background=struct('kind','gaussian_primitive','amplitude',physics.gaussianAmplitude, ...
        'centerX',physics.gaussianCenterX,'widthX',physics.gaussianWidthX, ...
        'widthY',physics.gaussianWidthY,'skewWeight',physics.gaussianSkewWeight, ...
        'tailCorrection',true);
else,background=struct('kind','none');end
transport=cfg.transport;transport.reconstruction='linear5';
transport.smoothSpeedCoefficient=.1;
options=struct('geometry',struct(),'solverOptions',struct('coreDegree',[48,48], ...
    'angularDegree',64,'polarDegree',16,'rOuter',1e10), ...
    'transport',transport,'cfl',.3,'maxDt',.2,'minDt',1e-10, ...
    'timeIntegrator','rk4','background',background);
names=fieldnames(overrides);
for k=1:numel(names)
    assert(isfield(options,names{k}),'ipm_cr:ExplicitOptions','Unknown explicit option %s.',names{k});
    options.(names{k})=overrides.(names{k});
end
g=ipm_cr.explicit.geometry(options.geometry);
frame=load(cfg.authorityFile,'frame');native=cfg;
native.grid.customX=frame.frame.x;native.grid.customY=frame.frame.y;
authority=struct('rescaling',cfg.scaling,'t0AuthorityFile',cfg.authorityFile);
referenceGrid=ipm_c.mesh.build(cfg.grid);
[cart,~,shared]=ipm_c.gauge.initialAuthority(native,authority,referenceGrid);
wallOps=cart.originalWallOps;
density=@(X,Y)ipm.field.initialDensity(struct('X',X,'Y',Y,'ny',size(X,1),'nx',size(X,2)),cfg.physics);
a=density(0,0);
rho=struct('core',density(g.core.X,g.core.Y)-a,'polar',density(g.polar.X,g.polar.Y)-a);
rho=ipm_cr.explicit.synchronize(rho,g);rho.core(1,g.core.originIndex)=0;
scale=cart.scale;scale.physicalTimeLow=0;scale.densityOrigin=a;
if strcmp(options.background.kind,'none')
    flow=ipm_cr.explicit.evaluate(rho,g,wallOps,options,shared.initialGaugeGuess);
else
    if isequaln(options.background,background)
        residual=struct('core',zeros(size(rho.core)),'polar',zeros(size(rho.polar)));
    else
        bc=ipm_cr.explicit.background(options.background,scale,g.core.X,g.core.Y);
        bp=ipm_cr.explicit.background(options.background,scale,g.polar.X,g.polar.Y);
        residual=struct('core',rho.core-bc.H,'polar',rho.polar-bp.H);
    end
    [flow,rho]=ipm_cr.explicit.evaluateResidual(residual,scale,g,wallOps,options,shared.initialGaugeGuess);
end
state=struct('rho',rho,'scale',scale,'flow',flow,'step',0, ...
    'referencePeak',shared.references.P0,'originalAuthority',shared.references);
if exist('residual','var'),state.residual=residual;end
end
