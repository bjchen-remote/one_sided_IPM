function [state,context]=initialize(cfg)
%INITIALIZE Centered analytic t0; the fixed wall authority is unchanged.
f=load(cfg.authorityFile,'frame');
native=cfg;native.grid.customX=f.frame.x;native.grid.customY=f.frame.y;
carrier=struct('rescaling',cfg.scaling,'t0AuthorityFile',cfg.authorityFile);
g=ipm_c.mesh.build(cfg.grid);
[cart,~,shared]=ipm_c.gauge.initialAuthority(native,carrier,g);
% Keep one authority; no physical or diagnostic arrays in the saved state.
context=struct('g',g,'cart',cart,'shared',shared,'cfg',cfg);
context.cart=rmfield(context.cart,{'wallRho','wallQ'});
cart=context.cart;
ops=ipm_c.spatial.build(cfg,cart,g);
ops.spatialOperators.remesh=@ipm_c.remesh.apply;
rho=ipm.field.initialDensity(struct('X',g.X,'Y',g.Y,'ny',g.nr,'nx',g.ntheta),cfg.physics);
scale=cart.scale;
scale.physicalTimeLow=0;
scale.densityOrigin=rho(1,1);rho=rho-scale.densityOrigin;rho(1,:)=0;
[rate,flow]=ipm.evolve.flow(rho,ops,scale);
state=struct('rho',rho,'flow',flow,'scale',scale,'ops',ops, ...
    'config',cfg,'step',0,'normalizedTime',0, ...
    'runMetadata',struct('architecture',cfg.architecture), ...
    'timeStep',ipm.evolve.emptyTimestep(cfg.time,scale,'initial'), ...
    'mass0',sum((rho+scale.densityOrigin).*ops.integrationWeights,'all'), ...
    'rhoRange0',scale.densityOrigin+[min(rho,[],'all'),max(rho,[],'all')]);
state.rhsCache=ipm.evolve.makeRhsCache(rho,rate,flow,ops,scale);
end
