function report=run_uniform_core_screen(outputPath)
%RUN_UNIFORM_CORE_SCREEN Native initial-state screen for near-uniform cores.
% Full F/Q and the independently known initial Gaussian source are compared.
% Fine axes evaluate the SAME analytic radial map at twice as many logical
% intervals, followed by full metric/derivative/Poisson geometry rebuilding.
% No transfer interpolation, time advance, or spatial qualification occurs.
assert(ischar(outputPath)||isstring(outputPath),'ipm_c:UniformCoreScreenOutput');
outputPath=char(outputPath);parent=fileparts(outputPath);
if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
cfg=ipm_c.config();cfg.grid.nr=321;cfg.grid.ntheta=257;cfg.grid.rmax=1e10;
cfg.grid.stencilWidth=7;cfg.grid.originRings=4;
cfg.runtime.threads=1;cfg.output.verbose=false;cfg.time.cfl=.4;
firstSpacings=[.00996755832309634,.03,.04,.05];
labels={'baseline_sinh','uniform_core_003','uniform_core_004','uniform_core_005'};
savedThreads=maxNumCompThreads(1);
threadCleanup=onCleanup(@()maxNumCompThreads(savedThreads)); %#ok<NASGU>
report=struct('kind','near_uniform_core_native_initial_spatial_screen_v1', ...
    'status','running','completed',false,'passed',false,'threads',1, ...
    'coarseNr',321,'coarseNtheta',257,'fineNr',641,'fineNtheta',513, ...
    'firstSpacings',firstSpacings,'rmax',1e10,'coreRadius',2,'coreTailFraction',.01, ...
    'radialRelativeNodeTolerance',1e-11,'angularAbsoluteNodeTolerance',eps(pi), ...
    'cases',struct([]),'candidateOverBaselineErrorRatios',struct([]), ...
    'timeAdvanced',0,'efficiencyMeasured',false,'accuracyThresholdPreselected',false, ...
    'accuracyQualificationClaimed',false,'referenceConvergenceEstablished',false, ...
    'interpretation',['Native full F/Q/source/velocity at analytic t0; two-level spatial differences ', ...
      'and exact Gaussian source errors are recorded without an accuracy threshold or throughput claim.']);
ipm_c.output.atomicJson(outputPath,report);
for group=1:numel(firstSpacings)
    report.activeCase=labels{group};report.status='running';
    try
        cfg.grid.nr=321;cfg.grid.ntheta=257;cfg.grid.firstSpacing=firstSpacings(group);
        cfg.grid.coreRadius=2;cfg.grid.coreTailFraction=.01;
        if group==1
            cfg.grid.radialMap='sinh';cfg.grid.originFitDegree=4;
        else
            cfg.grid.radialMap='uniform_core';cfg.grid.originFitDegree=6;
        end
        timer=tic;state=ipm_c.initialize(cfg);elapsed=toc(timer);
        g=state.ops.polarGrid;coarse=fields(state);cart=state.ops.cart;scale=state.scale;
        exact=ipm_c_legacy_fixture('background',cfg.physics,g.X,g.Y,scale);
        row=struct('label',labels{group},'firstSpacing',firstSpacings(group), ...
            'coarse',metadata(state,elapsed),'fine',struct(),'matchingNodes',struct(), ...
            'coarseVersusFine',struct(), ...
            'analyticSourceCoarse',region_errors(coarse.source,exact.rx,g.R<=2), ...
            'analyticSourceFine',struct(),'rates',struct(),'wallIPC',struct(), ...
            'completed',false);
        if strcmp(g.radialMap,'uniform_core')
            mapInfo=g.uniformCoreInfo;
            assert(all(isfield(mapInfo,{'a','alpha','power'})), ...
                'ipm_c:UniformCoreScreenMap','The geometry must expose its exact analytic map.');
        else
            mapInfo=struct('alpha',g.radialMapParameter);
        end
        coarseMap=g.radialMap;
        report.cases=append(report.cases,row);
        clear state exact g
        ipm_c.output.atomicJson(outputPath,report);

        s=linspace(0,1,641)';
        if strcmp(coarseMap,'uniform_core')
            rFine=mapInfo.a*s+(cfg.grid.rmax-mapInfo.a)* ...
                (sinh(mapInfo.alpha*s)/sinh(mapInfo.alpha)).^mapInfo.power;
        elseif mapInfo.alpha==0
            rFine=cfg.grid.rmax*s;
        else
            rFine=cfg.grid.rmax*sinh(mapInfo.alpha*s)/sinh(mapInfo.alpha);
        end
        rFine(1)=0;rFine(end)=cfg.grid.rmax;
        fineCfg=cfg;fineCfg.grid.nr=641;fineCfg.grid.ntheta=513;
        fineCfg.grid.firstSpacing=rFine(2);
        % The stub supplies uniform theta and the requested origin basis;
        % every radial field/metric/basis coefficient is rebuilt on rFine.
        stubOptions=fineCfg.grid;stubOptions.radialMap='sinh';
        timer=tic;stub=ipm_c.mesh.build(stubOptions);
        gFine=ipm_c.remesh.radialGeometry(stub,rFine);
        gFine.radialMap=[coarseMap,'_exact_nested_reference'];
        if strcmp(coarseMap,'uniform_core')
            gFine.uniformCoreInfo=mapInfo;
            gFine.uniformCoreInfo.analyticMapParametersSharedWithCoarse=true;
            gFine.uniformCoreInfo.parameterSourceNr=cfg.grid.nr;
            gFine.uniformCoreInfo.parameterSourceFirstSpacing=cfg.grid.firstSpacing;
            gFine.uniformCoreInfo.nr=641;
            gFine.uniformCoreInfo.firstSpacing=rFine(2);
        else
            gFine.radialMapParameter=mapInfo.alpha;
        end
        ops=ipm_c.spatial.build(fineCfg,cart,gFine);
        ops.spatialOperators.remesh=@ipm_c.remesh.apply;
        rho=ipm.field.initialDensity(ops,fineCfg.physics);
        [F,flow]=ipm.evolve.flow(rho,ops,scale);
        cache=ipm.evolve.makeRhsCache(rho,F,flow,ops,scale);
        state=struct('rho',rho,'ops',ops,'flow',flow,'scale',scale, ...
            'rhsCache',cache,'config',fineCfg);
        elapsed=toc(timer);fine=fields(state);
        row.fine=metadata(state,elapsed);
        exact=ipm_c_legacy_fixture('background',cfg.physics,gFine.X,gFine.Y,scale);
        row.analyticSourceFine=region_errors(fine.source,exact.rx,gFine.R<=2);
        radialMismatch=max(abs(fine.r(1:2:end)-coarse.r)./max(abs(coarse.r),realmin));
        angularMismatch=max(abs(fine.theta(1:2:end)-coarse.theta));
        row.matchingNodes=struct('radialRelativeDifferenceInfinity',radialMismatch, ...
            'angularAbsoluteDifferenceInfinity',angularMismatch, ...
            'radialStride',2,'angularStride',2,'interpolationUsed',false, ...
            'sameAnalyticMap',true,'sameOriginalWallAuthority',isequaln(ops.cart.originalWallOps,cart.originalWallOps));
        report.cases(group)=row;ipm_c.output.atomicJson(outputPath,report);
        assert(radialMismatch<=report.radialRelativeNodeTolerance&& ...
            angularMismatch<=report.angularAbsoluteNodeTolerance&& ...
            row.matchingNodes.sameOriginalWallAuthority, ...
            'ipm_c:UniformCoreScreenPhysicalNodes','Direct subsampling requires identical physical nodes and wall authority.');
        core=repmat(coarse.r<=2,1,numel(coarse.theta));
        for key={'F','Q','source','u1','u2'}
            name=key{1};fineField=fine.(name);
            row.coarseVersusFine.(name)=region_errors(coarse.(name),fineField(1:2:end,1:2:end),core);
        end
        velocityError=max(abs(coarse.u1-fine.u1(1:2:end,1:2:end)), ...
            abs(coarse.u2-fine.u2(1:2:end,1:2:end)));
        referenceVelocity=max(abs(fine.u1(1:2:end,1:2:end)),abs(fine.u2(1:2:end,1:2:end)));
        row.coarseVersusFine.velocity=region_error_magnitudes(velocityError,referenceVelocity,core);
        row.rates=struct('names',{{'c_l','c_omega','c_r'}},'coarse',coarse.rates, ...
            'fine',fine.rates,'absoluteDifference',abs(coarse.rates-fine.rates), ...
            'relativeDifference',abs(coarse.rates-fine.rates)./max(abs(fine.rates),realmin));
        row.wallIPC=struct('names',{{'I','P','C'}},'coarse',coarse.IPC, ...
            'fine',fine.IPC,'absoluteDifference',abs(coarse.IPC-fine.IPC));
        row.completed=true;report.cases(group)=row;
        if group>1
            ratio=case_ratios(row,report.cases(1));
            report.candidateOverBaselineErrorRatios=append(report.candidateOverBaselineErrorRatios,ratio);
        end
        clear state ops rho F flow cache stub gFine exact coarse fine fineField core velocityError referenceVelocity cart scale
        ipm_c.output.atomicJson(outputPath,report);
    catch problem
        report.status='failed';report.failure=struct('identifier',problem.identifier,'message',problem.message);
        ipm_c.output.atomicJson(outputPath,report);rethrow(problem);
    end
end
report.status='completed';report.completed=true;report.passed=true;
report=rmfield(report,'activeCase');ipm_c.output.atomicJson(outputPath,report);
end

function value=fields(state)
g=state.ops.polarGrid;f=state.flow;
value=struct('r',g.r,'theta',g.theta,'F',state.rhsCache.rhoRate,'Q',f.fullQ, ...
    'source',f.source,'u1',f.u1,'u2',f.u2,'rates',[f.c_l,f.c_omega,f.c_r], ...
    'IPC',[f.wallIntegral,f.globalWallPeak,f.wallCentroid]);
end

function value=metadata(state,elapsed)
g=state.ops.polarGrid;f=state.flow;a=f.ipmCInfo.gauge;
mapInfo=struct();if isfield(g,'uniformCoreInfo'),mapInfo=g.uniformCoreInfo;end
value=struct('nr',g.nr,'ntheta',g.ntheta,'firstSpacing',g.r(2), ...
    'radialMap',g.radialMap,'uniformCoreInfo',mapInfo, ...
    'originFitDegree',g.originFit.degree,'originFitRadius',g.originFit.scale, ...
    'originFitCondition',g.originFit.rankCondition,'coreSpacing',spacing(g), ...
    'setupSeconds',elapsed,'transportCflRate',f.generatorCflRate,'reactionCflRate',abs(f.c_omega), ...
    'dtCfl',state.config.time.cfl/(f.generatorCflRate+abs(f.c_omega)), ...
    'cflMaximumPosition',f.ipmCInfo.cflPosition,'poisson',f.poissonSolveInfo, ...
    'gauge',struct('rootResidual',a.rootResidual,'auditResidual',a.auditResidual, ...
        'rootTolerance',a.rootTolerance,'auditTolerance',a.auditTolerance, ...
        'rootCondition',a.rootCondition,'iterations',a.iterations,'evaluations',a.evaluations));
end

function result=spacing(g)
result=struct([]);gap=diff(g.r);total=numel(g.r);
for radius=[.5,1,2]
    inside=g.r(2:end)<=radius;d=gap(inside);n=nnz(g.r<=radius);
    if isempty(d),minimum=NaN;maximum=NaN;ratio=NaN;variation=NaN;
    else,minimum=min(d);maximum=max(d);ratio=maximum/minimum;variation=max(abs(d/g.r(2)-1));end
    row=struct('radius',radius,'nodeCountIncludingOrigin',n,'positiveRingCount',n-1, ...
        'nodeFraction',n/total,'completeIntervalCount',numel(d), ...
        'minimumGap',minimum,'maximumGap',maximum,'maximumToMinimumGapRatio',ratio, ...
        'maximumGapRelativeDeviationFromFirstSpacing',variation);
    result=append(result,row);
end
end

function result=region_errors(value,reference,core)
result=region_error_magnitudes(abs(value-reference),abs(reference),core);
end

function result=region_error_magnitudes(errorMagnitude,referenceMagnitude,core)
first=false(size(core));first(1:min(5,size(core,1)),:)=true;
origin=false(size(core));origin(1,:)=true;
result=struct('global',one_region(errorMagnitude,referenceMagnitude,true(size(core))), ...
    'coreRadiusAtMostTwo',one_region(errorMagnitude,referenceMagnitude,core), ...
    'origin',one_region(errorMagnitude,referenceMagnitude,origin), ...
    'originAndFirstFourRings',one_region(errorMagnitude,referenceMagnitude,first));
end

function result=one_region(errorMagnitude,referenceMagnitude,mask)
e=max(errorMagnitude(mask));r=max(referenceMagnitude(mask));
result=struct('absoluteInfinity',e,'referenceInfinity',r,'relativeInfinity',e/max(r,realmin));
end

function result=case_ratios(candidate,baseline)
result=struct('label',candidate.label,'cflStepGain',candidate.coarse.dtCfl/baseline.coarse.dtCfl, ...
    'analyticSourceCoarse',region_ratios(candidate.analyticSourceCoarse,baseline.analyticSourceCoarse), ...
    'analyticSourceFine',region_ratios(candidate.analyticSourceFine,baseline.analyticSourceFine), ...
    'coarseVersusFine',struct(), ...
    'rateAbsoluteDifference',candidate.rates.absoluteDifference./max(baseline.rates.absoluteDifference,realmin));
for key={'F','Q','source','u1','u2','velocity'}
    name=key{1};result.coarseVersusFine.(name)=region_ratios( ...
        candidate.coarseVersusFine.(name),baseline.coarseVersusFine.(name));
end
end

function result=region_ratios(candidate,baseline)
result=struct();
for key={'global','coreRadiusAtMostTwo','origin','originAndFirstFourRings'}
    name=key{1};c=candidate.(name);b=baseline.(name);
    result.(name)=struct('absoluteErrorRatio',c.absoluteInfinity/max(b.absoluteInfinity,realmin), ...
        'relativeErrorRatio',c.relativeInfinity/max(b.relativeInfinity,realmin));
end
end

function result=append(result,row)
if isempty(result),result=row;else,result(end+1)=row;end
end
