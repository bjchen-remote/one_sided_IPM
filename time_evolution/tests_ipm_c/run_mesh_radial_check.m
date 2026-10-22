function report=run_mesh_radial_check(outputPath)
%RUN_MESH_RADIAL_CHECK Native F/Q on exactly nested sinh/angle refinements.
% Each coarse axis keeps its sinh parameter when both interval counts double.
% Comparisons use direct subsampling at the same physical points; no transfer
% interpolation and no time advance occur. This is a two-level error screen,
% not a spatial convergence or trajectory qualification.
assert(ischar(outputPath)||isstring(outputPath),'ipm_c:RadialCheckOutput');
outputPath=char(outputPath);parent=fileparts(outputPath);
if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
cfg=ipm_c.config();cfg.grid.nr=321;cfg.grid.ntheta=257;
cfg.grid.rmax=1e10;cfg.grid.radialMap='sinh';cfg.grid.stencilWidth=7;
cfg.runtime.threads=1;cfg.output.verbose=false;cfg.time.cfl=.4;
firstSpacings=[.00996755832309634,.02];labels={'baseline','candidate'};
savedThreads=maxNumCompThreads(1);
threadCleanup=onCleanup(@()maxNumCompThreads(savedThreads)); %#ok<NASGU>
report=struct('kind','polar_initial_native_F_Q_nested_radial_angular_screen_v1', ...
    'status','running','completed',false,'passed',false,'threads',1, ...
    'coarseNr',321,'coarseNtheta',257,'fineNr',641,'fineNtheta',513, ...
    'firstSpacings',firstSpacings,'rmax',1e10,'coreMaximumRadius',2, ...
    'radialRelativeNodeTolerance',1e-11,'angularAbsoluteNodeTolerance',eps(pi), ...
    'cases',struct([]),'candidateOverBaselineErrorRatios',struct(), ...
    'timeAdvanced',0,'efficiencyMeasured',false,'accuracyThresholdPreselected',false, ...
    'accuracyQualificationClaimed',false,'referenceConvergenceEstablished',false, ...
    'interpretation',['Two-level full native spatial differences at exactly matching nodes; ', ...
       'all errors retained without a qualification threshold. Analytic source is an independent initial-data oracle.']);
ipm_c.output.atomicJson(outputPath,report);
for group=1:numel(firstSpacings)
    report.status='running';report.activeCase=labels{group};
    try
        cfg.grid.nr=321;cfg.grid.ntheta=257;cfg.grid.firstSpacing=firstSpacings(group);
        timer=tic;state=ipm_c.initialize(cfg);elapsed=toc(timer);
        coarse=fields(state);g=state.ops.polarGrid;alpha=g.radialMapParameter;
        exact=ipm_c_legacy_fixture('background',cfg.physics,g.X,g.Y,state.scale);
        core=g.R<=2;
        row=struct('label',labels{group},'firstSpacing',firstSpacings(group), ...
            'sinhParameter',alpha,'coarse',metadata(state,elapsed), ...
            'fine',struct(),'matchingNodes',struct(),'coarseVersusFine',struct(), ...
            'analyticSourceCoarse',region_errors(coarse.source,exact.rx,core), ...
            'analyticSourceFine',struct(),'rates',struct(),'wallIPC',struct(), ...
            'completed',false);
        report.cases=append(report.cases,row);
        clear state exact g core
        ipm_c.output.atomicJson(outputPath,report);

        cfg.grid.nr=641;cfg.grid.ntheta=513;
        cfg.grid.firstSpacing=cfg.grid.rmax*sinh(alpha/(cfg.grid.nr-1))/sinh(alpha);
        timer=tic;state=ipm_c.initialize(cfg);elapsed=toc(timer);
        fine=fields(state);g=state.ops.polarGrid;
        exact=ipm_c_legacy_fixture('background',cfg.physics,g.X,g.Y,state.scale);
        row.fine=metadata(state,elapsed);
        row.analyticSourceFine=region_errors(fine.source,exact.rx,g.R<=2);
        fineR=fine.r(1:2:end);fineTheta=fine.theta(1:2:end);
        radialMismatch=max(abs(fineR-coarse.r)./max(abs(coarse.r),realmin));
        angularMismatch=max(abs(fineTheta-coarse.theta));
        row.matchingNodes=struct('radialRelativeDifferenceInfinity',radialMismatch, ...
            'angularAbsoluteDifferenceInfinity',angularMismatch, ...
            'radialStride',2,'angularStride',2, ...
            'interpolationUsed',false);
        report.cases(group)=row;ipm_c.output.atomicJson(outputPath,report);
        assert(radialMismatch<=report.radialRelativeNodeTolerance&& ...
            angularMismatch<=report.angularAbsoluteNodeTolerance, ...
            'ipm_c:RadialCheckPhysicalNodes','The coarse and fine physical nodes must match before comparing fields.');
        core=repmat(coarse.r<=2,1,numel(coarse.theta));
        for key={'F','Q','source','u1','u2'}
            name=key{1};fineField=fine.(name);
            row.coarseVersusFine.(name)=region_errors(coarse.(name),fineField(1:2:end,1:2:end),core);
        end
        velocityError=max(abs(coarse.u1-fine.u1(1:2:end,1:2:end)), ...
            abs(coarse.u2-fine.u2(1:2:end,1:2:end)));
        referenceVelocity=max(abs(fine.u1(1:2:end,1:2:end)),abs(fine.u2(1:2:end,1:2:end)));
        row.coarseVersusFine.velocity=region_error_magnitudes(velocityError,referenceVelocity,core);
        row.rates=struct('names',{{'c_l','c_omega','c_r'}}, ...
            'coarse',coarse.rates,'fine',fine.rates, ...
            'absoluteDifference',abs(coarse.rates-fine.rates), ...
            'relativeDifference',abs(coarse.rates-fine.rates)./max(abs(fine.rates),realmin), ...
            'differenceDividedByMaximumOfOneAndReference', ...
            abs(coarse.rates-fine.rates)./max(1,abs(fine.rates)));
        row.wallIPC=struct('names',{{'I','P','C'}}, ...
            'coarse',coarse.IPC,'fine',fine.IPC, ...
            'absoluteDifference',abs(coarse.IPC-fine.IPC));
        row.completed=true;report.cases(group)=row;
        clear state g exact fine coarse fineField velocityError referenceVelocity core
        ipm_c.output.atomicJson(outputPath,report);
    catch problem
        report.status='failed';report.failure=struct('identifier',problem.identifier,'message',problem.message);
        ipm_c.output.atomicJson(outputPath,report);rethrow(problem);
    end
end
baseline=report.cases(1);candidate=report.cases(2);
ratios=struct('coarseVersusFine',struct(),'analyticSourceCoarse',struct(),'analyticSourceFine',struct(), ...
    'rateAbsoluteDifference',candidate.rates.absoluteDifference./max(baseline.rates.absoluteDifference,realmin));
for key={'F','Q','source','u1','u2','velocity'}
    name=key{1};ratios.coarseVersusFine.(name)=region_ratios( ...
        candidate.coarseVersusFine.(name),baseline.coarseVersusFine.(name));
end
ratios.analyticSourceCoarse=region_ratios(candidate.analyticSourceCoarse,baseline.analyticSourceCoarse);
ratios.analyticSourceFine=region_ratios(candidate.analyticSourceFine,baseline.analyticSourceFine);
report.candidateOverBaselineErrorRatios=ratios;
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
value=struct('nr',g.nr,'ntheta',g.ntheta,'firstSpacing',g.r(2), ...
    'originFitRadius',g.originFit.scale,'sinhParameter',g.radialMapParameter, ...
    'setupSeconds',elapsed,'transportCflRate',f.generatorCflRate,'reactionCflRate',abs(f.c_omega), ...
    'dtCfl',state.config.time.cfl/(f.generatorCflRate+abs(f.c_omega)), ...
    'cflMaximumPosition',f.ipmCInfo.cflPosition,'poisson',f.poissonSolveInfo, ...
    'gauge',struct('rootResidual',a.rootResidual,'auditResidual',a.auditResidual, ...
        'rootTolerance',a.rootTolerance,'auditTolerance',a.auditTolerance, ...
        'rootCondition',a.rootCondition,'iterations',a.iterations,'evaluations',a.evaluations));
end

function result=region_errors(value,reference,core)
result=region_error_magnitudes(abs(value-reference),abs(reference),core);
end

function result=region_error_magnitudes(errorMagnitude,referenceMagnitude,core)
first=false(size(core));first(1:min(5,size(core,1)),:)=true;
result=struct('global',one_region(errorMagnitude,referenceMagnitude,true(size(core))), ...
    'coreRadiusAtMostTwo',one_region(errorMagnitude,referenceMagnitude,core), ...
    'originAndFirstFourRings',one_region(errorMagnitude,referenceMagnitude,first));
end

function result=one_region(errorMagnitude,referenceMagnitude,mask)
e=max(errorMagnitude(mask));r=max(referenceMagnitude(mask));
result=struct('absoluteInfinity',e,'referenceInfinity',r,'relativeInfinity',e/max(r,realmin));
end

function result=region_ratios(candidate,baseline)
result=struct();
for key={'global','coreRadiusAtMostTwo','originAndFirstFourRings'}
    name=key{1};c=candidate.(name);b=baseline.(name);
    result.(name)=struct('absoluteErrorRatio',c.absoluteInfinity/max(b.absoluteInfinity,realmin), ...
        'relativeErrorRatio',c.relativeInfinity/max(b.relativeInfinity,realmin));
end
end

function result=append(result,row)
if isempty(result),result=row;else,result(end+1)=row;end
end
