function report=run_mesh_cfl_screen(outputPath)
%RUN_MESH_CFL_SCREEN Initial-state geometry/CFL screen, without time advance.
% Each first-ring radius has its own 513-angle reference on IDENTICAL radial
% nodes.  Its subsamples compare complete native F/Q/source/velocity without
% interpolation. This isolates angular resolution; radial accuracy is not
% established by this screen. The exact Gaussian source is a test oracle.
assert(ischar(outputPath)||isstring(outputPath),'ipm_c:MeshScreenOutput');
outputPath=char(outputPath);parent=fileparts(outputPath);
if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
cfg=ipm_c.config();baseFirstSpacing=.00996755832309634;
cfg.grid.nr=321;cfg.grid.rmax=1e10;cfg.grid.radialMap='sinh';
cfg.grid.stencilWidth=7;cfg.runtime.threads=1;cfg.output.verbose=false;
angles=[257,129,65];factors=[1,2,4];referenceAngles=513;
savedThreads=maxNumCompThreads(1);
threadCleanup=onCleanup(@()maxNumCompThreads(savedThreads)); %#ok<NASGU>
report=struct('kind','polar_initial_geometry_CFL_and_angular_reference_screen_v1', ...
    'status','running','completed',false,'passed',false, ...
    'nr',321,'rmax',1e10,'baseFirstSpacing',baseFirstSpacing, ...
    'angleCounts',angles,'firstSpacingFactors',factors, ...
    'referenceAngleCount',referenceAngles,'coreMaximumRadius',2, ...
    'threads',1,'references',struct([]),'cases',struct([]), ...
    'timeAdvanced',0,'efficiencyMeasured',false, ...
    'accuracyThresholdPreselected',false,'accuracyQualificationClaimed',false, ...
    'radialReferenceRefinementPerformed',false, ...
    'interpretation',['Angular subsampling comparisons and an exact initial Gaussian source oracle only; ', ...
      'native F/Q radial refinement and equal-final-time efficiency remain separate work.']);
ipm_c.output.atomicJson(outputPath,report);baseline=[];
for group=1:numel(factors)
    cfg.grid.firstSpacing=baseFirstSpacing*factors(group);
    cfg.grid.ntheta=referenceAngles;
    try
        timer=tic;state=ipm_c.initialize(cfg);referenceSetup=toc(timer);
        reference=fields(state);
        g=state.ops.polarGrid;
        exact=ipm_c_legacy_fixture('background',cfg.physics,g.X,g.Y,state.scale);
        core=g.R<=2;
        referenceRecord=metadata(state,referenceSetup,factors(group));
        referenceRecord.analyticSource=region_errors(reference.source,exact.rx,core);
        report.references=append(report.references,referenceRecord);
        clear state g exact core
        ipm_c.output.atomicJson(outputPath,report);
    catch problem
        report.status='reference_failed';report.failure=error_record(problem);
        ipm_c.output.atomicJson(outputPath,report);rethrow(problem);
    end
    for angle=1:numel(angles)
        cfg.grid.ntheta=angles(angle);
        try
            timer=tic;state=ipm_c.initialize(cfg);elapsed=toc(timer);
            g=state.ops.polarGrid;value=fields(state);
            stride=(referenceAngles-1)/(g.ntheta-1);
            assert(stride==fix(stride)&&stride>=1,'ipm_c:MeshScreenNestedAngles');
            index=1:stride:referenceAngles;
            assert(isequal(reference.r,g.r)&& ...
                max(abs(reference.theta(index)-g.theta))<=8*eps(pi), ...
                'ipm_c:MeshScreenSamePhysicalNodes');
            core=g.R<=2;exact=ipm_c_legacy_fixture('background',cfg.physics,g.X,g.Y,state.scale);
            row=metadata(state,elapsed,factors(group));
            row.referenceAngleCount=referenceAngles;row.referenceStride=stride;
            row.analyticSource=region_errors(value.source,exact.rx,core);
            row.angularReference=struct();
            for key={'F','Q','source','u1','u2'}
                name=key{1};fineField=reference.(name);fine=fineField(:,index);
                row.angularReference.(name)=region_errors(value.(name),fine,core);
            end
            velocityError=max(abs(value.u1-reference.u1(:,index)), ...
                abs(value.u2-reference.u2(:,index)));
            referenceVelocity=max(abs(reference.u1(:,index)),abs(reference.u2(:,index)));
            row.angularReference.velocity=region_error_magnitudes(velocityError,referenceVelocity,core);
            row.rates=struct('names',{{'c_l','c_omega','c_r'}}, ...
                'values',value.rates,'reference',reference.rates, ...
                'absoluteDifference',abs(value.rates-reference.rates), ...
                'relativeDifference',abs(value.rates-reference.rates)./max(abs(reference.rates),realmin), ...
                'differenceDividedByMaximumOfOneAndReference', ...
                abs(value.rates-reference.rates)./max(1,abs(reference.rates)));
            row.wallIPC=struct('names',{{'I','P','C'}},'values',value.IPC, ...
                'reference',reference.IPC,'absoluteDifference',abs(value.IPC-reference.IPC));
            if isempty(baseline),baseline=row;end
            row.relativeToBaseline=baseline_comparison(row,baseline);
            row.success=true;row.errorIdentifier='';row.errorMessage='';
            report.cases=append(report.cases,row);
            clear state g value exact fine fineField velocityError referenceVelocity core
            ipm_c.output.atomicJson(outputPath,report);
        catch problem
            report.status='case_failed';report.failure=error_record(problem);
            report.failedCase=struct('firstSpacingFactor',factors(group),'ntheta',angles(angle));
            ipm_c.output.atomicJson(outputPath,report);rethrow(problem);
        end
    end
    clear reference
end
report.status='completed';report.completed=true;report.passed=true;
ipm_c.output.atomicJson(outputPath,report);
end

function v=fields(state)
g=state.ops.polarGrid;flow=state.flow;
v=struct('r',g.r,'theta',g.theta,'F',state.rhsCache.rhoRate, ...
    'Q',flow.fullQ,'source',flow.source,'u1',flow.u1,'u2',flow.u2, ...
    'rates',[flow.c_l,flow.c_omega,flow.c_r], ...
    'IPC',[flow.wallIntegral,flow.globalWallPeak,flow.wallCentroid]);
end

function row=metadata(state,elapsed,factor)
g=state.ops.polarGrid;flow=state.flow;
transport=flow.generatorCflRate;reaction=abs(flow.c_omega);
row=struct('nr',g.nr,'ntheta',g.ntheta,'firstSpacingFactor',factor, ...
    'firstSpacing',g.r(2),'originFitRadius',g.originFit.scale, ...
    'independentDofs',1+(g.nr-1)*g.ntheta, ...
    'setupSeconds',elapsed,'transportCflRate',transport,'reactionCflRate',reaction, ...
    'totalCflRate',transport+reaction,'cfl',state.config.time.cfl, ...
    'dtCfl',state.config.time.cfl/(transport+reaction), ...
    'cflMaximumPosition',flow.ipmCInfo.cflPosition, ...
    'poissonResidual',flow.poissonResidual, ...
    'gaugeResidual',max(abs(flow.gaugeResidual)), ...
    'factorCacheRebuilt',flow.poissonSolveInfo.factorCacheRebuilt);
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
result=struct('absoluteInfinity',e,'referenceInfinity',r, ...
    'relativeInfinity',e/max(r,realmin));
end

function comparison=baseline_comparison(row,baseline)
comparison=struct('baselineFirstSpacing',baseline.firstSpacing, ...
    'baselineAngleCount',baseline.ntheta,'cflStepGain',row.dtCfl/baseline.dtCfl, ...
    'setupSecondsRatio',row.setupSeconds/baseline.setupSeconds, ...
    'setupIsNotPerStepThroughput',true, ...
    'analyticSourceGlobalErrorRatio',ratio(row.analyticSource.global.absoluteInfinity, ...
        baseline.analyticSource.global.absoluteInfinity), ...
    'analyticSourceCoreErrorRatio',ratio(row.analyticSource.coreRadiusAtMostTwo.absoluteInfinity, ...
        baseline.analyticSource.coreRadiusAtMostTwo.absoluteInfinity), ...
    'angularReferenceErrorRatios',struct());
for key={'F','Q','source','u1','u2','velocity'}
    name=key{1};current=row.angularReference.(name);old=baseline.angularReference.(name);
    comparison.angularReferenceErrorRatios.(name)=struct( ...
        'global',ratio(current.global.absoluteInfinity,old.global.absoluteInfinity), ...
        'coreRadiusAtMostTwo',ratio(current.coreRadiusAtMostTwo.absoluteInfinity,old.coreRadiusAtMostTwo.absoluteInfinity), ...
        'originAndFirstFourRings',ratio(current.originAndFirstFourRings.absoluteInfinity,old.originAndFirstFourRings.absoluteInfinity));
end
end

function result=ratio(numerator,denominator)
result=numerator/max(denominator,realmin);
end

function result=error_record(problem)
result=struct('identifier',problem.identifier,'message',problem.message);
end

function result=append(result,row)
if isempty(result),result=row;else,result(end+1)=row;end
end
