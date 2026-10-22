function report=inspect_ipm_cr_explicit_initial(outputDir,options)
%INSPECT_IPM_CR_EXPLICIT_INITIAL Read-only spatial/CFL initial-state screen.
% The original analytic datum and immutable wall authority are shared by all
% levels. No accepted time step, remesh or new gauge reference is generated.
% The high-resolution oracle uses an independently evaluated analytic source
% and continuum advective RHS with the ORIGINAL discrete wall functionals.
% It is a finite-resolution comparator, not an exact velocity or trajectory.
if nargin<2,options=struct();end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'time_evolution'));
if nargin<1||isempty(outputDir)
    outputDir=fullfile(root,'test','results','explicit_initial_inspection');
end
assert(isstruct(options)&&isscalar(options),'ipm_cr:InitialInspectionOptions', ...
    'Use scalar initial-inspection options.');
includeOracle=option(options,'includeOracle',true);
includeCombined=option(options,'includeCombined',true);
includeConfigDefaults=option(options,'includeConfigDefaults',false);
pinSources=option(options,'pinSources',true);
for flag={includeOracle,includeCombined,includeConfigDefaults,pinSources}
    validateattributes(flag{1},{'logical'},{'scalar'});
end
if ~isfolder(outputDir),mkdir(outputDir);end
path=fullfile(outputDir,'initial_inspection.json');
assert(~isfile(path),'ipm_cr:InitialInspectionPreserve','Use a fresh numerical inspection directory.');
started=tic;cfg=ipm_c.config();physics=cfg.physics;
fineGeometry=struct('coreSpacing',.015, ...
    'grid',struct('nr',641,'ntheta',256,'firstSpacing',.015));
finePoisson=struct('coreDegree',[64,64],'angularDegree',64,'polarDegree',16,'rOuter',1e10);
plan=struct('name',{'default','finer_transport','finer_poisson'}, ...
    'overrides',{struct(),struct('geometry',fineGeometry),struct('solverOptions',finePoisson)});
if includeCombined
    plan(end+1)=struct('name','finer_transport_and_poisson', ...
        'overrides',struct('geometry',fineGeometry,'solverOptions',finePoisson));
end
report=struct('kind','explicit_hybrid_initial_spatial_and_physical_CFL_inspection_v1', ...
    'status','initializing','tau',0,'timeAdvanced',false,'remeshApplied',false, ...
    'originalAuthorityPreserved',true,'wallReferencesReset',false, ...
    'physics',physics,'authorityFile',cfg.authorityFile, ...
    'authoritySha256',ipm_c.output.sha256(cfg.authorityFile),'comparisonCfl',.3, ...
    'explicitRecords',struct([]),'nativeRecords',struct([]),'reference',struct('available',false), ...
    'sourcePins',struct([]),'sourcePinsUnchanged',false,'elapsedSeconds',0, ...
    'analyticSourceReferenceRequested',includeOracle, ...
    'accuracyCertified',false,'equalAccuracySpeedupClaimed',false, ...
    'qualification','Initial spatial and time-step diagnostics only; discretizations and errors differ, so dt ratios are not equal-accuracy speedups.');
if pinSources,report.sourcePins=source_pins();end
persist(path,report);firstWall=[];firstAuthority=[];reference=[];
try
    for k=1:numel(plan)
        timer=tic;
        [state,g,wallOps,resolved]=ipm_cr.explicit.initialize(plan(k).overrides);
        prepareSeconds=toc(timer);
        if isempty(firstWall)
            firstWall=wallOps;firstAuthority=state.originalAuthority;
            report.originalReferences=state.originalAuthority;
            report.originalWallNodes=numel(wallOps.x);
            if includeOracle
                try
                    [reference,report.reference]=make_reference(physics,wallOps,options);
                catch referenceProblem
                    % Preserve the independently useful source/CFL screen if
                    % the optional high-resolution oracle cannot qualify.
                    report.reference=struct('available',false, ...
                        'errorIdentifier',referenceProblem.identifier, ...
                        'errorMessage',getReport(referenceProblem,'extended','hyperlinks','off'));
                    includeOracle=false;
                end
            end
        else
            assert(isequaln(firstWall,wallOps)&&isequaln(firstAuthority,state.originalAuthority), ...
                'ipm_cr:InitialInspectionAuthority','A refinement replaced the original wall operators or references.');
        end
        qErrors=source_errors(state.flow.q,g,physics);
        row=struct('name',plan(k).name,'overrides',plan(k).overrides, ...
            'resolvedOptions',resolved,'initializationSeconds',prepareSeconds, ...
            'geometry',geometry_summary(g),'sourceErrors',qErrors, ...
            'flow',explicit_summary(state.flow,g),'poisson',state.flow.poisson, ...
            'initialTimeStep',dt_summary(state.flow.rate,resolved.maxDt), ...
            'oracleComparison',struct('available',false));
        if includeOracle
            row.oracleComparison=compare_explicit(state,g,reference,physics);
        end
        report.explicitRecords=append(report.explicitRecords,row);
        report.elapsedSeconds=toc(started);persist(path,report);
        fprintf('CR INITIAL %-29s rate %.6g dt(.3) %.6g core q %.3e annular q %.3e\n', ...
            row.name,state.flow.rate,row.initialTimeStep.nominalCflLimitedDt, ...
            qErrors.core.activeNearInfinity,qErrors.polar.activeNearInfinity);
        clear state g
    end

    % The saved tau200 experiment uses 385-by-128 transport. Vary only its
    % first radial spacing to expose the actual origin/angular CFL penalty.
    nativePlan=struct('name',{'saved_baseline_385x128_h06','saved_baseline_385x128_h03'}, ...
        'nr',{385,385},'ntheta',{128,128},'spacing',{.06,.03});
    if includeConfigDefaults
        nativePlan(end+1)=struct('name','config_default_nodes_h06', ...
            'nr',cfg.grid.nr,'ntheta',cfg.grid.ntheta,'spacing',.06);
        nativePlan(end+1)=struct('name','config_default_nodes_h03', ...
            'nr',cfg.grid.nr,'ntheta',cfg.grid.ntheta,'spacing',.03);
    end
    for k=1:numel(nativePlan)
        a=nativePlan(k);nativeCfg=cfg;nativeCfg.grid.nr=a.nr;nativeCfg.grid.ntheta=a.ntheta;
        nativeCfg.grid.firstSpacing=a.spacing;nativeCfg.remesh.minimumFirstSpacing=a.spacing;
        nativeCfg.time.cfl=.3;nativeCfg.runtime.threads=1;nativeCfg.output.verbose=false;
        timer=tic;[state,context]=ipm_c.initialize(nativeCfg);seconds=toc(timer);
        assert(isequaln(context.cart.originalWallOps,firstWall), ...
            'ipm_cr:InitialInspectionAuthority','The native comparator changed the original wall authority.');
        flow=state.flow;rate=flow.generatorCflRate+abs(flow.c_omega);
        [~,position]=ipm_c.spatial.transportRate(flow.transportU1,flow.transportU2,context.g);
        q=flow.source;F=state.rhsCache.rhoRate;
        row=struct('name',a.name,'gridOptions',nativeCfg.grid,'initializationSeconds',seconds, ...
            'nr',context.g.nr,'ntheta',context.g.ntheta,'independentDofCount',1+(context.g.nr-1)*context.g.ntheta, ...
            'firstPositiveRadius',context.g.r(2),'fullTransportRate',flow.generatorCflRate, ...
            'reactionRate',abs(flow.c_omega),'fullRate',rate,'cflPosition',position, ...
            'rates',[flow.c_l,flow.c_omega,flow.c_r],'gaugeCenteredRates', ...
            [flow.c_l,flow.c_omega,flow.generatorTranslationRate], ...
            'Uref',flow.generatorReferenceVelocity,'gaugeResidual',flow.gaugeResidual, ...
            'initialTimeStep',dt_summary(rate,nativeCfg.time.maxDt), ...
            'sourceNearInfinity',max(abs(q(context.g.r<=3,:)- ...
                source(context.g.X(context.g.r<=3,:),context.g.Y(context.g.r<=3,:),physics)),[],'all'), ...
            'oracleComparison',struct('available',false));
        if includeOracle
            row.oracleComparison=compare_native(state,context.g,reference,physics);
        end
        row.hybridDefaultOverNativeNominalDt=report.explicitRecords(1).initialTimeStep.nominalCflLimitedDt/ ...
            row.initialTimeStep.nominalCflLimitedDt;
        row.timeStepRatioIsEqualAccuracySpeedup=false;
        report.nativeRecords=append(report.nativeRecords,row);
        report.elapsedSeconds=toc(started);persist(path,report);
        fprintf('NATIVE INITIAL %-29s rate %.6g dt(.3) %.6g limiter r %.6g theta %.6g\n', ...
            a.name,rate,row.initialTimeStep.nominalCflLimitedDt,position.r,position.theta);
        clear state context q F
    end
    if pinSources
        report.sourcePinsUnchanged=isequaln(report.sourcePins,source_pins());
        assert(report.sourcePinsUnchanged,'ipm_cr:InitialInspectionSources','A pinned source changed during inspection.');
    end
    report.status='completed_initial_spatial_and_CFL_screen';
    if report.analyticSourceReferenceRequested&&~includeOracle
        report.status='completed_source_CFL_screen_analytic_reference_unavailable';
    end
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(started);persist(path,report);rethrow(problem);
end
report.elapsedSeconds=toc(started);persist(path,report);
end

function [ref,summary]=make_reference(p,wallOps,options)
cfg=struct('coreDegree',[96,96],'angularDegree',128,'polarDegree',24, ...
    'rOuter',1e10,'requireAccuracy',true);
if isfield(options,'referenceSolverOptions')
    names=fieldnames(options.referenceSolverOptions);
    for k=1:numel(names),cfg.(names{k})=options.referenceSolverOptions.(names{k});end
end
timer=tic;solver=ipm_cr.Solver(cfg);buildSeconds=toc(timer);
timer=tic;[solution,info]=solver.solve(@(x,y)source(x,y,p));solveSeconds=toc(timer);
x=wallOps.x(:)';wall=density(x,zeros(size(x)),p);a=density(0,0,p);W=wall-a;
q=ipm_c.gauge.wallDerivative(W,wallOps.Dx);
integral=ipm.evolve.wallGradientIntegral(q,wallOps,wallOps.rescaling);
peak=ipm.evolve.globalSplineMaximum(q,wallOps);
centroid=ipm.evolve.wallGradientCentroid(q,wallOps);
m=struct('I',integral.value,'P',peak.value,'C',centroid.value, ...
    'integral',integral,'peak',peak,'centroid',centroid);
wallVelocity=solver.evaluate(solution,x,zeros(size(x)));
atReference=solver.evaluate(solution,m.C,0);Uref=atReference.u1;
qx=source(x,zeros(size(x)),p);
B=-(wallVelocity.u1-Uref).*qx;L=-x.*qx;T=-qx;
observe=@(f)rate_observations(f,m,wallOps);
matrix=[observe(L),observe(T),observe(W)];
target=[0;0;wallOps.rescaling.peakPositionRestoringRate*(m.C-wallOps.rescaling.referenceWallGradientCentroid)];
rates=matrix\(-observe(B)-target);
actual=B+rates(1)*L+rates(2)*T+rates(3)*W;
residual=observe(actual)+target;
assert(all(isfinite(rates))&&max(abs(residual))<=1e-8, ...
    'ipm_cr:InitialReferenceGauge','The independent analytic-source gauge failed.');
atOrigin=solver.evaluate(solution,0,0);
F0=-((atOrigin.u1-Uref)+rates(2))*source(0,0,p);
ref=struct('solver',solver,'solution',solution,'Uref',Uref, ...
    'c_l',rates(1),'delta',rates(2),'c_omega',rates(3),'c_r',rates(2)-Uref, ...
    'originAmplitude',a,'F0',F0,'measurement',m);
summary=struct('available',true,'solverOptions',cfg,'buildSeconds',buildSeconds, ...
    'solveSeconds',solveSeconds,'poisson',info,'rates',[ref.c_l,ref.c_omega,ref.c_r], ...
    'gaugeCenteredRates',[ref.c_l,ref.c_omega,ref.delta],'Uref',Uref, ...
    'originalWallMeasurement',[m.I,m.P,m.C],'gaugeResidual',residual, ...
    'gaugeMatrixCondition',cond(matrix),'transportOrigin',F0, ...
    'sourceDefinition','Independent analytic derivative of the original Gaussian primitive; no native source interpolation.', ...
    'rhsDefinition','Continuum advective density RHS evaluated with the high-resolution Poisson field, then observed by the original discrete wall I/P/C functionals.', ...
    'exactVelocityKnown',false,'referenceRefinementCertified',false, ...
    'sameSemiDiscreteWenoOperator',false,'continuumAccuracyCertified',false);
end

function v=rate_observations(field,m,wallOps)
dq=ipm_c.gauge.wallDerivative(field,wallOps.Dx);
v=[m.integral.weights*dq'/m.I; ...
    m.peak.weights*dq(m.peak.indices)'/m.P; ...
    m.centroid.weights*dq(m.centroid.indices)'];
end

function comparison=compare_explicit(state,g,ref,p)
f=state.flow;comparison=struct('available',true);
for name={'core','polar'}
    key=name{1};patch=g.(key);selected=patch.active&hypot(patch.X,patch.Y)<=3;
    comparison.(key)=compare_samples(patch.X(selected),patch.Y(selected), ...
        f.q.(key)(selected),f.F.(key)(selected),f.physical.(key).u1(selected), ...
        f.physical.(key).u2(selected),f.Uref,ref,p);
end
comparison.rates=[f.c_l,f.c_omega,f.c_r]-[ref.c_l,ref.c_omega,ref.c_r];
comparison.gaugeCenteredRates=[f.c_l,f.c_omega,f.delta]-[ref.c_l,ref.c_omega,ref.delta];
comparison.UrefDifference=f.Uref-ref.Uref;
comparison.description='Actual WENO RHS/source and own-Uref-centered velocity versus independent analytic-source/continuum-advection reference at identical active native nodes, r<=3.';
end

function comparison=compare_native(state,g,ref,p)
f=state.flow;comparison=struct('available',true);
near=hypot(g.X,g.Y)<=3;core=abs(g.X)<=1.11&g.Y<=1.11;
for k=1:2
    if k==1,selected=core;name='core';else,selected=near;name='near';end
    comparison.(name)=compare_samples(g.X(selected),g.Y(selected), ...
        f.source(selected),state.rhsCache.rhoRate(selected),f.u1(selected),f.u2(selected), ...
        f.generatorReferenceVelocity,ref,p);
end
comparison.rates=[f.c_l,f.c_omega,f.c_r]-[ref.c_l,ref.c_omega,ref.c_r];
comparison.gaugeCenteredRates=[f.c_l,f.c_omega,f.generatorTranslationRate]-[ref.c_l,ref.c_omega,ref.delta];
comparison.UrefDifference=f.generatorReferenceVelocity-ref.Uref;
comparison.description='Original polar full RHS/source and own-Uref-centered velocity versus the same analytic-source reference; original native discretization is not exact truth.';
end

function error=compare_samples(x,y,q,F,u1,u2,Uref,ref,p)
count=numel(x);exactU1=zeros(count,1);exactU2=exactU1;
for first=1:2048:count
    indices=first:min(first+2047,count);
    v=ref.solver.evaluate(ref.solution,x(indices),y(indices));
    exactU1(indices)=v.u1;exactU2(indices)=v.u2;
end
qx=source(x,y,p);R=density(x,y,p);qy=-2*y/p.gaussianWidthY^2.*R;
W=R-ref.originAmplitude;
Fref=-((exactU1-ref.Uref)+ref.c_l*x+ref.delta).*qx- ...
    (exactU2+ref.c_l*y).*qy+ref.c_omega*W-ref.F0;
du=[abs((u1-Uref)-(exactU1-ref.Uref));abs(u2-exactU2)];
error=struct('samples',count,'maximumRadius',max(hypot(x,y)), ...
    'sourceQAbsoluteInfinity',max(abs(q-qx)), ...
    'UrefCenteredVelocityAbsoluteInfinity',max(du), ...
    'centeredDensityRhsAbsoluteInfinity',max(abs(F-Fref)), ...
    'referenceVelocityInfinity',max([abs(exactU1-ref.Uref);abs(exactU2)]), ...
    'referenceRhsInfinity',max(abs(Fref)));
end

function e=source_errors(q,g,p)
e=struct();
for name={'core','polar'}
    key=name{1};patch=g.(key);exact=source(patch.X,patch.Y,p);
    near=patch.active&hypot(patch.X,patch.Y)<=3;difference=abs(q.(key)-exact);
    e.(key)=struct('activeNearInfinity',max(difference(near)), ...
        'activeCompleteInfinity',max(difference(patch.active)), ...
        'completeNativeInfinity',max(difference,[],'all'),'activeNearSamples',nnz(near));
end
end

function value=density(x,y,p)
z=(x-p.gaussianCenterX)/p.gaussianWidthX;w=p.gaussianSkewWeight;
primitive=(1-w)*erfc(-z)+w*erfc(-(z-.8)/1.25);
value=p.gaussianAmplitude*p.gaussianWidthX*sqrt(pi)/2*primitive.*exp(-(y/p.gaussianWidthY).^2);
end

function q=source(x,y,p)
z=(x-p.gaussianCenterX)/p.gaussianWidthX;w=p.gaussianSkewWeight;
q=p.gaussianAmplitude*((1-w)*exp(-z.^2)+w/1.25*exp(-((z-.8)/1.25).^2)).*exp(-(y/p.gaussianWidthY).^2);
end

function s=dt_summary(rate,maxDt)
s=struct('fullTransportPlusReactionRate',rate,'comparisonCfl',.3, ...
    'unclippedCflDt',.3/max(rate,realmin),'maxDt',maxDt, ...
    'nominalCflLimitedDt',min(maxDt,.3/max(rate,realmin)), ...
    'withPoint95SafetyMargin',min(maxDt,.95*.3/max(rate,realmin)), ...
    'actualAcceptedStep',false,'equalAccuracySpeedupClaimed',false);
end

function s=explicit_summary(f,g)
s=struct('rates',[f.c_l,f.c_omega,f.c_r],'gaugeCenteredRates',[f.c_l,f.c_omega,f.delta], ...
    'Uref',f.Uref,'transportOrigin',f.transportOrigin,'transportRate',f.transportRate, ...
    'reactionRate',abs(f.c_omega),'fullRate',f.rate,'cflPosition',f.cflPosition, ...
    'gaugeResidual',f.gaugeResidual,'I',f.measurement.I,'P',f.measurement.P,'C',f.measurement.C, ...
    'coreActiveQRateInfinity',max(abs(f.Q.core(g.core.active))), ...
    'polarActiveQRateInfinity',max(abs(f.Q.polar(g.polar.active))), ...
    'coreCompleteNativeQRateInfinity',max(abs(f.Q.core),[],'all'), ...
    'polarCompleteNativeQRateInfinity',max(abs(f.Q.polar),[],'all'));
end

function s=geometry_summary(g)
s=struct('coreNx',g.core.nx,'coreNy',g.core.ny,'coreDx',g.core.dx,'coreDy',g.core.dy, ...
    'polarNr',g.polar.nr,'polarNtheta',g.polar.ntheta,'polarInnerRadius',g.polar.r(1), ...
    'polarInnerActiveRadius',g.polar.innerActiveRadius,'outerRadius',g.polar.r(end), ...
    'activeDofCount',g.activeDofCount,'originIsCartesian',true,'polarOriginPresent',false);
end

function rows=source_pins()
names={'ipm_cr.explicit.initialize','ipm_cr.explicit.geometry','ipm_cr.explicit.synchronize', ...
    'ipm_cr.explicit.differentiate','ipm_cr.explicit.transport','ipm_cr.explicit.poisson', ...
    'ipm_cr.explicit.interpolate','ipm_cr.explicit.evaluate','ipm_cr.Solver', ...
    'ipm_c.config','ipm_c.initialize','ipm_c.spatial.evaluate','ipm_c.spatial.poisson', ...
    'ipm_c.gauge.initialAuthority','ipm_c.gauge.solve','ipm_c.gauge.measure'};
p=mfilename('fullpath');if ~isfile(p),p=[p,'.m'];end
rows=struct('path',p,'sha256',ipm_c.output.sha256(p));
for k=1:numel(names)
    p=which(names{k});rows=append(rows,struct('path',p,'sha256',ipm_c.output.sha256(p)));
end
end

function persist(path,report)
temporary=[path,'.writing'];fid=fopen(temporary,'w');
assert(fid>=0,'ipm_cr:InitialInspectionReport','Cannot write the initial numerical screen.');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,path,'f');assert(ok,'ipm_cr:InitialInspectionReport','%s',message);
end

function s=append(s,row)
if isempty(s),s=row;else,s(end+1)=row;end
end

function value=option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
