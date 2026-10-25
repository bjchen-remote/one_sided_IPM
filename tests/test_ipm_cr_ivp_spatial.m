function report=test_ipm_cr_ivp_spatial(outputDir,options)
%TEST_IPM_CR_IVP_SPATIAL Three independent original-datum IPM trajectories.
% Reconstruct primary W/q/Q from the analytic background and native eta.
% Full-array reinterpolation is separately retained as a display diagnostic.
% No exterior query is deleted and no spatial order is inferred from a
% localized manufacture. The same global cloud is used on all three grids.
if nargin<2,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','ivp_spatial');end
o=resolve(options);outputDir=char(outputDir);
assert(startsWith(outputDir,filesep),'ipm_cr:IVPOutput');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'ivp_spatial.json')),'ipm_cr:IVPPreserve','Choose a fresh output directory.');
timer=tic;plan=struct('h',{.06,.03,.015},'nr',{193,385,769},'ntheta',{65,129,257});
factor=2^o.refinementStartLevel;
for j=1:3
    plan(j).h=plan(j).h/factor;
    plan(j).nr=(plan(j).nr-1)*factor+1;
    plan(j).ntheta=(plan(j).ntheta-1)*factor+1;
end
geometryKeys=cell(1,3);hotspots=struct([]);
for j=1:3
    g=ipm_cr.mesh.build(geometry_options(plan(j),o.radialRedistribution));
    geometryKeys{j}=g.geometryKey;
    hotspots=append(hotspots,bare_nullspace_hotspot(g,plan(j)));
end
clear g
queries=query_cloud(hotspots);fixedQueries=~isempty(o.queryCloudFile);
if fixedQueries
    reference=load(o.queryCloudFile,'report');queries=reference.report.queryCloud;
    assert(all(isfield(queries,{'x','y','near','wall','hotspot','interfaceBand','remoteWallStrip'})), ...
        'ipm_cr:IVPQueryCloud','Reference must contain the complete unmasked comparison cloud.');
end
rows=struct([]);samples=cell(1,3);
report=struct('kind','original_datum_three_grid_short_IPM_spatial_self_convergence_v1', ...
    'status','running','options',o,'levels',plan,'records',rows, ...
    'queryCloud',queries,'bareNullspaceHotspots',hotspots,'comparisons',struct([]), ...
    'orders',struct([]),'rateOrders',struct([]),'scalarOrders',struct([]), ...
    'timeHalfStep',struct(),'globalFourthOrderObserved',false, ...
    'timeContributionScreenPassed',false,'fullIVPShortWindowQualified',false, ...
    'tau200Qualified',false,'sourceInventory',source_inventory(root), ...
    'sourcesUnchanged',false,'elapsedSeconds',0, ...
    'samplingPolicy','W=analytic H+interp eta; q=analytic Gx+interp synchronized Dx eta; Q=analytic Htx+interp synchronized Dx eta_rhs. Velocities are sampled from actual native Poisson returns; Uref remains separate.', ...
    'rawReinterpolationRole','Direct local interpolation of full native W/q/Q is retained separately as a display/reconstruction diagnostic.', ...
    'limitations',{{'The global query cloud spans the entire domain and resolved/underresolved wall strips, but is a finite sampling set, not a rigorous global infinity bound.', ...
        'No tail mask, post-hoc shift, scale alignment, source clipping or Richardson extrapolation is used.', ...
        'A short original-datum trajectory is distinct from a full 200-tau spatial convergence certificate.', ...
        'Spatial differences, order failures and time contamination are retained even when they fail fourth order.'}});
persist(outputDir,report);authority=[];authorityWall=[];finestInitial=[];finestMethod=[];finestGrid=[];finestWall=[];finestContext=[];
try
    for j=1:3
        method=struct('geometry',geometry_options(plan(j),o.radialRedistribution),'timeIntegrator','rk4','cfl',o.cfl, ...
            'maxDt',o.maxDt,'solverOptions',o.solverOptions);
        [initial,ctx]=ipm_cr.init(method);g=ctx.mesh;wallOps=ctx.wall;method=ctx.opts;
        assert(isequaln(g.geometryKey,geometryKeys{j}),'ipm_cr:IVPGeometry');
        assert(isfield(initial,'eta'),'ipm_cr:IVPRepresentation','Use the analytic-background production path.');
        if isempty(authority)
            authority=initial.originalAuthority;authorityWall=wallOps;
            if ~fixedQueries,queries=append_authority_wall(queries,wallOps.x);end
            report.queryCloud=queries;
        else
            assert(isequaln(authority,initial.originalAuthority)&&isequaln(authorityWall,wallOps), ...
                'ipm_cr:IVPAuthority','All levels must retain the original authority and original wall operators.');
        end
        [state,run]=evolve(initial,g,wallOps,method,o.tauEnd, ...
            fullfile(outputDir,sprintf('level_%d',j)),o.checkpointed,max(0,.9*o.wallSecondsTotal-toc(timer)),ctx);
        if ~run.completed
            report.status='partial_wall_budget_before_three_levels';report.partialLevel=j;report.partialRun=run;
            report.sourcesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
            report.elapsedSeconds=toc(timer);persist(outputDir,report);return;
        end
        [samples{j},diagnostic]=ipm_cr.diag.compare(state,g,method,queries);
        row=struct('level',plan(j),'run',run,'scalars',samples{j}.scalars, ...
            'rates',samples{j}.rates,'samplingDiagnostic',diagnostic, ...
            'actualGeometry',struct('coreSize',[g.core.ny,g.core.nx], ...
                'polarSize',[g.polar.nr,g.polar.ntheta],'innerHaloRadius',g.polar.innerHaloRadius, ...
                'innerActiveRadius',g.polar.innerActiveRadius,'activeDofCount',g.activeDofCount), ...
            'nativeOverlapAudit',ipm_cr.diag.overlap(state,g), ...
            'originalAuthority',initial.originalAuthority);
        rows=append(rows,row);report.records=rows;
        report.sampledStates=samples;report.elapsedSeconds=toc(timer);persist(outputDir,report);
        if j==3,finestInitial=initial;finestMethod=method;finestGrid=g;finestWall=wallOps;finestContext=ctx;end
    end
    comparisons=struct([]);
    for j=1:2
        comparisons=append(comparisons,compare_samples(samples{j},samples{j+1},queries));
    end
    report.comparisons=comparisons;
    names={'W','q','F','Q','u1Centered','u2'};orders=struct([]);passes=true;
    for name=names
        key=name{1};coarse=comparisons(1).fields.(key).global.absoluteInfinity;
        fine=comparisons(2).fields.(key).global.absoluteInfinity;
        order=log2(coarse/fine);floorValue=1024*eps(max(1,max(abs(samples{3}.(key)))));
        aboveFloor=coarse>floorValue&&fine>floorValue;
        pass=aboveFloor&&order>=4;passes=passes&&pass;
        orders=append(orders,struct('quantity',key,'globalDifferences',[coarse,fine], ...
            'globalObservedOrder',order,'nearObservedOrder',log2( ...
            comparisons(1).fields.(key).near.absoluteInfinity/comparisons(2).fields.(key).near.absoluteInfinity), ...
            'roundoffScreen',floorValue,'bothDifferencesAboveRoundoff',aboveFloor, ...
            'nearRoundoff',max(coarse,fine)<=1e-10*max(1,max(abs(samples{3}.(key)))), ...
            'relativeDiagnosticThreshold',1e-10, ...
            'atLeastFourthObserved',pass));
    end
    report.orders=orders;report.globalFourthOrderObserved=passes;
    report.rateOrders=scalar_orders(comparisons(1).rateSignedDifferences, ...
        comparisons(2).rateSignedDifferences,samples{3}.rates);
    report.scalarOrders=scalar_orders(comparisons(1).scalarSignedDifferences, ...
        comparisons(2).scalarSignedDifferences,samples{3}.scalars);
    report.rateAndScalarFourthOrderObserved=all([report.rateOrders.atLeastFourthObserved])&& ...
        all([report.scalarOrders.atLeastFourthObserved]);
    report.status='spatial_comparisons_complete_time_check_pending';
    report.elapsedSeconds=toc(timer);persist(outputDir,report);
    if o.timeHalfCheck
        if toc(timer)>=.9*o.wallSecondsTotal
            report.status='partial_wall_budget_before_time_half';
            report.timeHalfStep=struct('performed',false,'completed',false,'reason','ninety_percent_total_wall_budget', ...
                'timeErrorNegligibleClaimed',false);
            report.sourcesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
            report.elapsedSeconds=toc(timer);persist(outputDir,report);return;
        end
        halfMethod=finestMethod;
        if strcmp(o.timeHalfMode,'cfl')
            refinementFactor=o.timeRefinementFactor;halfMethod.cfl=o.cfl*refinementFactor;
        else,refinementFactor=.5;halfMethod.maxDt=o.maxDt/2;end
        [half,halfRun]=evolve(finestInitial,finestGrid,finestWall,halfMethod,o.tauEnd, ...
            fullfile(outputDir,'half'),o.checkpointed,max(0,.9*o.wallSecondsTotal-toc(timer)),finestContext);
        if ~halfRun.completed
            report.status='partial_wall_budget_during_time_half';
            report.timeHalfStep=struct('performed',true,'completed',false,'mode',o.timeHalfMode, ...
                'refinementFactor',refinementFactor,'run',halfRun, ...
                'timeErrorNegligibleClaimed',false);
            report.sourcesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
            report.elapsedSeconds=toc(timer);persist(outputDir,report);return;
        end
        [halfSample,~]=ipm_cr.diag.compare(half,finestGrid,halfMethod,queries);
        comparison=compare_samples(samples{3},halfSample,queries);ratios=struct();timePassed=true;
        for name=names
            key=name{1};ratio=comparison.fields.(key).global.absoluteInfinity/ ...
                max(comparisons(2).fields.(key).global.absoluteInfinity,realmin);
            ratios.(key)=ratio;timePassed=timePassed&&ratio<=.1;
        end
        ratios.rates=abs(comparison.rateSignedDifferences)./ ...
            max(abs(comparisons(2).rateSignedDifferences),realmin);
        ratios.scalars=abs(comparison.scalarSignedDifferences)./ ...
            max(abs(comparisons(2).scalarSignedDifferences),realmin);
        estimated=ratios;keys=fieldnames(estimated);
        for k=1:numel(keys),estimated.(keys{k})=estimated.(keys{k})/(1-refinementFactor^4);end
        gateRatios=ratios;if strcmp(o.timeHalfMode,'cfl'),gateRatios=estimated;end
        timePassed=true;
        for k=1:numel(keys),timePassed=timePassed&&all(gateRatios.(keys{k})<=.1,'all');end
        report.timeHalfStep=struct('performed',true,'completed',true,'mode',o.timeHalfMode,'run',halfRun,'comparison',comparison, ...
            'refinementFactor',refinementFactor,'timeDifferenceOverFinestSpatialDifference',ratios, ...
            'fourthOrderEstimatedTimeErrorOverFinestSpatialDifference',estimated, ...
            'fourthOrderEstimateAssumesNominalStepRatio',true,'gateUsesFourthOrderEstimate',strcmp(o.timeHalfMode,'cfl'), ...
            'thresholdFraction',.1,'screenPassed',timePassed);
        report.timeContributionScreenPassed=timePassed;
    else,report.timeHalfStep=struct('performed',false,'timeErrorNegligibleClaimed',false);end
    report.fullIVPShortWindowQualified=report.globalFourthOrderObserved&& ...
        report.rateAndScalarFourthOrderObserved&&report.timeContributionScreenPassed;
    report.sourcesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
    assert(report.sourcesUnchanged,'ipm_cr:IVPSources','Numerical source inventory changed.');
    if report.fullIVPShortWindowQualified,report.status='completed_short_IVP_fourth_order_screens';
    else,report.status='completed_short_IVP_measurements_not_qualified';end
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function o=resolve(options)
o=struct('tauEnd',.3,'maxDt',.001,'cfl',.3,'timeHalfCheck',true,'timeHalfMode','max_dt','timeRefinementFactor',.5, ...
    'checkpointed',false,'wallSecondsTotal',Inf,'refinementStartLevel',0, ...
    'radialRedistribution',struct(),'queryCloudFile','', ...
    'solverOptions',struct('coreDegree',[64,64],'angularDegree',80,'polarDegree',20, ...
        'rOuter',1e10,'requireAccuracy',true));
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(o))),'ipm_cr:IVPOptions');
names=fieldnames(options);for j=1:numel(names),o.(names{j})=options.(names{j});end
validateattributes(o.tauEnd,{'numeric'},{'scalar','finite','positive'});
validateattributes(o.maxDt,{'numeric'},{'scalar','finite','positive','<=',.2});
validateattributes(o.cfl,{'numeric'},{'scalar','finite','positive','<=',.3});
validateattributes(o.timeHalfCheck,{'logical'},{'scalar'});
validateattributes(o.checkpointed,{'logical'},{'scalar'});
validateattributes(o.wallSecondsTotal,{'numeric'},{'scalar','real','positive'});
assert(isfinite(o.wallSecondsTotal)||o.wallSecondsTotal==Inf,'ipm_cr:IVPWallBudget');
o.timeHalfMode=char(o.timeHalfMode);assert(ismember(o.timeHalfMode,{'max_dt','cfl'}),'ipm_cr:IVPHalfMode');
validateattributes(o.timeRefinementFactor,{'numeric'},{'scalar','real','finite','positive','<',1});
validateattributes(o.refinementStartLevel,{'numeric'},{'scalar','finite','integer','nonnegative'});
assert(isstruct(o.radialRedistribution)&&isscalar(o.radialRedistribution),'ipm_cr:IVPRadialRedistribution');
o.queryCloudFile=char(o.queryCloudFile);
assert(isempty(o.queryCloudFile)||(startsWith(o.queryCloudFile,filesep)&&isfile(o.queryCloudFile)), ...
    'ipm_cr:IVPQueryCloud','Use an absolute existing MAT reference report.');
end

function [state,run]=evolve(initial,g,wall,method,T,outputFolder,checkpointed,wallBudget,ctx)
timer=tic;state=initial;maximumCfl=0;maximumGauge=0;count=0;minimumDt=Inf;maximumDt=0;retries=0;
if checkpointed&&wallBudget>0
    finiteBudget=wallBudget;if isinf(finiteBudget),finiteBudget=realmax;end
    [summary,state,actualGrid,actualWall]=run_ipm_cr_explicit(struct('outputDir',outputFolder, ...
        'tauEnd',T,'wallSeconds',finiteBudget,'explicitOptions',method,'resume',true, ...
        'checkpointEveryTau',min(1,T),'checkpointWallSeconds',60,'logEveryTau',min(1,T)));
    assert(isequaln(actualGrid.geometryKey,g.geometryKey)&&isequaln(actualWall,wall)&& ...
        isequaln(state.originalAuthority,initial.originalAuthority),'ipm_cr:IVPCheckpointAuthority');
    assert((summary.steps==0&&~isfinite(summary.maximumAcceptedStageCfl))|| ...
        summary.maximumAcceptedStageCfl<=method.cfl*(1+128*eps),'ipm_cr:IVPCfl');
    data=load(fullfile(outputFolder,'summary.mat'),'history','historyColumns');
    dt=data.history(:,strcmp(data.historyColumns,'dt'));dt=dt(isfinite(dt)&dt>0);
    if ~isempty(dt),minimumDt=min(dt);maximumDt=max(dt);end
    run=struct('tau',state.scale.canonicalTime,'steps',summary.steps,'seconds',summary.invocationWallSeconds, ...
        'maximumStageCfl',summary.maximumAcceptedStageCfl,'maximumStageGaugeResidual',summary.maximumRecordedGaugeResidual, ...
        'minimumAcceptedDt',minimumDt,'maximumAcceptedDt',maximumDt,'configuredMaxDt',method.maxDt, ...
        'configuredCfl',method.cfl,'retries',summary.totalStageCflRetries,'fullExplicitCflRequired',true, ...
        'checkpointed',true,'outputDir',outputFolder,'runnerStatus',summary.status, ...
        'acceptedCheckpoint',summary.acceptedCheckpoint,'completed',state.scale.canonicalTime==T);
    return;
end
while state.scale.canonicalTime<T
    if toc(timer)>=wallBudget,break;end
    [state,info]=ipm_cr.time.step(state,g,wall,method,T,ctx);
    assert(strcmp(info.timeIntegrator,'rk4'),'ipm_cr:IVPMethod');
    maximumCfl=max(maximumCfl,info.realizedMaximumStageCfl);
    maximumGauge=max(maximumGauge,info.maximumStageGaugeResidual);count=count+1;
    minimumDt=min(minimumDt,info.dt);maximumDt=max(maximumDt,info.dt);retries=retries+info.retries;
    if mod(count,100)==0,fprintf('IVP SPACE h=%.5g tau=%.6g/%g steps=%d CFL=%.6g\n',g.core.dx,state.scale.canonicalTime,T,count,maximumCfl);end
end
assert(maximumCfl<=method.cfl*(1+128*eps),'ipm_cr:IVPCfl');
run=struct('tau',state.scale.canonicalTime,'steps',count,'seconds',toc(timer), ...
    'maximumStageCfl',maximumCfl,'maximumStageGaugeResidual',maximumGauge, ...
    'minimumAcceptedDt',minimumDt,'maximumAcceptedDt',maximumDt, ...
    'configuredMaxDt',method.maxDt,'configuredCfl',method.cfl,'retries',retries,'fullExplicitCflRequired',true, ...
    'checkpointed',false,'completed',state.scale.canonicalTime==T);
end

function comparison=compare_samples(a,b,q)
comparison=struct('fields',struct());
for name={'W','q','Q','F','u1Centered','u2'}
    key=name{1};comparison.fields.(key)=region_metrics(a.(key)-b.(key),q);
end

comparison.rateNames=a.rates.names;
comparison.rateSignedDifferences=a.rates.values-b.rates.values;
comparison.scalarNames=a.scalars.names;
comparison.scalarSignedDifferences=a.scalars.values-b.scalars.values;
comparison.scalarSignedDifferences(3)=(a.scalars.clockHigh-b.scalars.clockHigh)+ ...
    (a.scalars.clockLow-b.scalars.clockLow);
end

function rows=scalar_orders(coarse,fine,reference)
rows=struct([]);
for j=1:numel(reference.names)
    order=log2(abs(coarse(j))/abs(fine(j)));floorValue=1024*eps(max(1,abs(reference.values(j))));
    aboveFloor=abs(coarse(j))>floorValue&&abs(fine(j))>floorValue;
    rows=append(rows,struct('quantity',reference.names{j},'signedDifferences',[coarse(j),fine(j)], ...
        'observedOrder',order,'roundoffScreen',floorValue,'bothDifferencesAboveRoundoff',aboveFloor, ...
        'nearRoundoff',max(abs([coarse(j),fine(j)]))<=1e-10*max(1,abs(reference.values(j))), ...
        'relativeDiagnosticThreshold',1e-10,'atLeastFourthObserved',aboveFloor&&order>=4));
end
end

function info=region_metrics(d,q)
assert(all(isfinite(d)),'ipm_cr:IVPQueries');
info=struct('global',metric(d,true(size(d)),q),'near',metric(d,q.near,q), ...
    'interfaceBand',metric(d,q.interfaceBand,q),'wall',metric(d,q.wall,q), ...
    'remoteWallStrip',metric(d,q.remoteWallStrip,q),'bareNullspaceHotspots',metric(d,q.hotspot,q));
end

function info=metric(d,mask,q)
indices=find(mask);[value,k]=max(abs(d(indices)));j=indices(k);
info=struct('sampleCount',numel(indices),'absoluteInfinity',value,'rms',sqrt(mean(d(indices).^2)), ...
    'maximumLocation',struct('X',q.x(j),'Y',q.y(j),'r',hypot(q.x(j),q.y(j))));
end

function q=query_cloud(hotspots)
rn=3*((0:48)'+.371)/49;theta=pi*((0:96)+.413)/97;[R,A]=ndgrid(rn,theta);
x=R(:).*cos(A(:));y=R(:).*sin(A(:));
rg=exp(log(.001)+((0:96)'+.217)/97*log(1e13));[R,A]=ndgrid(rg,theta);
x=[x;R(:).*cos(A(:))];y=[y;R(:).*sin(A(:))];
wallX=exp(log(.001)+((0:127)'+.317)/128*log(.99e13));
[X,Y]=ndgrid([-flipud(wallX);0;wallX],[0,.1,.5,1,2]);
x=[x;X(:);-1e10;1e10;0];y=[y;Y(:);0;0;1e10];
hot=false(size(x));
for j=1:numel(hotspots)
    p=hotspots(j);xx=[p.X;-p.X;.9*p.X;1.1*p.X];yy=[p.Y;p.Y;p.Y;p.Y];
    x=[x;xx];y=[y;yy];hot=[hot;true(size(xx))]; %#ok<AGROW>
end
r=hypot(x,y);
q=struct('x',x,'y',y,'near',r<=3,'wall',y==0,'hotspot',hot, ...
    'interfaceBand',abs(r-.81)<=.12| ...
        (abs(abs(x)-1.11)<=.12&y<=1.23)|(abs(y-1.11)<=.12&abs(x)<=1.23), ...
    'remoteWallStrip',r>3&y<=2,'fullRadiusRange',[min(r),max(r)], ...
    'allQueriesRetained',true,'exteriorMaskApplied',false,'globalInfinityBoundClaimed',false);
end

function info=bare_nullspace_hotspot(g,level)
f=struct('core',exp(-g.core.Y.^2),'polar',exp(-g.polar.Y.^2));d=ipm_cr.mesh.diff(f,g);
indices=find(g.polar.active);[error,k]=max(abs(d.polar.x(indices)));i=indices(k);
info=struct('level',level,'X',g.polar.X(i),'Y',g.polar.Y(i),'r',hypot(g.polar.X(i),g.polar.Y(i)), ...
    'rawSourceInfinity',error,'role','Query placement only; this old bare-operator defect is not treated as the new source path.');
end

function q=append_authority_wall(q,wallX)
x=wallX(:);assert(all(abs(x)<=1e10),'ipm_cr:IVPAuthorityWall');
q.x=[q.x;x];q.y=[q.y;zeros(size(x))];q.near=[q.near;abs(x)<=3];
q.wall=[q.wall;true(size(x))];q.hotspot=[q.hotspot;false(size(x))];
q.interfaceBand=[q.interfaceBand;abs(abs(x)-.81)<=.12|abs(abs(x)-1.11)<=.12];
q.remoteWallStrip=[q.remoteWallStrip;abs(x)>3];
end

function o=geometry_options(level,redistribution)
o=struct('coreSpacing',level.h,'grid',struct('nr',level.nr,'ntheta',level.ntheta,'firstSpacing',level.h));
if ~isempty(fieldnames(redistribution)),o.radialRedistribution=redistribution;end
end

function inventory=source_inventory(root)
files=dir(fullfile(root,'src','+ipm_cr','**','*.m'));paths={mfilename('fullpath')};
if ~endsWith(paths{1},'.m'),paths{1}=[paths{1},'.m'];end
for j=1:numel(files),paths{end+1}=fullfile(files(j).folder,files(j).name);end %#ok<AGROW>
paths{end+1}=fullfile(root,'src','+ipm_cr','+op','wenoMex.cpp');
native=which('ipm_cr.op.wenoMex');if ~isempty(native),paths{end+1}=native;end
paths=sort(unique(paths));inventory=struct([]);
for j=1:numel(paths)
    f=dir(paths{j});inventory=append(inventory,struct('path',paths{j},'bytes',f.bytes,'sha256',ipm_cr.io.sha(paths{j})));
end
end

function persist(folder,report)
temporary=fullfile(folder,'ivp_spatial.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'ivp_spatial.mat'),'f');assert(ok,'ipm_cr:IVPSave','%s',message);
temporary=fullfile(folder,'ivp_spatial.json.writing');fid=fopen(temporary,'w');assert(fid>=0,'ipm_cr:IVPSave');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'ivp_spatial.json'),'f');assert(ok,'ipm_cr:IVPSave','%s',message);
end

function a=append(a,b)
if isempty(a),a=b;else,a(end+1)=b;end
end
