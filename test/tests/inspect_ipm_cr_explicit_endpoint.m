function report=inspect_ipm_cr_explicit_endpoint(options)
%INSPECT_IPM_CR_EXPLICIT_ENDPOINT Strict Poisson refinement at one accepted rho.
% No time advance, density transfer, gauge-reference reset or changed WENO
% operator is performed. OPTIONS.checkpointFolder is required explicitly;
% expectedTau defaults to 200. The recorded spectral configuration is
% followed by [64,64]/64/16 and [96,96]/128/24, all requireAccuracy=true.
% States carrying eta use evaluateResidual at the saved scale/background;
% other states use evaluate. Reconstructed complete rho must be identical.
% Recorded files, current package files and actual MATLAB resolutions must
% match the checkpoint's SHA identity. A fresh outputDir is optional.
% This diagnoses sensitivity of small Q to the Poisson discretization; it
% does not establish long-time spatial accuracy or an exact velocity.

if nargin<1,options=struct();end
assert(isstruct(options)&&isscalar(options),'ipm_cr:EndpointOptions');
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'test','tests'),fullfile(root,'time_evolution'));
allowed={'checkpointFolder','expectedTau','outputDir'};
assert(all(ismember(fieldnames(options),allowed)),'ipm_cr:EndpointOptions', ...
    'Use only checkpointFolder, expectedTau and outputDir.');
assert(isfield(options,'checkpointFolder')&& ...
    (ischar(options.checkpointFolder)||(isstring(options.checkpointFolder)&&isscalar(options.checkpointFolder))), ...
    'ipm_cr:EndpointFolder','Supply an explicit absolute checkpointFolder.');
folder=options.checkpointFolder;
folder=char(folder);expectedTau=option(options,'expectedTau',200);
validateattributes(expectedTau,{'numeric'},{'scalar','real','finite','positive'});
assert(startsWith(folder,filesep),'ipm_cr:EndpointFolder','checkpointFolder must be absolute.');
[cp,pointer]=load_endpoint(folder,expectedTau);
[sourceBindings,resolvedFunctions]=verify_sources(cp.sourceInventory,root);
if isfield(options,'outputDir')
    outputDir=char(options.outputDir);
else
    parent=fullfile(root,'test','results');[~,token]=fileparts(tempname(parent));
    outputDir=fullfile(parent,['explicit_endpoint_poisson_', ...
        char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')),'_',token]);
end
assert(startsWith(outputDir,filesep)&&~isfolder(outputDir)&&~isfile(outputDir), ...
    'ipm_cr:EndpointFreshOutput','Choose a fresh absolute endpoint-inspection output directory.');
mkdir(outputDir);timer=tic;
state=cp.state;g=cp.geometry;wallOps=cp.wallOps;method=cp.explicitOptions;
hasResidual=isfield(state,'residual');
if hasResidual
    sourceDefinition=['Fixed saved eta and scalar scale, including the compensated physical clock. ', ...
        'The same synchronized native D_X eta plus analytic background G_X is used at every level; ', ...
        'the analytic background is evaluated directly at each spectral grid and includes any saved physical-time tail.'];
    evaluationPath='ipm_cr.explicit.evaluateResidual';
else
    sourceDefinition=['The same saved centered rho is synchronized and differentiated on the same native FD transport grids. ', ...
        'Only interpolation to the different spectral grids changes source reconstruction.'];
    evaluationPath='ipm_cr.explicit.evaluate';
end
report=struct('kind','fixed_state_strict_Poisson_refinement_endpoint_inspection_v2', ...
    'status','running','outputDir',outputDir,'checkpoint',pointer,'expectedTau',expectedTau, ...
    'actualTau',state.scale.canonicalTime,'acceptedStep',state.step, ...
    'recordedSourceHashesVerified',true,'currentPackageSourceHashesVerified',true, ...
    'resolvedPackageSourceHashesVerified',true,'sourceBindings',sourceBindings, ...
    'resolvedEvaluationFunctions',resolvedFunctions,'evaluationPath',evaluationPath, ...
    'residualBackgroundRepresentation',hasResidual,'timeAdvanced',0,'densityTransferred',false, ...
    'transportGeometryChanged',false,'transportOptionsChanged',false, ...
    'wallReferencesReset',false,'requireAccuracy',true,'accuracyGatesRelaxed',false, ...
    'fixedNativeSourceDefinition',sourceDefinition, ...
    'sourceReconstructionExactClaimed',false,'thresholdsSpecified',false, ...
    'smallQSpatiallyCertified',false,'longTimeSpatialAccuracyCertified',false, ...
    'globalAbsoluteVelocityCertified',false,'continuumAccuracyCertified',false, ...
    'regions','All active core nodes; separately each active patch and their combined duplicate degrees of freedom within r<=3.', ...
    'limitations',{{'All levels share one fixed finite-resolution state, its saved background/scale and the same native source definition.', ...
        'Small algebraic, gauge or overlap residuals do not certify the exact Poisson velocity or the exterior tail.', ...
        'Spectral source sampling/reconstruction changes with each Poisson grid and is not asserted exact.', ...
        'A single endpoint Poisson refinement does not test accumulated transport, ghost exchange or time errors.', ...
        'Only the actual final refined level, not an exact continuum solution, is used for pairwise differences.'}});
report.originalWallAuthority=cp.authority;
report.fixedScale=state.scale;
if hasResidual,report.backgroundModel=method.background;end
inspectionPath=[mfilename('fullpath'),'.m'];
report.inspectionSource=struct('path',inspectionPath,'sha256',ipm_c.output.sha256(inspectionPath));
report.transportGeometry=struct('coreSize',size(g.core.X),'polarSize',size(g.polar.X), ...
    'coreSpacing',[g.core.dx,g.core.dy],'activeDofCount',g.activeDofCount, ...
    'outerRadius',g.polar.r(end),'geometryKey',g.geometryKey);
report.recordedEndpoint=flow_summary(state.flow,g,method);
report.recordedEndpointOverlapAudit=ipm_cr.explicit.audit(state,g);
ids={'coarse_fresh','core64_angle64_panel16','core96_angle128_panel24'};
configs=cell(1,3);configs{1}=method.solverOptions;
configs{2}=configs{1};configs{2}.coreDegree=[64,64];
configs{2}.angularDegree=64;configs{2}.polarDegree=16;
configs{3}=configs{1};configs{3}.coreDegree=[96,96];
configs{3}.angularDegree=128;configs{3}.polarDegree=24;
for j=1:3
    configs{j}.rOuter=g.polar.r(end);configs{j}.requireAccuracy=true;
end
report.requestedPoissonLevels=struct('coarse',configs{1},'fine',configs{2},'finest',configs{3});
report.levels=struct([]);report.comparisons=struct();write_report(outputDir,report);
flows=cell(1,3);success=false(1,3);
for j=1:3
    actual=method;actual.solverOptions=configs{j};levelTimer=tic;
    entry=struct('id',ids{j},'solverOptions',configs{j},'status','running', ...
        'evaluateSeconds',NaN,'flowSummary',struct(),'poisson',struct(), ...
        'overlapAudit',struct(),'densityReconstructionIdenticalToRecorded',false, ...
        'densityReconstructionAbsoluteDifference',NaN,'nativeSourceIdenticalToRecorded',false, ...
        'nativeSourceAbsoluteDifferenceToRecorded',NaN,'failure',struct());
    try
        if hasResidual
            [flow,reconstructed]=ipm_cr.explicit.evaluateResidual(state.residual,state.scale, ...
                g,wallOps,actual,[state.flow.c_l,state.flow.c_r]);
        else
            [flow,reconstructed]=ipm_cr.explicit.evaluate(state.rho,g,wallOps,actual, ...
                [state.flow.c_l,state.flow.c_r]);
        end
        entry.evaluateSeconds=toc(levelTimer);
        entry.densityReconstructionIdenticalToRecorded=isequaln(reconstructed,state.rho);
        entry.densityReconstructionAbsoluteDifference=max([ ...
            max(abs(reconstructed.core-state.rho.core),[],'all'), ...
            max(abs(reconstructed.polar-state.rho.polar),[],'all')]);
        assert(entry.densityReconstructionIdenticalToRecorded, ...
            'ipm_cr:EndpointDensityChanged', ...
            'State-aware reconstruction changed saved rho (maximum difference %.17g).', ...
            entry.densityReconstructionAbsoluteDifference);
        entry.flowSummary=flow_summary(flow,g,actual);entry.poisson=flow.poisson;
        entry.nativeSourceIdenticalToRecorded=isequaln(flow.q,state.flow.q);
        entry.nativeSourceAbsoluteDifferenceToRecorded=max([ ...
            max(abs(flow.q.core-state.flow.q.core),[],'all'), ...
            max(abs(flow.q.polar-state.flow.q.polar),[],'all')]);
        auditState=state;auditState.flow=flow;
        entry.overlapAudit=ipm_cr.explicit.audit(auditState,g);
        entry.status='strict_evaluation_completed';flows{j}=flow;success(j)=true;
    catch problem
        entry.status='strict_evaluation_failed';entry.evaluateSeconds=toc(levelTimer);
        entry.failure=struct('identifier',problem.identifier,'message',problem.message, ...
            'stack',problem.stack,'accuracyGateRelaxed',false, ...
            'diagnosticScope','Original exception from the strict full evaluation; no failed solve is rerun with weakened accuracy requirements.');
    end
    if isempty(report.levels),report.levels=entry;else,report.levels(end+1)=entry;end %#ok<AGROW>
    report.elapsedSeconds=toc(timer);write_report(outputDir,report);
end
if success(1)
    report.recordedVsFreshCoarse=compare_flows(state.flow,flows{1},g,method);
end
if success(1)&&success(2)
    report.comparisons.coarseVsFine=compare_flows(flows{1},flows{2},g,method);
end
if success(1)&&success(3)
    report.comparisons.coarseVsFinest=compare_flows(flows{1},flows{3},g,method);
end
if success(2)&&success(3)
    report.comparisons.fineVsFinest=compare_flows(flows{2},flows{3},g,method);
end
report.completedStrictLevels=ids(success);report.failedStrictLevels=ids(~success);
if all(success),report.status='completed_strict_finite_resolution_comparison';
else,report.status='completed_with_strict_level_failures';end
report.elapsedSeconds=toc(timer);
report.acceptedCheckpointFileUnchanged=strcmp(ipm_c.output.sha256(pointer.path),pointer.sha256);
[endBindings,endResolved]=verify_sources(cp.sourceInventory,root);
report.sourceIdentitiesUnchanged=isequaln(sourceBindings,endBindings)&&isequaln(resolvedFunctions,endResolved);
assert(report.sourceIdentitiesUnchanged,'ipm_cr:EndpointSourceChanged', ...
    'Recorded, current or resolved numerical source identities changed during inspection.');
write_report(outputDir,report);
fprintf('CR ENDPOINT tau=%.17g strictLevels=%d/%d elapsed=%.3fs output=%s\n', ...
    state.scale.canonicalTime,nnz(success),numel(success),report.elapsedSeconds,outputDir);
end

function [cp,pointer]=load_endpoint(folder,expectedTau)
pointer=jsondecode(fileread(fullfile(folder,'latest.json')));
assert(isstruct(pointer)&&all(isfield(pointer,{'kind','path','sha256','sequence','tau'}))&& ...
    strcmp(pointer.kind,'ipm_cr_explicit_accepted_pointer_v1')&& ...
    any(strcmp(pointer.path,{fullfile(folder,'accepted_0.mat'),fullfile(folder,'accepted_1.mat')}))&& ...
    strcmp(ipm_c.output.sha256(pointer.path),pointer.sha256), ...
    'ipm_cr:EndpointCheckpoint','The actual accepted pointer or SHA identity is invalid.');
assert(pointer.tau==expectedTau,'ipm_cr:EndpointTime', ...
    'Require accepted tau %.17g; actual pointer is %.17g.',expectedTau,pointer.tau);
loaded=load(pointer.path,'checkpoint');cp=loaded.checkpoint;
assert(strcmp(ipm_c.output.sha256(pointer.path),pointer.sha256)&& ...
    strcmp(cp.kind,'ipm_cr_explicit_accepted_checkpoint_v1')&&cp.acceptedOnly&& ...
    cp.sequence==pointer.sequence&&cp.state.step==pointer.sequence&& ...
    cp.state.scale.canonicalTime==expectedTau&&isequaln(cp.rho,cp.state.rho), ...
    'ipm_cr:EndpointCheckpoint','The accepted endpoint payload changed or disagrees with its pointer.');
end

function [bindings,resolvedFunctions]=verify_sources(inventory,root)
assert(isstruct(inventory)&&~isempty(inventory)&&all(isfield(inventory,{'path','sha256'})), ...
    'ipm_cr:EndpointSourceInventory','The checkpoint must contain a numerical source SHA inventory.');
bindings=struct([]);
for j=1:numel(inventory)
    recorded=char(inventory(j).path);expected=char(inventory(j).sha256);
    assert(isfile(recorded)&&strcmp(ipm_c.output.sha256(recorded),expected), ...
        'ipm_cr:EndpointSourceChanged','Recorded numerical source changed: %s.',recorded);
    relative=package_relative_path(recorded);
    entry=struct('recordedPath',recorded,'sha256',expected,'relativePackagePath',relative, ...
        'currentPath','','currentSha256','','resolvedSymbol','','resolvedPath','','resolvedSha256','');
    if ~isempty(relative)
        current=fullfile(root,strrep(relative,'/',filesep));
        assert(isfile(current)&&strcmp(ipm_c.output.sha256(current),expected), ...
            'ipm_cr:EndpointActiveSourceChanged', ...
            'Current package source differs from the accepted checkpoint: %s.',current);
        entry.currentPath=canonical_path(current);entry.currentSha256=expected;
        symbol=package_symbol(relative);entry.resolvedSymbol=symbol;
        if ~isempty(symbol)
            resolved=which(symbol);
            assert(~isempty(resolved)&&isfile(resolved), ...
                'ipm_cr:EndpointSourceResolution','Cannot resolve recorded package symbol %s.',symbol);
            resolved=canonical_path(resolved);
            assert(strcmp(resolved,entry.currentPath)&&strcmp(ipm_c.output.sha256(resolved),expected), ...
                'ipm_cr:EndpointSourceResolution', ...
                'Actual MATLAB resolution for %s differs from the verified current package: %s.',symbol,resolved);
            entry.resolvedPath=resolved;entry.resolvedSha256=expected;
        end
    end
    if isempty(bindings),bindings=entry;else,bindings(end+1)=entry;end %#ok<AGROW>
end
symbols={'ipm_cr.explicit.evaluate','ipm_cr.explicit.evaluateResidual', ...
    'ipm_cr.explicit.poisson','ipm_cr.Solver'};
resolvedFunctions=struct([]);
for j=1:numel(symbols)
    index=find(strcmp({bindings.resolvedSymbol},symbols{j}));
    assert(numel(index)==1,'ipm_cr:EndpointSourceInventory', ...
        'The checkpoint must bind exactly one actual implementation of %s.',symbols{j});
    b=bindings(index);
    entry=struct('symbol',symbols{j},'path',b.resolvedPath,'sha256',b.resolvedSha256, ...
        'recordedPath',b.recordedPath,'currentPath',b.currentPath);
    if isempty(resolvedFunctions),resolvedFunctions=entry;else,resolvedFunctions(end+1)=entry;end %#ok<AGROW>
end
end

function relative=package_relative_path(path)
path=strrep(char(path),'\','/');relative='';
prefixes={'test/+ipm_cr/','time_evolution/+ipm/','time_evolution/+ipm_c/'};
for j=1:numel(prefixes)
    prefix=prefixes{j};positions=strfind(['/',path],['/',prefix]);
    if isempty(positions),continue;end
    relative=path(positions(end):end);
    assert(~any(strcmp(strsplit(relative,'/'),'..')), ...
        'ipm_cr:EndpointSourceInventory','Invalid package-relative source path.');
    return;
end
end

function symbol=package_symbol(relative)
[folder,name,extension]=fileparts(relative);symbol='';
if ~strcmp(extension,'.m')&&~startsWith(extension,'.mex'),return;end
parts=strsplit(strrep(folder,'\','/'),'/');packages=parts(startsWith(parts,'+'));
if isempty(packages),return;end
packages=cellfun(@(part)part(2:end),packages,'UniformOutput',false);
symbol=strjoin([packages,{name}],'.');
end

function path=canonical_path(path)
file=java.io.File(char(path));path=char(file.getCanonicalPath());
end

function s=flow_summary(f,g,method)
core=g.core.active;nearCore=core&hypot(g.core.X,g.core.Y)<=3;
nearPolar=g.polar.active&hypot(g.polar.X,g.polar.Y)<=3;
s=struct('rates',[f.c_l,f.c_omega,f.c_r],'rateNames',{{'c_l','c_omega','c_r'}}, ...
    'delta',f.delta,'Uref',f.Uref,'transportOrigin',f.transportOrigin, ...
    'transportCflRate',f.transportRate,'totalCflRate',f.rate, ...
    'endpointDtAtConfiguredCfl',method.cfl/max(f.rate,realmin), ...
    'endpointDtIsNotAnAllStageStep',true,'cflPosition',f.cflPosition, ...
    'gaugeResidual',f.gaugeResidual,'wallIPC',[f.measurement.I,f.measurement.P,f.measurement.C], ...
    'QInfinityActiveCore',max(abs(f.Q.core(core))), ...
    'QInfinityActiveNearR3',max([max(abs(f.Q.core(nearCore))),max(abs(f.Q.polar(nearPolar)))]), ...
    'QInfinityAllActive',max([max(abs(f.Q.core(core))),max(abs(f.Q.polar(g.polar.active)))]), ...
    'FInfinityActiveCore',max(abs(f.F.core(core))), ...
    'FInfinityActiveNearR3',max([max(abs(f.F.core(nearCore))),max(abs(f.F.polar(nearPolar)))]));
end

function result=compare_flows(a,b,g,method)
result=struct();near=struct();
near.core=g.core.active&hypot(g.core.X,g.core.Y)<=3;
near.polar=g.polar.active&hypot(g.polar.X,g.polar.Y)<=3;
result.coreActive=region_metrics(a,b,g,struct('core',g.core.active));
result.nearR3.core=region_metrics(a,b,g,struct('core',near.core));
result.nearR3.polar=region_metrics(a,b,g,struct('polar',near.polar));
result.nearR3.combined=region_metrics(a,b,g,near);
names={'c_l','c_omega','c_r','delta','Uref','transportOrigin','transportRate','rate'};
av=zeros(1,numel(names));bv=av;
for j=1:numel(names),av(j)=a.(names{j});bv(j)=b.(names{j});end
result.ratesAndCfl=struct('names',{names},'firstLevel',av,'secondLevel',bv, ...
    'signedDifference',av-bv,'mixedRelativeDifference',abs(av-bv)./max(1,abs(bv)));
result.endpointDtAtConfiguredCfl=[method.cfl/max(a.rate,realmin),method.cfl/max(b.rate,realmin)];
result.wallIPCSignedDifference=[a.measurement.I-b.measurement.I, ...
    a.measurement.P-b.measurement.P,a.measurement.C-b.measurement.C];
result.firstGaugeResidual=a.gaugeResidual;result.secondGaugeResidual=b.gaugeResidual;
result.timeAdvanced=0;result.transportGeometryIdentical=true;
end

function result=region_metrics(a,b,g,masks)
names={'q','F','Q','rawU1','rawU2','UrefCenteredU1'};result=struct();
for j=1:numel(names)
    name=names{j};[av,x,y,patch]=collect(a,g,masks,name);bv=collect(b,g,masks,name);
    result.fields.(name)=metric(av,bv,x,y,patch);
end
au1=collect(a,g,masks,'rawU1');au2=collect(a,g,masks,'rawU2');
bu1=collect(b,g,masks,'rawU1');bu2=collect(b,g,masks,'rawU2');
result.rawVelocity=vector_metric(au1,au2,bu1,bu2,x,y,patch);
result.UrefCenteredVelocity=vector_metric(au1-a.Uref,au2,bu1-b.Uref,bu2,x,y,patch);
result.nativeSamples=numel(av);result.activeOverlapIncludedForEverySelectedPatch=true;
end

function [values,x,y,patch]=collect(flow,g,masks,name)
values=[];x=[];y=[];patch=[];
for item={'core','polar'}
    key=item{1};if ~isfield(masks,key),continue;end
    mask=masks.(key);assert(any(mask,'all'),'ipm_cr:EndpointComparisonRegion');
    switch name
        case {'q','F','Q'},field=flow.(name).(key);
        case 'rawU1',field=flow.physical.(key).u1;
        case 'rawU2',field=flow.physical.(key).u2;
        case 'UrefCenteredU1',field=flow.physical.(key).u1-flow.Uref;
    end
    assert(all(isfinite(field(mask))),'ipm_cr:EndpointComparisonFinite');
    values=[values;field(mask)];x=[x;g.(key).X(mask)];y=[y;g.(key).Y(mask)]; %#ok<AGROW>
    patch=[patch;(1+strcmp(key,'polar'))*ones(nnz(mask),1)]; %#ok<AGROW>
end
end

function m=metric(first,second,x,y,patch)
d=first-second;[maximum,index]=max(abs(d));
m=scale_record(maximum,max(abs(first)),max(abs(second)));
m.rms=sqrt(mean(d.^2));m.signedDifferenceAtMaximum=d(index);
m.firstValueAtMaximum=first(index);m.secondValueAtMaximum=second(index);
m.maximumLocation=location(x,y,patch,index);
end

function m=vector_metric(ax,ay,bx,by,x,y,patch)
d=hypot(ax-bx,ay-by);[maximum,index]=max(d);
m=scale_record(maximum,max(hypot(ax,ay)),max(hypot(bx,by)));
m.rms=sqrt(mean(d.^2));m.maximumLocation=location(x,y,patch,index);
end

function m=scale_record(error,firstScale,secondScale)
m=struct('absoluteInfinity',error,'firstInfinityScale',firstScale, ...
    'secondInfinityScale',secondScale,'normalizedByFirstScale',error/max(firstScale,realmin), ...
    'normalizedBySecondScale',error/max(secondScale,realmin), ...
    'mixedRelativeInfinity',error/max([1,firstScale,secondScale]));
end

function p=location(x,y,patch,index)
name='core';if patch(index)==2,name='polar';end
p=struct('patch',name,'X',x(index),'Y',y(index),'r',hypot(x(index),y(index)));
end

function write_report(folder,report)
temporary=fullfile(folder,'endpoint.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'endpoint.mat'),'f');
assert(ok,'ipm_cr:EndpointReport','%s',message);
temporary=fullfile(folder,'endpoint.json.writing');fid=fopen(temporary,'w');
assert(fid>=0,'ipm_cr:EndpointReport');cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'endpoint.json'),'f');
assert(ok,'ipm_cr:EndpointReport','%s',message);
end

function value=option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
