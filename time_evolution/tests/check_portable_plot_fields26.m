function report=check_portable_plot_fields26(outputRoot,inputCheckpoint)
% One native SDIRK2 step from the same signed CP, output on/off, pure MATLAB.
assert(startsWith(outputRoot,filesep) && ~isfolder(outputRoot));
packageRoot=fileparts(fileparts(mfilename('fullpath')));
assert(startsWith(which('ipm.solve'),[packageRoot,filesep]));
assert(isempty(which('ipm.field.weno5LocalMex')));
manifest=jsondecode(fileread(fullfile(packageRoot,'MANIFEST.json')));
verify_files(packageRoot,manifest.bundleRows);
inputHash=local_sha(inputCheckpoint);cp=ipm.output.readCheckpoint(inputCheckpoint);
requestedCanonicalAdvance=.2;
saved=cp.payload.state;assert(strcmp(saved.config.time.timeIntegrator,'sdirk2_alexander_practical_v2'));
assert(strcmp(saved.config.physics.initialCondition,'localized_gaussian_primitive'), ...
    'Use this new26 Gaussian acceptance checkpoint, not a legacy25 IVP.');
mkdir(outputRoot);timer=tic;
report=struct('kind','portable_plot_fields26_actual_one_step_io_on_off_checkpoint', ...
    'status','running','passed',false,'pureMatlab',true,'goalHoldPassed',false, ...
    'runtimeTreeSha256',manifest.runtimeTreeSha256,'inputCheckpoint',inputCheckpoint, ...
    'inputCheckpointSha256',inputHash,'sourceStep',saved.step, ...
    'sourceTau',saved.scale.canonicalTime, ...
    'requestedCanonicalAdvance',requestedCanonicalAdvance,'sharedDatasetTouched',false);
persist(outputRoot,report);
try
    results=cell(1,2);checkpoints=cell(1,2);states=cell(1,2);
    for k=1:2
        label='off';if k==1,label='on';end
        folder=fullfile(outputRoot,label);mkdir(folder);
        metadata=saved.config.output.caseMetadata;
        metadata.fullRhoXRateTarget=1e-5;metadata.fullRhoXRateHoldTau=.5;
        metadata.fullRhoXRateMonitorFile=fullfile(folder,'monitor.csv');
        overrides=struct('finalTime',saved.scale.canonicalTime+requestedCanonicalAdvance, ...
            'maxSteps',saved.step+1,'verbose',false,'saveResults',false, ...
            'makePlots',false,'livePlot',false,'writeVideo',false,'caseMetadata',metadata, ...
            'checkpoint',struct('enabled',true,'file',fullfile(folder,'accepted.mat'), ...
                'every',Inf,'atExit',true));
        policy=struct('enabled',k==1,'file',fullfile(folder,'data.mat'),'frameEvery',.01);
        results{k}=ipm.solve(overrides,inputCheckpoint,policy);
        checkpoints{k}=last_checkpoint(folder);
        [states{k},~,~]=ipm.output.restoreCheckpoint(checkpoints{k});
        assert(states{k}.step==saved.step+1 && states{k}.scale.canonicalTime>saved.scale.canonicalTime);
        assert(abs(states{k}.scale.canonicalTime- ...
            (saved.scale.canonicalTime+requestedCanonicalAdvance))<=16*eps(1), ...
            'I/O must use the same actual default h=.2 step in both branches.');
        assert(states{k}.ops.rescaling.referenceWallIntegral==saved.rescaling.referenceWallIntegral && ...
            states{k}.ops.rescaling.referenceWallGradientCentroid==saved.rescaling.referenceWallGradientCentroid);
    end
    assert(isequaln(results{1}.state,results{2}.state) && ...
        isequaln(results{1}.grid,results{2}.grid) && isequaln(results{1}.scale,results{2}.scale), ...
        'I/O changed the actual accepted numerical result.');
    assert(isequaln(states{1}.rho,states{2}.rho) && isequaln(states{1}.scale,states{2}.scale));
    assert(isequaln(states{1}.flow.source,states{2}.flow.source) && ...
        isequaln(states{1}.flow.psi,states{2}.flow.psi) && ...
        isequaln(states{1}.flow.u1,states{2}.flow.u1) && isequaln(states{1}.flow.u2,states{2}.flow.u2));
    assert(isequaln(states{1}.rho,results{1}.state.rho) && ...
        isequaln(states{1}.flow.source,results{1}.state.omega));
    assert(~isfile(fullfile(outputRoot,'off','data.mat')));
    dataFile=fullfile(outputRoot,'on','data.mat');data=ipm.output.loadData(dataFile);
    assert(isequaln(data.residual_inf,data.fullRateNorm(:)) && iscolumn(data.residual_inf));
    tau=states{1}.scale.canonicalTime;selected=ipm_select_time(data,tau,'exact');
    assert(isequaln(selected.vx,selected.transport_u1) && ...
        isequaln(selected.vy,selected.transport_u2) && ...
        isequaln(selected.vmag,hypot(selected.vx,selected.vy)));
    full=max(abs(selected.rho_xt),[],'all');j=find(data.tau==tau & data.step==states{1}.step,1,'last');
    assert(~isempty(j) && data.residual_inf(j)==full);
    restoredQ=states{1}.rhsCache.rhoRate*states{1}.ops.Dx';
    qRestoreError=max(abs(restoredQ-selected.rho_xt),[],'all');
    assert(qRestoreError<=128*eps(max(1,full)));
    vxFormula=states{1}.flow.u1+states{1}.flow.c_l*states{1}.ops.X+states{1}.flow.c_r;
    vyFormula=states{1}.flow.u2+states{1}.flow.c_l*states{1}.ops.Y;
    vxError=max(abs(vxFormula-selected.vx),[],'all');vyError=max(abs(vyFormula-selected.vy),[],'all');
    assert(all(abs(vxFormula-selected.vx)<=128*eps(max(1,abs(vxFormula))),'all') && ...
        all(abs(vyFormula-selected.vy)<=128*eps(max(1,abs(vyFormula))),'all'));
    packedFile=fullfile(outputRoot,'packed.mat');ipm.output.packData(dataFile,packedFile);
    packed=ipm.output.loadData(packedFile);packed=ipm_select_time(packed,tau,'exact');
    assert(isequaln(packed.residual_inf,selected.residual_inf) && ...
        isequaln(packed.vx,selected.vx) && isequaln(packed.vy,selected.vy) && ...
        isequaln(packed.vmag,selected.vmag) && isequaln(packed.rho_xt,selected.rho_xt));
    verify_files(packageRoot,manifest.bundleRows);assert(strcmp(local_sha(inputCheckpoint),inputHash));
    report.actualTau=tau;report.actualPhysicalTime=states{1}.scale.physicalTime;
    report.actualStep=states{1}.step;report.actualCanonicalAdvance=tau-saved.scale.canonicalTime;
    report.onOffAcceptedStateBitwiseEqual=true;report.checkpointRestoredStatesBitwiseEqual=true;
    report.checkpoints=checkpoints;report.dataset=dataFile;report.packedDataset=packedFile;
    report.fullQRate=full;report.residualInfAliasExact=true;
    report.compositeVelocityAliasesExact=true;report.vmagDefinitionExact=true;
    report.compositeFormulaErrors=[vxError,vyError];report.fullRhsRestoreRoundoffError=qRestoreError;
    report.originalI0=saved.rescaling.referenceWallIntegral;
    report.originalC0=saved.rescaling.referenceWallGradientCentroid;
    report.originalReferencesRetained=true;report.allCandidateBytesUnchanged=true;
    report.status='one_real_step_io_bitmatch_checkpoint_alias_pack_passed_not_goal_hold';report.passed=true;
catch e
    report.status='failed';report.errorIdentifier=e.identifier;
    report.errorMessage=getReport(e,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputRoot,report);rethrow(e)
end
report.elapsedSeconds=toc(timer);persist(outputRoot,report);
fprintf('PORTABLE_PLOT_FIELDS26 passed=%d elapsed=%.3fs\n',report.passed,report.elapsedSeconds);
end
function p=last_checkpoint(folder)
lines=splitlines(string(fileread(fullfile(folder,'checkpoint_manifest.jsonl'))));
lines=lines(strlength(strtrim(lines))>0);row=jsondecode(lines(end));p=row.checkpointFile;
end
function verify_files(root,rows)
for k=1:numel(rows),assert(strcmp(local_sha(fullfile(root,rows(k).path)),rows(k).sha256));end
end
function value=local_sha(path)
fid=fopen(path,'rb');assert(fid>0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
function persist(folder,report)
save(fullfile(folder,'report.mat'),'report','-v7.3');fid=fopen(fullfile(folder,'report.json'),'w');assert(fid>0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
end
