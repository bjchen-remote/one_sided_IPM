function result = resume_same_equation(checkpoint,outputRoot,finalTau)
%RESUME_SAME_EQUATION Continue a main.m checkpoint in a NEW directory.
% Only terminal/output controls are changed. The initial condition, box,
% mesh, C/gauge choices, transport, and PDE remain frozen in the checkpoint.
% Use main('resume',checkpoint,newFolder,finalTau).
arguments
    checkpoint (1,:) char
    outputRoot (1,:) char
    finalTau (1,1) double {mustBePositive}
end
assert(isfile(checkpoint) && startsWith(checkpoint,filesep));
assert(startsWith(outputRoot,filesep) && ...
    ~isfolder(outputRoot) && ~isfile(outputRoot), ...
    'Choose a NEW absolute output directory.');
workflowRoot=fileparts(mfilename('fullpath'));
solverRoot=fileparts(fileparts(workflowRoot));
addpath(solverRoot,'-begin');
assert(startsWith(which('ipm.solve'),[solverRoot,filesep]));

% Verify the signed accepted state before making a new output directory.
threads=maxNumCompThreads;
maxNumCompThreads(1);
try
    cp=ipm.output.readCheckpoint(checkpoint);
catch errorInfo
    maxNumCompThreads(threads);
    rethrow(errorInfo);
end
maxNumCompThreads(threads);
state=cp.payload.state;
startTau=state.scale.canonicalTime;
assert(finalTau>startTau,'finalTau must exceed the accepted start tau.');
monitorFile=fullfile(outputRoot,'monitor_full_absolute_rate.csv');
resultFile=fullfile(outputRoot,'result_or_target.mat');
checkpointFile=fullfile(outputRoot,'checkpoints','checkpoint.mat');
metadata=state.config.output.caseMetadata;
assert(isfield(metadata,'fullRhoXRateTarget') && ...
    isfield(metadata,'fullRhoXRateHoldTau') && ...
    isfield(metadata,'fullRhoXRateMonitorFile'), ...
    'Continue a checkpoint created through main.m.');
metadata.continuationExperiment='main_same_equation_continuation';
metadata.continuationSource=checkpoint;
metadata.previousMonitorFile=metadata.fullRhoXRateMonitorFile;
metadata.fullRhoXRateMonitorFile=monitorFile;
metadata.targetHoldPreservedFromCheckpoint=true;
overrides=struct('finalTime',finalTau, ...
    'maxSteps',max(state.config.time.maxSteps,state.step+1000000), ...
    'saveResults',true,'verbose',true, ...
    'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'resultFile',resultFile,'caseMetadata',metadata, ...
    'checkpoint',struct('enabled',true,'file',checkpointFile, ...
        'every',0.5,'atExit',true));
mkdir(outputRoot);
save(fullfile(outputRoot,'continuation_parameters.mat'), ...
    'checkpoint','overrides','-v7.3');
maxNumCompThreads(2);
if isfield(metadata,'sourceCflRel')
    setenv('IPM_SOURCE_CFL_REL',char(metadata.sourceCflRel));
else
    % Historical main.m checkpoints used this experiment before its
    % setting became explicit checkpoint metadata. Preserve that trajectory.
    setenv('IPM_SOURCE_CFL_REL','1e-10');
end
diary(fullfile(outputRoot,'run.log'));
diaryGuard=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('output=%s\n',outputRoot);
result=ipm.solve(overrides,checkpoint);
end
