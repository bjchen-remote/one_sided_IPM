function [checkpoint,fileName] = writeCheckpoint(checkpoint,requestedFile)
%IPM.OUTPUT.WRITECHECKPOINT Atomically install an immutable checkpoint.

checkpoint = ipm.output.readCheckpoint(checkpoint);
if ~(ischar(requestedFile) || ...
        (isstring(requestedFile) && isscalar(requestedFile)))
    error('ipm:CheckpointPath','Checkpoint output path must be text.');
end
requestedFile = char(requestedFile);
[directory,stem,extension] = fileparts(requestedFile);
if isempty(stem)
    error('ipm:CheckpointPath','Checkpoint output path must name a file.');
end
if isempty(extension)
    extension = '.mat';
end
if isempty(directory)
    directory = '.';
elseif ~isfolder(directory)
    mkdir(directory);
end
caseId = checkpoint.payload.state.runMetadata.caseId;
safeCaseId = regexprep(char(caseId),'[^A-Za-z0-9_.-]','_');
step = checkpoint.payload.state.step;
candidate = fullfile(directory,sprintf('%s_%s_step%010d%s', ...
    stem,safeCaseId,step,extension));
fileName = ipm.output.uniquePath(candidate,caseId);
checkpoint.storageFile = string(fileName);

temporaryFile = [tempname(directory),'.mat'];
temporaryCleanup = onCleanup(@()delete_temporary_file(temporaryFile));
save(temporaryFile,'checkpoint','-v7.3');
[moved,message] = movefile(temporaryFile,fileName);
if ~moved
    error('ipm:CheckpointMove', ...
        'Could not install checkpoint %s: %s',fileName,message);
end
append_manifest(checkpoint,fileName,directory);
end

function append_manifest(checkpoint,fileName,directory)
state = checkpoint.payload.state;
record = struct('caseId',state.runMetadata.caseId, ...
    'createdUtc',checkpoint.createdUtc, ...
    'schemaVersion',checkpoint.schemaVersion, ...
    'checkpointFile',fileName,'step',state.step, ...
    'normalizedTime',state.normalizedTime, ...
    'canonicalTime',state.scale.canonicalTime, ...
    'physicalTime',state.scale.physicalTime, ...
    'grid',[numel(state.x),numel(state.y)]);
manifestFile = fullfile(directory,'checkpoint_manifest.jsonl');
fileId = fopen(manifestFile,'a');
if fileId < 0
    warning('ipm:CheckpointManifestWrite', ...
        'Saved the checkpoint but could not append %s.',manifestFile);
    return;
end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,'%s\n',jsonencode(record));
end

function delete_temporary_file(fileName)
if isfile(fileName)
    try
        delete(fileName);
    catch exception
        fprintf(2,'IPM warning: could not remove checkpoint temp %s: %s\n', ...
            fileName,exception.message);
    end
end
end
