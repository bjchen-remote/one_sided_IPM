function result = write(result,output)
%IPM.OUTPUT.WRITE Persist one versioned result when requested.

if nargin < 2 || isempty(output)
    output = result.config.output;
end
if ~output.saveResults
    result.metadata.resultFile = "";
    return;
end
requestedFile = output.resultFile;
if strlength(result.metadata.resultFile) > 0
    requestedFile = char(result.metadata.resultFile);
end
actualFile = ipm.output.uniquePath(requestedFile,result.metadata.caseId);
result.metadata.resultFile = string(actualFile);
resultDirectory = fileparts(actualFile);
if isempty(resultDirectory)
    resultDirectory = '.';
elseif ~isfolder(resultDirectory)
    mkdir(resultDirectory);
end

% Build the MAT-file beside its destination, then install it only after
% serialization succeeds.  A failed or interrupted save cannot leave a
% partially written configured result file.
temporaryFile = [tempname(resultDirectory),'.mat'];
temporaryCleanup = onCleanup(@()delete_temporary_file(temporaryFile));
save(temporaryFile,'result','-v7.3');
[moved,message] = movefile(temporaryFile,actualFile);
if ~moved
    error('ipm:ResultMove', ...
        'Could not install saved result %s: %s',actualFile,message);
end

manifestFile = fullfile(resultDirectory,'manifest.jsonl');
record = struct('caseId',result.metadata.caseId, ...
    'createdAt',result.metadata.createdAt, ...
    'schemaVersion',result.schemaVersion, ...
    'resultFile',actualFile, ...
    'geometry',result.scale.geometry, ...
    'grid',[numel(result.grid.x),numel(result.grid.y)], ...
    'steps',result.state.steps,'stopReason',result.state.stopReason);
fileId = fopen(manifestFile,'a');
if fileId < 0
    warning('ipm:ManifestWrite', ...
        'Saved the result but could not append %s.',manifestFile);
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
        fprintf(2,'IPM warning: could not remove temporary result %s: %s\n', ...
            fileName,exception.message);
    end
end
end
