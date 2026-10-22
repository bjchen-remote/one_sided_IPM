function fileName = uniquePath(requestedFile,caseId)
%IPM.OUTPUT.UNIQUEPATH Preserve an existing artifact by choosing a suffix.

requestedFile = char(requestedFile);
if isempty(requestedFile)
    fileName = requestedFile;
    return;
end
if isfolder(requestedFile)
    error('ipm:OutputPathDirectory', ...
        'Configured output path is a directory: %s',requestedFile);
end
if ~isfile(requestedFile)
    fileName = requestedFile;
    return;
end

[directory,stem,extension] = fileparts(requestedFile);
safeCaseId = regexprep(char(caseId),'[^A-Za-z0-9_.-]','_');
baseName = sprintf('%s_%s',stem,safeCaseId);
fileName = located_name(directory,[baseName,extension]);
suffix = 2;
while isfile(fileName) || isfolder(fileName)
    fileName = located_name(directory, ...
        sprintf('%s_%d%s',baseName,suffix,extension));
    suffix = suffix+1;
end
end

function fileName = located_name(directory,name)
if isempty(directory)
    fileName = name;
else
    fileName = fullfile(directory,name);
end
end
