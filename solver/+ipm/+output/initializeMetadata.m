function metadata = initializeMetadata(config,ops)
%IPM.OUTPUT.INITIALIZEMETADATA Create one identity and reserve output names.

scaling = config.scaling;
output = config.output;

timestamp = datetime('now','TimeZone','local');
createdTimestamp = timestamp;
createdTimestamp.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX';
caseTimestamp = timestamp;
caseTimestamp.Format = 'yyyyMMdd''T''HHmmssSSS';
[~,token] = fileparts(tempname);
caseId = sprintf('%s_%s_%s_%dx%d_%s',char(caseTimestamp), ...
    scaling.rescalingMode,ops.dynamicScaleGeometry,ops.nx,ops.ny,token);

resultFile = '';
if output.saveResults
    resultFile = ipm.output.uniquePath(output.resultFile,caseId);
end
videoFile = '';
if output.writeVideo
    videoFile = ipm.output.uniquePath(output.videoFile,caseId);
end
metadata = struct('caseId',caseId, ...
    'createdAt',string(createdTimestamp),'solver','IPM', ...
    'caseMetadata',output.caseMetadata, ...
    'gaugeContract',ipm.evolve.gaugeContract(config,ops), ...
    'resultFile',string(resultFile),'videoFile',string(videoFile));
end
