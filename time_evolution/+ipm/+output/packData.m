function data = packData(inputFile,outputFile)
%IPM.OUTPUT.PACKDATA Produce one portable MAT with all native frame payloads.
% Daily runs keep immutable frames separately. Packing is explicit and costs
% memory proportional to the selected archive, never part of an RK step.
assert(~strcmp(char(inputFile),char(outputFile)), ...
    'ipm:DatasetPackPath','Use a different output file from the live catalog.');
saved=load(inputFile,'data'); data=saved.data;
assert(data.schemaVersion==2,'ipm:DatasetSchema','Expected schema-2 data.');
data=ipm.output.withResidualNorm(data);
[ok,attributes]=fileattrib(inputFile); assert(ok);
data.meta.baseFolder=fileparts(attributes.Name);
for k=1:numel(data.frames)
    frame=ipm.output.loadFrame(data,k);
    data.frames(k).payload=frame;
    data.frames(k).file='';
end
data.meta.storage='packed_native_frames';
data.meta.packedFrom=attributes.Name;
data.meta.datasetFile=char(outputFile);
data.meta.baseFolder=fileparts(char(outputFile));
folder=fileparts(char(outputFile)); if isempty(folder), folder='.'; end
if ~isfolder(folder), mkdir(folder); end
assert(~isfile(outputFile),'ipm:DatasetPackExists','Choose a new archive file.');
temporary=[tempname(folder),'.mat'];
guard=onCleanup(@()delete_if_exists(temporary));
save(temporary,'data','-v7.3');
[ok,message]=movefile(temporary,outputFile); assert(ok,'%s',message);
end

function delete_if_exists(file)
if isfile(file), delete(file); end
end
