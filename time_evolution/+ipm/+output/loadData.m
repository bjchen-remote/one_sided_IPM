function data = loadData(file)
%IPM.OUTPUT.LOADDATA Load the small catalog and its latest native field.
%   data=ipm.output.loadData('/run/data.mat'); plot(data.tau,data.cr)
%   data=ipm.output.selectTime(data,3.5); mesh(data.x,data.y,data.rho_xt)
assert(isfile(file),'ipm:DatasetMissing','Dataset does not exist: %s',file);
saved=load(file,'data'); data=saved.data;
assert(data.schemaVersion==2,'ipm:DatasetSchema','Expected schema-2 data.');
data=ipm.output.withResidualNorm(data);
[ok,attributes]=fileattrib(file); assert(ok);
data.meta.datasetFile=attributes.Name;
data.meta.baseFolder=fileparts(attributes.Name);
if ~isempty(data.frames)
    data=ipm.output.selectTime(data,numel(data.frames),'index');
end
end
