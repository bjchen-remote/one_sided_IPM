function inventory=test_inventory(paths)
%TEST_INVENTORY Record explicit test dependencies; no numerical work.
% Preserve caller order. The helper records itself so direct test calls are
% also protected when the shared reporting/enumeration code changes.
self=[mfilename('fullpath'),'.m'];paths=[reshape(paths,1,[]),{self}];
paths=unique(paths,'stable');inventory=struct('path',{},'bytes',{},'sha256',{});
for k=1:numel(paths)
    path=char(paths{k});assert(~isempty(path)&&isfile(path),'ipm_cr:TestSource','Missing test dependency %s.',path);
    item=dir(path);inventory(end+1)=struct('path',path,'bytes',item.bytes,'sha256',ipm_cr.io.sha(path)); %#ok<AGROW>
end
end
