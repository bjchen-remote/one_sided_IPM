function audit_accepted(checkpoint,outputRoot,minTau,stride)
%AUDIT_ACCEPTED Recompute an accepted fixed-grid state independently.
% Use main('audit',checkpoint,newFolder).
arguments
    checkpoint (1,:) char
    outputRoot (1,:) char
    minTau (1,1) double {mustBeFinite,mustBeNonnegative} = 0
    stride (1,1) double {mustBeInteger,mustBePositive} = 1
end
assert(startsWith(checkpoint,filesep) && isfile(checkpoint));
assert(startsWith(outputRoot,filesep) && ...
    ~isfolder(outputRoot) && ~isfile(outputRoot), ...
    'Use a new absolute output folder for each audit.');
workflowRoot=fileparts(mfilename('fullpath'));
solverRoot=fileparts(fileparts(workflowRoot));
names={'IPM_AUDIT_CODE_ROOT','IPM_AUDIT_CHECKPOINT', ...
    'IPM_AUDIT_OUTPUT_DIR','IPM_AUDIT_MIN_TAU','IPM_AUDIT_STRIDE'};
old=cellfun(@getenv,names,'UniformOutput',false);
restore=onCleanup(@()restore_environment(names,old)); %#ok<NASGU>
setenv(names{1},solverRoot);
setenv(names{2},checkpoint);
setenv(names{3},outputRoot);
setenv(names{4},num2str(minTau,17));
setenv(names{5},num2str(stride));
run(fullfile(solverRoot,'audit_method','audit_checkpoint.m'));
end

function restore_environment(names,values)
for k=1:numel(names)
    setenv(names{k},values{k});
end
end
