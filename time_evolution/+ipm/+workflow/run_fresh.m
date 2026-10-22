function result = run_fresh(outputRoot,s)
%RUN_FRESH Start a NEW short validation, preserving every existing result.
% Called by main.m with its settings structure.
arguments
    outputRoot (1,:) char
    s (1,1) struct
end
assert(startsWith(outputRoot,filesep), ...
    'Use an absolute output directory to make provenance unambiguous.');
assert(~isfolder(outputRoot) && ~isfile(outputRoot), ...
    'Choose a NEW output directory; existing results are never overwritten.');
[opts,meshAudit]=ipm.workflow.build_fresh_options(s,outputRoot);
workflowRoot=fileparts(mfilename('fullpath'));
solverRoot=fileparts(fileparts(workflowRoot));
addpath(solverRoot,'-begin');
assert(startsWith(which('ipm.solve'),[solverRoot,filesep]));
ipm.config.resolve(opts); % Reject incompatible C/gauge choices before I/O.
maxNumCompThreads(s.runtime.threads);
setenv('IPM_SOURCE_CFL_REL',num2str(s.runtime.sourceCflRel,17));
mkdir(outputRoot);
save(fullfile(outputRoot,'settings_used.mat'),'s','opts','meshAudit','-v7.3');
diary(fullfile(outputRoot,'run.log'));
diaryGuard=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('output=%s\n',outputRoot);
result=ipm.solve(opts);
end
