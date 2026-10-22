restoredefaultpath;maxNumCompThreads(1);setenv('IPM_SOURCE_CFL_REL','');
testRoot=fileparts(mfilename('fullpath'));root=fileparts(testRoot);validation=fileparts(root);
addpath(testRoot,'-begin');addpath(root,'-begin');
parentReport=jsondecode(fileread(fullfile(validation,'complete_sdirk26_portable_acceptance_02','report.json')));
manifest=jsondecode(fileread(fullfile(root,'MANIFEST.json')));
assert(parentReport.passed && isfile(parentReport.checkpointFile) && ...
    strcmp(parentReport.timeIntegrator,'sdirk2_alexander_practical_v2') && ...
    strcmp(parentReport.initialCondition,'localized_gaussian_primitive') && ...
    strcmp(parentReport.runtimeTreeSha256,manifest.runtimeTreeSha256));
outputRoot=fullfile(validation,'complete_sdirk26_portable_io_acceptance_02');
report=check_portable_plot_fields26(outputRoot,parentReport.checkpointFile);
assert(report.passed && ~report.goalHoldPassed);
