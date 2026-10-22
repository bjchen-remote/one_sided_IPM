% Prepared acceptance control; do not count it as an executed test.
restoredefaultpath;maxNumCompThreads(1);setenv('IPM_SOURCE_CFL_REL','');
testRoot=fileparts(mfilename('fullpath'));packageRoot=fileparts(testRoot);
addpath(testRoot,'-begin');addpath(packageRoot,'-begin');
assert(startsWith(which('ipm.solve'),[packageRoot,filesep]));
assert(isempty(which('ipm.field.weno5LocalMex')));
outputRoot=fullfile(fileparts(packageRoot),'complete_sdirk26_portable_acceptance_02');
assert(~isfolder(outputRoot),'The acceptance output must be a new directory.');
report=check_portable_sdirk26(outputRoot,true);
assert(report.passed && report.realTrajectorySteps==2 && ~report.goalHoldPassed);
