% Server job entry: copy this file to a job folder, then run its absolute path.
% This script loads the package from the fixed release path below.
% IPM_PREFLIGHT_ONLY=1 resolves configuration without a Poisson solve.
% IPM_VERIFY_ONLY=1 runs bounded numerical checks ON THE SERVER.
% IPM_SMOKE_ONLY=1 runs one fresh evolution step ON THE SERVER.
% IPM_PORTABLE_ACCEPTANCE=1 runs this candidate's pure MATLAB SDIRK2 acceptance.
jobRoot = fileparts(mfilename('fullpath'));
releaseRoot=jobRoot; % self-contained deployment: +ipm is beside this script
if ~isfolder(fullfile(releaseRoot,'+ipm'))
    releaseRoot='/data/user/hd58131/ipm/one_sided_ipm_wall_integral_20261007';
end
assert(isfolder(fullfile(releaseRoot,'+ipm')) && ...
    isfile(fullfile(releaseRoot,'+ipm','+evolve','wallIntegralGauge.m')), ...
    'ipm:ServerRelease','Deploy the complete new package at releaseRoot.');
restoredefaultpath;
addpath(releaseRoot,'-begin');
assert(startsWith(which('ipm.solve'),[releaseRoot,filesep]), ...
    'ipm:ServerRelease','A different +ipm package shadows this release.');

% ===== Settings to edit for a NEW run =====
s.grid.nodesX = 513;       % number of X nodes; must be odd
s.grid.nodesY = 321;       % number of Y nodes
s.box.left = 1e10;         % X runs from -left to right
s.box.right = 1e10;
s.box.top = 1e10;          % Y runs from 0 to top

% Localized positive R_X: Gaussian mixture in X times exp(-(Y/sy)^2).
% R is its analytic X primitive with R(-infinity,Y)=0.
s.initial.kind = 'localized_gaussian_primitive';
s.initial.gaussianAmplitude = .7;
s.initial.gaussianCenterX = 0;
s.initial.gaussianWidthX = .75;
s.initial.gaussianWidthY = 1;
s.initial.gaussianSkewWeight = .2; % zero selects one Gaussian
s.initial.power = 8; % inherited initial mesh template only
s.initial.cutoffB = 8; % used only by the older cutoff profile

% Keep wall Rx integral and global spline amplitude; restore the Rx^2 centroid.
s.gauge.length = 'wall_gradient_integral';
s.gauge.amplitude = 'wall_omega_global_spline_maximum';
s.gauge.position = 'wall_gradient_centroid';
s.gauge.contract = 'wall_integral_centroid_practical_v1';
s.gauge.assembly = 'net_velocity_density_centered_advective'; % explicit density-gradient WENO assembly
s.gauge.peakPositionRestoringRate = 0.5;
s.gauge.wallIntegralRadius = 1; % keep integral of R_X(X,0) over [-L,L]

% Run length and stopping target.
s.time.finalTau = Inf;     % stop on the full R_X rate target below
s.stop.fullAbsoluteRhoXRate = 1e-5;
s.stop.holdTau = 0.5;
% ===== End of settings normally edited =====

s = fixed_settings(s);
s.server.releaseRoot = releaseRoot;
outputRoot = strtrim(getenv('IPM_OUTPUT_DIR'));
if isempty(outputRoot)
    outputRoot = fullfile(jobRoot, ...
        ['output_wall_integral_',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
end
assert(startsWith(outputRoot,filesep), ...
    'ipm:ServerOutput','IPM_OUTPUT_DIR must be a new absolute directory.');
restartCheckpoint = strtrim(getenv('IPM_RESTART_CHECKPOINT'));
if strcmp(getenv('IPM_SMOKE_ONLY'),'1')
    assert(isempty(restartCheckpoint),'ipm:ServerSmoke', ...
        'Smoke mode starts a fresh run; it must not change an old checkpoint.');
    s.time.maxSteps = 1;
end
[opts,meshAudit] = ipm.workflow.build_fresh_options(s,outputRoot); %#ok<NASGU>
config = ipm.config.resolve(opts); %#ok<NASGU>
fprintf('Release: %s\nJob: %s\nOutput: %s\n',releaseRoot,jobRoot,outputRoot);
fprintf('Fresh-run grid: %d x %d; wall R_X integral radius L=%.8g\n', ...
    s.grid.nodesX,s.grid.nodesY,s.gauge.wallIntegralRadius);
if strcmp(getenv('IPM_PREFLIGHT_ONLY'),'1')
    fprintf('Preflight only: no output directory, Poisson solve or PDE step.\n');
    return
end
maxNumCompThreads(s.runtime.threads);
setenv('IPM_SOURCE_CFL_REL',num2str(s.runtime.sourceCflRel,17));
if strcmp(getenv('IPM_PORTABLE_ACCEPTANCE'),'1')
    assert(isempty(restartCheckpoint),'ipm:ServerTestRestart', ...
        'The portable acceptance starts fresh and writes its own test checkpoints.');
    addpath(fullfile(releaseRoot,'tests'),'-begin');
    acceptance = check_portable_sdirk26(outputRoot,true);
    assert(acceptance.passed,'ipm:ServerTestFailed', ...
        'The portable SDIRK2 acceptance did not pass; inspect its report.');
    return
end
if strcmp(getenv('IPM_VERIFY_ONLY'),'1')
    report = ipm.verify();
    disp(report);
    return
end
if isempty(restartCheckpoint)
    % run_fresh refuses existing output and saves settings/log/checkpoints.
    result = ipm.workflow.run_fresh(outputRoot,s);
else
    assert(startsWith(restartCheckpoint,filesep) && isfile(restartCheckpoint), ...
        'ipm:ServerRestart','Use an existing absolute checkpoint path.');
    cp = ipm.output.readCheckpoint(restartCheckpoint);
    savedAssembly='signed_generators';
    if isfield(cp.payload.state.config.scaling,'gaugeTransportAssembly')
        savedAssembly=cp.payload.state.config.scaling.gaugeTransportAssembly;
    end
    assert(strcmp(savedAssembly,s.gauge.assembly),'ipm:ServerRestartAssembly', ...
        'This checkpoint has a different numerical assembly; use an explicit new-experiment conversion or start fresh.');
    assert(strcmp(cp.payload.state.config.time.timeIntegrator,s.time.integrator), ...
        'ipm:ServerRestartIntegrator','This experimental entry resumes only its selected time integrator.');
    assert(strcmp(cp.payload.state.config.scaling.lengthGauge, ...
        'wall_gradient_integral') && ...
        strcmp(cp.payload.state.config.scaling.scalingContract, ...
        'wall_integral_centroid_practical_v1'), ...
        'ipm:ServerRestartGauge', ...
        'This entry resumes only the same practical integral/centroid contract; start fresh for a different gauge.');
    fprintf('Resume: %s; frozen L=%.8g\n',restartCheckpoint, ...
        cp.payload.state.config.scaling.wallIntegralRadius);
    clear cp
    result = ipm.workflow.resume_same_equation( ...
        restartCheckpoint,outputRoot,s.time.finalTau);
end
fprintf('Finished tau=%.8g; stop=%s\n', ...
    result.state.canonicalTime,result.state.stopReason);

function s = fixed_settings(s)
% Numerical defaults. The experiment-specific choices are above.
s.grid.adaptiveRemesh = true;
s.grid.spatialOrder = 'high_order';
s.grid.transferScheme = 'high_order';
s.grid.adaptiveLevels = [0.1,0.2,0.5];
s.transport.inflowDensityDegree = 1;
s.transport.scheme = 'weno5_fd';
s.transport.alphaMode = 'local';
s.boundary.farMode = 'green';
s.boundary.transportMode = 'open';
s.time.integrator = 'sdirk2_alexander_practical_v2';
s.time.cfl = 0.5;
s.time.maxDt = 0.2;
s.time.physicalLimit = Inf;
s.time.maxSteps = Inf;
s.time.outputEvery = 0.1;
s.time.checkpointEvery = 0.5;
s.poisson.tolerance = 1e-5;
s.poisson.solver = 'separable';
s.poisson.greenSourceTolerance = 1e-12;
s.poisson.greenMaxSources = 48000;
s.runtime.threads = 2;
s.runtime.sourceCflRel = []; % full-domain CFL; no source-mask experiment
s.output.storeSnapshots = false;
s.output.livePlot = false;
end
