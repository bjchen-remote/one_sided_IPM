function varargout = main(action,varargin)
%MAIN Edit the settings below, then run main in MATLAB.
% main('check') checks the mesh, initial density and C conditions without solving.
% main('run','/absolute/output/folder') starts a fresh run.
% main('resume',checkpoint,newFolder,finalTau) continues an accepted run.

% ===== Settings to edit for a NEW run =====
s.grid.nodesX = 513;       % number of X nodes; must be odd
s.grid.nodesY = 321;       % number of Y nodes
s.box.left = 1e10;         % X runs from -left to right
s.box.right = 1e10;
s.box.top = 1e10;          % Y runs from 0 to top

% Localized Gaussian R_X; R is its analytic primitive with left density zero.
s.initial.kind = 'localized_gaussian_primitive';
s.initial.gaussianAmplitude = .7;
s.initial.gaussianCenterX = 0;
s.initial.gaussianWidthX = .75;
s.initial.gaussianWidthY = 1;
s.initial.gaussianSkewWeight = .2;
s.initial.power = 8;
s.initial.cutoffB = 8;     % only used for the ..._cutoff variant

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
root = fileparts(mfilename('fullpath'));
solverRoot = root;
addpath(solverRoot,'-begin');
assert(startsWith(which('ipm.solve'),[solverRoot,filesep]), ...
    'The solver package is not available from this folder.');
if nargin < 1 || isempty(action), action = 'run'; end
assert(nargout<=1,'main returns at most one result.');
action = validatestring(char(action),{'check','smoke','run','resume'});
result = [];
switch action
    case 'check'
        assert(isempty(varargin),'main(''check'') takes no other arguments.');
        ipm.workflow.self_check(s);
    case {'run','smoke'}
        assert(numel(varargin)<=1,'Pass at most one output folder.');
        if isempty(varargin)
            outputRoot = fullfile(root,'runs', ...
                ['run_',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
        else
            outputRoot = char(varargin{1});
        end
        if strcmp(action,'smoke')
            s.time.maxSteps = 1;
        end
        ipm.workflow.self_check(s);
        result = ipm.workflow.run_fresh(outputRoot,s);
    case 'resume'
        assert(numel(varargin)==3, ...
            'Use main(''resume'',checkpoint,newFolder,finalTau).');
        cp=ipm.output.readCheckpoint(char(varargin{1}));
        assert(strcmp(cp.payload.state.config.scaling.gaugeTransportAssembly,s.gauge.assembly), ...
            'ipm:MainRestartAssembly', ...
            'This entry resumes its selected assembly; use the development same-equation workflow for a historical operator.');
        clear cp
        result = ipm.workflow.resume_same_equation(char(varargin{1}), ...
            char(varargin{2}),varargin{3});
end
if nargout > 0
    varargout{1} = result;
end
end

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
 
