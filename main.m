function varargout = main(action,varargin)
%MAIN Edit the settings below, then run main in MATLAB.
% main('check') checks the mesh, initial density and C conditions without solving.
% main('run','/absolute/output/folder') starts a fresh run.
% main('resume',checkpoint,newFolder,finalTau) continues an accepted run.
% main('audit',checkpoint,newFolder) audits an accepted checkpoint.

% ===== Settings to edit for a NEW run =====
s.grid.nodesX = 401;       % number of X nodes; must be odd
s.grid.nodesY = 321;       % number of Y nodes
s.box.left = 1e5;          % X runs from -left to right
s.box.right = 1e5;
s.box.top = 1e5;           % Y runs from 0 to top

% Initial density: log(1+max(X+a,0)^k/(1+Y^k))/k, a=(k-1)^(1/k).
s.initial.kind = 'shifted_degenerate_primitive';
s.initial.power = 8;
s.initial.cutoffB = 8;     % only used for the ..._cutoff variant

% Three C/gauge conditions; the solver computes c_l, c_omega and c_r.
s.gauge.length = 'wall_gradient_curvature';
s.gauge.amplitude = 'wall_omega_quadratic_peak';
s.gauge.position = 'gradient_peak';
s.gauge.contract = 'peak_curvature_restoring_v1';
s.gauge.peakPositionRestoringRate = 0.5;
s.gauge.peakCurvatureRadius = 0.5;

% Run length and stopping target.
s.time.finalTau = 20;
s.stop.fullAbsoluteRhoXRate = 1e-4;
s.stop.holdTau = 0.5;
% ===== End of settings normally edited =====

s = fixed_settings(s);
root = fileparts(mfilename('fullpath'));
solverRoot = fullfile(root,'solver');
addpath(solverRoot,'-begin');
assert(startsWith(which('ipm.solve'),[solverRoot,filesep]), ...
    'The solver package is not available from this folder.');
if nargin < 1 || isempty(action), action = 'run'; end
assert(nargout<=1,'main returns at most one result.');
action = validatestring(char(action),{'check','smoke','run','resume','audit'});
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
        result = ipm.workflow.resume_same_equation(char(varargin{1}), ...
            char(varargin{2}),varargin{3});
    case 'audit'
        assert(numel(varargin)>=2 && numel(varargin)<=4, ...
            'Use main(''audit'',checkpoint,newFolder[,minTau,stride]).');
        args = varargin;
        args{1} = char(args{1});
        args{2} = char(args{2});
        ipm.workflow.audit_accepted(args{:});
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
s.time.integrator = 'ssprk3';
s.time.cfl = 0.5;
s.time.maxDt = 0.01;
s.time.physicalLimit = 40;
s.time.maxSteps = 100000;
s.time.outputEvery = 0.1;
s.time.checkpointEvery = 0.5;
s.poisson.tolerance = 1e-5;
s.poisson.solver = 'separable';
s.poisson.greenSourceTolerance = 1e-8;
s.poisson.greenMaxSources = 24000;
s.runtime.threads = 4;
s.runtime.sourceCflRel = 1e-10;
end
 