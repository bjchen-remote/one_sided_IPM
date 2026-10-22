function result=main_c(action,folder,value)
%MAIN_C Direct polar Poisson with CFL-limited explicit SSPRK3 evolution.
if nargin<1,action='check';end
root=fileparts(mfilename('fullpath'));addpath(root);
cfg=ipm_c.config();
switch action
    case 'test'
        addpath(fullfile(root,'tests_ipm_c'),fullfile(root,'tests'));result=run_direct_components();
    case 'check'
        [state,context]=ipm_c.initialize(cfg);
        result=struct('grid',[context.g.nr,context.g.ntheta], ...
            'references',context.shared.references,'scale',state.scale, ...
            'method',cfg.time.timeIntegrator,'primaryState','total_rho', ...
            'explicitStep','ipm.evolve.stepRk(ssprk3)','advance','ipm_c.evolve.advance', ...
            'poisson','order-4 cubic B-spline Galerkin, FFT and cached radial LU','stepRetryUsed',false, ...
            'exponentialAuxiliarySystemUsed',false,'timeAdvanced',false);
    case {'pilot','run'}
        assert(nargin>=2,'ipm_c:OutputRequired','Supply a new absolute result folder.');
        cfg.targetTau=1e-4;if nargin>=3,cfg.targetTau=value;end
        result=ipm_c.solve(cfg,folder);
    case 'resume'
        assert(nargin>=3,'ipm_c:ResumeArguments');
        cp=ipm_c.output.loadCheckpoint(folder);
        cfg=cp.cfg;cfg.targetTau=value;
        result=ipm_c.solve(cfg,folder);
    otherwise
        error('ipm_c:Action','Use check, test, pilot, run or resume. Read IPM_C.md for validation scope.');
end
end
