function opts = activeCase()
%IPM.CONFIG.ACTIVECASE No-argument practical half-plane production settings.
%   Use the same analytic multiscale axes and options builder as main.m.
%   This prepares configuration only: no PDE solve or file write is started.
%   Default artifacts receive distinct names in the existing current folder.

s.grid = struct('nodesX',513,'nodesY',321,'adaptiveRemesh',true, ...
    'spatialOrder','high_order','transferScheme','high_order', ...
    'adaptiveLevels',[.1,.2,.5],'maximumCellRatio',6);
s.box = struct('left',1e10,'right',1e10,'top',1e10);
s.initial = struct('kind','shifted_degenerate_primitive','power',8,'cutoffB',8);
s.gauge = struct('length','wall_gradient_integral', ...
    'assembly','net_velocity_density_centered_advective', ...
    'amplitude','wall_omega_global_spline_maximum', ...
    'position','wall_gradient_centroid','contract','wall_integral_centroid_practical_v1', ...
    'peakPositionRestoringRate',.5,'wallIntegralRadius',1);
s.time = struct('finalTau',Inf,'integrator','ssprk3','cfl',.5,'maxDt',.01, ...
    'physicalLimit',Inf,'maxSteps',Inf,'outputEvery',.1,'checkpointEvery',.5);
s.stop = struct('fullAbsoluteRhoXRate',1e-5,'holdTau',.5);
s.transport = struct('inflowDensityDegree',1,'scheme','weno5_fd','alphaMode','local');
s.boundary = struct('farMode','green','transportMode','open');
s.poisson = struct('tolerance',1e-5,'solver','separable', ...
    'greenSourceTolerance',1e-12,'greenMaxSources',48000);
s.runtime = struct('threads',2,'sourceCflRel',[]);
s.output = struct('storeSnapshots',false,'livePlot',false);

outputRoot = pwd;
[opts,~] = ipm.workflow.build_fresh_options(s,outputRoot);
[~,token] = fileparts(tempname);
stem = ['ipm_practical_',token];
opts.resultFile = fullfile(outputRoot,[stem,'_result.mat']);
opts.checkpoint.file = fullfile(outputRoot,[stem,'_checkpoint.mat']);
opts.caseMetadata.fullRhoXRateMonitorFile = fullfile(outputRoot,[stem,'_full_rate.csv']);
opts.caseMetadata.study = 'practical_default_active_case';
end
