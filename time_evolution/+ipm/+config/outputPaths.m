function paths = outputPaths(opts)
%IPM.CONFIG.OUTPUTPATHS Derive the default result and video names in one place.

anisotropic = strcmpi(string(opts.dynamicScaleGeometry),'anisotropic');
doubleOdd = strcmpi(string(opts.symmetryMode),'double_odd_omega');
dynamic = strcmpi(string(opts.rescalingMode),'dynamic');
if anisotropic && doubleOdd
    resultName = 'ipm_double_odd_anisotropic_dynamic.mat';
    videoName = 'ipm_double_odd_anisotropic_dynamic.mp4';
elseif anisotropic
    resultName = 'ipm_halfplane_anisotropic_dynamic.mat';
    videoName = 'ipm_halfplane_anisotropic_dynamic.mp4';
elseif doubleOdd && dynamic
    resultName = 'ipm_double_odd_dynamic.mat';
    videoName = 'ipm_double_odd_dynamic.mp4';
elseif doubleOdd
    resultName = 'ipm_double_odd_physical.mat';
    videoName = 'ipm_double_odd_physical.mp4';
elseif dynamic
    resultName = 'ipm_halfplane_dynamic.mat';
    videoName = 'ipm_dynamic_evolution.mp4';
else
    resultName = 'ipm_halfplane_physical.mat';
    videoName = 'ipm_physical_evolution.mp4';
end
paths = struct('resultFile',fullfile('result','production',resultName), ...
    'videoFile',fullfile('result','production',videoName));
end
