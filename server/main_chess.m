% Upload this MATLAB script as the CHESS job file. Match the UI CPU allocation.
% CHESS supplies the job directory. The release remains at its absolute path.
releaseRoot = '/data/user/hd58131/ipm/ipm_cr_chess_configurable';
caseName = 'examples/bspline4.json';
targetTau = 400;
wallSeconds = 64*3600;
threads = 16;  % Set [] to use the scheduler allocation.
k = 8;          % Set [] before using the individual counts below.
nx = [];        % Cartesian X1 nodes, including ghosts. Must be odd.
ny = [];        % Cartesian X2 nodes, including ghosts.
nr = [];        % Annular radial nodes, including inward ghosts.
ntheta = [];    % Angular nodes, including both wall endpoints.
% Leave any count empty to use the selected case's mesh rule.
radialRatioLimit = [];  % Optional advisory r step ratio, e.g. 1.2.

overrides = grid_overrides(k,nx,ny,nr,ntheta,radialRatioLimit);
configuredRoot = strtrim(getenv('IPM_CHESS_RELEASE_ROOT'));
if ~isempty(configuredRoot), releaseRoot = configuredRoot; end
jobRoot = strtrim(getenv('IPM_CHESS_JOB_DIR'));
if isempty(jobRoot), jobRoot = fileparts(mfilename('fullpath')); end
configuredThreads = strtrim(getenv('IPM_CHESS_THREADS'));
if ~isempty(configuredThreads), threads = str2double(configuredThreads); end
if isempty(threads), threads = allocated_threads(); end
if strcmp(getenv('IPM_SMOKE_ONLY'),'1')
    targetTau = .002;
    wallSeconds = 180;
end
run_chess_job(releaseRoot,jobRoot,caseName,targetTau,wallSeconds,threads,overrides);

function run_chess_job(releaseRoot,jobRoot,caseName,targetTau,wallSeconds,threads,overrides)
releaseRoot = absolute_path(releaseRoot,'Release');
jobRoot = absolute_path(jobRoot,'CHESS job');
assert(isfolder(jobRoot),'ipm_cr:ChessJobDirectory','Missing CHESS job directory %s.',jobRoot);
assert(~strcmp(jobRoot,releaseRoot)&&~startsWith(jobRoot,[releaseRoot,filesep]), ...
    'ipm_cr:ChessJobDirectory','Submit the script in a CHESS job directory or set IPM_CHESS_JOB_DIR.');
validateattributes(threads,{'numeric'},{'scalar','real','finite','integer','positive'});
check_allocation(threads);

runtime = fullfile(jobRoot,'runtime');
output = fullfile(jobRoot,'native');
caseFile = fullfile(jobRoot,'case.json');
recordFile = fullfile(jobRoot,'chess_run.json');
fprintf('Release %s\nCHESS job %s\nOutput %s\nRequested MATLAB threads %d\n', ...
    releaseRoot,jobRoot,output,threads);
if strcmp(getenv('IPM_PREFLIGHT_ONLY'),'1')
    restoredefaultpath;
    if isfolder(runtime)
        addpath(fullfile(runtime,'src'),'-begin');
        opts = jsondecode(fileread(caseFile));
        verify_grid(opts,overrides);
    else
        addpath(fullfile(releaseRoot,'src'),'-begin');
        opts = read_case(releaseRoot,caseName,overrides);
    end
    [~,method] = ipm_cr.config(opts,'run');
    g = ipm_cr.mesh.build(method.geometry);
    fprintf('Mesh nodes nx=%d ny=%d nr=%d ntheta=%d (including ghosts)\n', ...
        g.core.nx,g.core.ny,g.polar.nr,g.polar.ntheta);
    if isfield(g.core,'stair')
        radial = g.axisSource.interfaceRadial;
        fprintf('Maximum radial step ratio %.8f (warning threshold %.3f)\n', ...
            radial.maxWholeAxisStepRatio,radial.adjacentStepRatioLimit);
    end
    return
end

file = java.io.RandomAccessFile(fullfile(jobRoot,'chess_owner.lock'),'rw');
channel = file.getChannel();
lock = channel.tryLock();
assert(~isempty(lock),'ipm_cr:ChessOwner','A MATLAB process already owns this CHESS job.');
owner = onCleanup(@()release_owner(lock,channel,file)); %#ok<NASGU>

restoredefaultpath;
if ~isfolder(runtime)
    assert(isfile(fullfile(releaseRoot,'src','+ipm_cr','run.m'))&& ...
        isfile(fullfile(releaseRoot,caseName)), ...
        'ipm_cr:ChessRelease','Upload the solver and case to %s.',releaseRoot);
    assert(~isfile(fullfile(output,'latest.json')), ...
        'ipm_cr:ChessRuntime','The runtime for an existing checkpoint is missing.');
    opts = read_case(releaseRoot,caseName,overrides);
    if ~isfield(opts,'run'), opts.run = struct(); end
    opts.run.tauEnd = targetTau;
    opts.run.wallSeconds = wallSeconds;
    opts.run.outputDir = output;
    opts.run.resume = true;
    mkdir(runtime);
    for name = {'src','data','main.m'}
        copy_required(fullfile(releaseRoot,name{1}),fullfile(runtime,name{1}));
    end
    mkdir(fullfile(runtime,'tools'));
    copy_required(fullfile(releaseRoot,'tools','build_ipm_cr_native.m'), ...
        fullfile(runtime,'tools','build_ipm_cr_native.m'));
    addpath(fullfile(runtime,'src'),fullfile(runtime,'tools'),'-begin');
    ipm_cr.config(opts,'run');
    native = build_ipm_cr_native(fullfile(jobRoot,'native_build'));
    ipm_cr.io.identity('json',caseFile,opts);
    record = struct('kind','ipm_cr_chess_run_v1','release',releaseRoot,'sourceCase',caseName, ...
        'threads',threads,'source',ipm_cr.io.identity('source'),'caseSha256', ...
        ipm_cr.io.identity('sha',caseFile),'nativeBuildSha256',native.binarySha256);
    ipm_cr.io.identity('json',recordFile,record);
else
    assert(isfile(recordFile)&&isfile(caseFile),'ipm_cr:ChessRuntime', ...
        'An earlier job setup is incomplete. Inspect native_build before submitting a new job.');
    addpath(fullfile(runtime,'src'),'-begin');
    record = jsondecode(fileread(recordFile));
    currentSource = ipm_cr.io.identity('source');
    assert(strcmp(record.kind,'ipm_cr_chess_run_v1')&& ...
        strcmp(record.caseSha256,ipm_cr.io.identity('sha',caseFile))&& ...
        strcmp(record.source.format,currentSource.format)&& ...
        isequaln(record.source.files(:),currentSource.files(:)), ...
        'ipm_cr:ChessRuntime','Job case or frozen runtime differs from the saved job.');
    opts = jsondecode(fileread(caseFile));
    verify_grid(opts,overrides);
    opts.run.outputDir = output;
end
rehash;
kernel = fullfile(runtime,'src','+ipm_cr','+mesh',['sampleMex.',mexext]);
assert(strcmp(which('ipm_cr.mesh.sampleMex'),kernel), ...
    'ipm_cr:ChessKernel','The job runtime must supply the platform sampling MEX.');
summary = ipm_cr.run(opts,threads);
record.lastThreads = summary.threads;
record.latestNodeCounts = struct('nx',summary.coreSize(2),'ny',summary.coreSize(1), ...
    'nr',summary.polarSize(1),'ntheta',summary.polarSize(2));
ipm_cr.io.identity('json',recordFile,record);
fprintf('Mesh nodes nx=%d ny=%d nr=%d ntheta=%d (including ghosts)\n', ...
    summary.coreSize(2),summary.coreSize(1),summary.polarSize(1),summary.polarSize(2));
fprintf('Finished tau=%.12g, status=%s, MATLAB threads=%d\n', ...
    summary.tau,summary.status,summary.threads);
end

function overrides = grid_overrides(k,nx,ny,nr,ntheta,radialRatioLimit)
overrides = struct();
names = {'nx','ny','nr','ntheta'};
values = {nx,ny,nr,ntheta};
counts = struct();
for i = 1:numel(names)
    name = names{i};
    value = configured_number(['IPM_CHESS_',upper(name)],values{i});
    if isempty(value), continue; end
    validateattributes(value,{'double'},{'real','scalar','finite','integer','>=',7}, ...
        'main_chess',name);
    counts.(name) = value;
end
if ~isempty(fieldnames(counts)), overrides.nodeCounts = counts; end
k = configured_number('IPM_CHESS_K',k);
if ~isempty(k)
    validateattributes(k,{'double'},{'real','scalar','finite','integer','>=',5,'<=',51});
    assert(isempty(fieldnames(counts)),'ipm_cr:ChessGridLevel', ...
        'Choose k or the individual grid counts.');
    overrides.k = k;
end
radialRatioLimit = configured_number('IPM_CHESS_RADIAL_RATIO_LIMIT',radialRatioLimit);
if ~isempty(radialRatioLimit)
    validateattributes(radialRatioLimit,{'double'},{'real','scalar','finite','>',1});
    overrides.radialRatioLimit = radialRatioLimit;
end
end

function value = configured_number(variable,value)
configured = strtrim(getenv(variable));
if ~isempty(configured), value = str2double(configured); end
end

function opts = read_case(releaseRoot,caseName,overrides)
path = fullfile(releaseRoot,caseName);
assert(isfile(path),'ipm_cr:ChessRelease','Missing release case %s.',path);
opts = jsondecode(fileread(path));
if isempty(fieldnames(overrides)), return; end
if ~isfield(opts,'method'), opts.method = struct(); end
if ~isfield(opts.method,'geometry')
    opts.method.geometry = struct('preset','medium','staircase',true);
end
opts.method.geometry = apply_grid(opts.method.geometry,overrides);
end

function geometry = apply_grid(geometry,overrides)
for name = fieldnames(overrides).'
    key = name{1};
    if strcmp(key,'nodeCounts')
        for axis = fieldnames(overrides.nodeCounts).'
            geometry.nodeCounts.(axis{1}) = overrides.nodeCounts.(axis{1});
        end
    else
        geometry.(key) = overrides.(key);
    end
end
end

function verify_grid(opts,overrides)
if isempty(fieldnames(overrides)), return; end
assert(isfield(opts,'method')&&isfield(opts.method,'geometry')&& ...
    isequaln(opts.method.geometry,apply_grid(opts.method.geometry,overrides)), ...
    'ipm_cr:ChessFrozenGrid','To change grid settings, submit a new CHESS job directory.');
end

function threads = allocated_threads()
limits = [];
for key = {'SLURM_CPUS_PER_TASK','PBS_NCPUS'}
    value = strtrim(getenv(key{1}));
    if ~isempty(value), limits(end+1) = str2double(value); end %#ok<AGROW>
end
assert(~isempty(limits),'ipm_cr:ChessAllocation', ...
    'No scheduler CPU count is visible. Set threads or IPM_CHESS_THREADS to the CHESS UI allocation.');
threads = min(limits);
end

function path = absolute_path(value,label)
assert(~isempty(value)&&java.io.File(value).isAbsolute(), ...
    'ipm_cr:ChessAbsolutePath','%s directory must be an absolute path.',label);
path = char(java.io.File(value).getCanonicalPath());
end

function check_allocation(threads)
for name = {'SLURM_CPUS_PER_TASK','PBS_NCPUS'}
    value = strtrim(getenv(name{1}));
    if ~isempty(value)
        assert(str2double(value)>=threads,'ipm_cr:ChessAllocation', ...
            '%s=%s cannot provide %d MATLAB threads.',name{1},value,threads);
    end
end
if isunix && ~ismac
    [status,count] = system('nproc');
    if status==0
        assert(str2double(strtrim(count))>=threads,'ipm_cr:ChessAllocation', ...
            'The job exposes %s CPUs, requested %d.',strtrim(count),threads);
    end
end
end

function copy_required(source,destination)
[ok,message] = copyfile(source,destination);
assert(ok,'ipm_cr:ChessCopy','%s',message);
end

function release_owner(lock,channel,file)
lock.release();
channel.close();
file.close();
end
