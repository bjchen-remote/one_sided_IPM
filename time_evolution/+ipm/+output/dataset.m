function stream = dataset(action,varargin)
%IPM.OUTPUT.DATASET Publish a compact catalogue and immutable native frames.
%   STREAM = DATASET('open',STATE,LOG,POLICY,RESTARTINPUT) prepares an
%   observation-only writer. POLICY has enabled, file and frameEvery fields.
%   STREAM = DATASET('record',STREAM,STATE,LOG,FORCE,REASON,ANALYSIS) records
%   one accepted state. No solver state, checkpoint, or numerical option is
%   changed. A frame journal is durable before its catalogue is published.
%   STREAM = DATASET('seed',PACKEDFILE,POLICY) copies a schema-2 packed
%   archive into a new native journal. POLICY.file is required. It performs
%   no RHS/operator reconstruction; clear STREAM to release its writer lock
%   before continuing the solver. Repeating the same interrupted seed is safe.

switch lower(char(action))
    case 'open'
        stream = open_stream(varargin{:});
    case 'record'
        stream = record_state(varargin{:});
    case 'seed'
        stream = seed_archive(varargin{:});
    otherwise
        error('ipm:DatasetAction','Unknown dataset action "%s".',char(action));
end
end

function stream = seed_archive(packedFile,policy)
if nargin<2 || ~isstruct(policy) || ~isscalar(policy) || ...
        ~isfield(policy,'file') || strlength(string(policy.file))==0
    error('ipm:DatasetSeedPolicy','Seeding requires an explicit new policy.file.');
end
if ~isfield(policy,'enabled'), policy.enabled = true; end
if ~isequal(policy.enabled,true)
    error('ipm:DatasetSeedPolicy','Seeding requires enabled=true.');
end
sourceFile = absolute_path(char(packedFile));
targetFile = absolute_path(char(policy.file));
if strcmp(sourceFile,targetFile)
    error('ipm:DatasetSeedPath','The native catalogue must differ from its read-only packed source.');
end
sourceHash = file_hash(sourceFile);
loaded = load(sourceFile,'data');
if ~isfield(loaded,'data') || ~isstruct(loaded.data) || ~isscalar(loaded.data)
    error('ipm:DatasetSeedArchive','The packed source must contain one scalar data structure.');
end
archive = loaded.data;
clear loaded
if ~all(isfield(archive,{'schemaVersion','meta','diagnostics','frames','tau','step', ...
        'fullRateNorm','coreRateNorm'})) || archive.schemaVersion~=2 || ...
        ~all(isfield(archive.meta,{'caseId','storage','config','runMetadata','sourceCommonNames'})) || ...
        ~strcmp(char(archive.meta.storage),'packed_native_frames') || ...
        isempty(archive.frames) || ~all(isfield(archive.frames,'payload'))
    error('ipm:DatasetSeedArchive','Seeding requires a complete schema-2 packed_native_frames archive.');
end
history = archive.diagnostics;
if ~isfield(history,'common') || ~isstruct(history.common)
    error('ipm:DatasetSeedArchive','The packed archive lacks diagnostics.common.');
end
aliases = archive.meta.sourceCommonNames;
if ~isstruct(aliases) || ~isscalar(aliases)
    error('ipm:DatasetSeedArchive','sourceCommonNames must map raw history names to plotting series.');
end
names = fieldnames(aliases);
for k = 1:numel(names)
    raw = names{k}; target = char(aliases.(raw));
    if ~isfield(archive,target) || numel(archive.(target))~=numel(archive.tau)
        error('ipm:DatasetSeedArchive','Cannot restore raw history field %s from data.%s.',raw,target);
    end
    if isfield(history.common,raw) && ...
            ~isequaln(history.common.(raw)(:),archive.(target)(:))
        error('ipm:DatasetSeedArchive','Raw history alias %s disagrees with its plotting series.',raw);
    end
    history.common.(raw) = archive.(target)(:);
end
log = struct('history',history);
state = struct('config',archive.meta.config,'runMetadata',archive.meta.runMetadata, ...
    'scale',struct('canonicalTime',archive.tau(end)),'step',archive.step(end));
if ~strcmp(char(archive.meta.caseId),char(state.runMetadata.caseId))
    error('ipm:DatasetSeedIdentity','Packed archive caseId differs from its runtime metadata.');
end
if isfile(targetFile)
    previous = load(targetFile,'data');
    if ~isfield(previous,'data') || ~isfield(previous.data,'meta') || ...
            ~isfield(previous.data.meta,'seedArchiveSha256') || ...
            ~strcmp(char(previous.data.meta.seedArchiveSha256),sourceHash)
        error('ipm:DatasetSeedExists', ...
            'Seed target already exists and was not published from this archive: %s.',targetFile);
    end
    clear previous
end
stream = open_stream(state,log,policy,[]);
markerFile = fullfile(stream.frameFolder,'seed_archive.mat');
seed = struct('schemaVersion',1,'caseId',stream.caseId, ...
    'fingerprint',stream.fingerprint,'archiveFile',sourceFile, ...
    'archiveSha256',sourceHash);
if isfile(markerFile)
    previous = load(markerFile,'seed');
    if ~isfield(previous,'seed') || ...
            ~all(isfield(previous.seed,{'caseId','fingerprint','archiveSha256'})) || ...
            ~strcmp(char(previous.seed.caseId),seed.caseId) || ...
            ~strcmp(char(previous.seed.fingerprint),seed.fingerprint) || ...
            ~strcmp(char(previous.seed.archiveSha256),seed.archiveSha256)
        error('ipm:DatasetSeedIdentity','Partial seed journal belongs to another source archive: %s.',markerFile);
    end
else
    if ~isempty(stream.frames)
        error('ipm:DatasetSeedExists','A new seed cannot adopt an existing non-seed frame journal.');
    end
    temporaryFile = [tempname(stream.frameFolder),'.mat'];
    cleanup = onCleanup(@()delete_temporary(temporaryFile));
    save(temporaryFile,'seed','-v7.3');
    [ok,message] = movefile(temporaryFile,markerFile);
    if ~ok, error('ipm:DatasetSeedMove','Cannot install seed identity: %s',message); end
end
stream.seedProvenance = seed;
stream.seedProvenance.sourceCommonNames = aliases;
stream.seedProvenance.copyDefinition = ...
    'Bit-preserving copy of packed accepted fields; no PDE, RHS, derivative, interpolation, or operator reconstruction';
metadata = archive.meta;
removed = intersect(fieldnames(metadata),{'config','runMetadata','sourceCommonNames', ...
    'datasetFile','baseFolder'});
stream.seedProvenance.sourceMetadata = rmfield(metadata,removed);
if isfield(archive,'references'), stream.references = archive.references; end
stream.references.seedArchive = struct('file',sourceFile,'sha256',sourceHash, ...
    'kind','read_only_packed_native_seed');
if isfield(archive,'analysis')
    stream.analysis.seedArchive=struct('tau',archive.tau, ...
        'fieldTau',reshape([archive.frames.tau],[],1), ...
        'publishedThroughTau',archive.tau(end),'data',archive.analysis);
end
stream = import_extras(stream,archive);
count = numel(archive.tau);
if any([numel(archive.step),numel(archive.fullRateNorm),numel(archive.coreRateNorm)]~=count)
    error('ipm:DatasetSeedArchive','Packed scalar observation vectors are misaligned.');
end
recoveredObservations = stream.observations;
stream.observations = empty_observations();
for k = 1:count
    if isfinite(archive.fullRateNorm(k))
        stream.observations = merge_observation(stream.observations,archive.tau(k), ...
            archive.step(k),archive.fullRateNorm(k),archive.coreRateNorm(k));
    end
end
% The packed archive's first scalar measurements remain authoritative even
% when resuming a partial journal whose sparse frame norms differ by a last
% bit. Check compatibility without replacing those original measurements.
for k = 1:numel(recoveredObservations.tau)
    stream.observations = merge_observation(stream.observations, ...
        recoveredObservations.tau(k),recoveredObservations.step(k), ...
        recoveredObservations.fullRateNorm(k),recoveredObservations.coreRateNorm(k));
end
standard = {'schemaVersion','tau','t','step','x','y','rho','w','rho_t','rho_xt', ...
    'rho_y','psi','psiLow','u1','u2','transport_u1','transport_u2', ...
    'Cx','Cy','Comega','Xshift','remeshCount','fullRateNorm','coreRateNorm','caseId'};
frameMetadata = repmat(struct('tau',0,'step',0,'extras',struct()),numel(archive.frames),1);
for k = 1:numel(archive.frames)
    payload = archive.frames(k).payload;
    if ~isstruct(payload) || ~isscalar(payload) || ~all(isfield(payload,standard))
        error('ipm:DatasetSeedFrame','Packed frame %d lacks the standard accepted-frame contract.',k);
    end
    frame = struct();
    for j = 1:numel(standard), frame.(standard{j}) = payload.(standard{j}); end
    assert_frame_key(frame,archive.frames(k).tau,archive.frames(k).step);
    if ~strcmp(char(frame.caseId),stream.caseId)
        error('ipm:DatasetSeedIdentity','Packed frame %d belongs to another case.',k);
    end
    index = history_index(frame.tau,frame.step,log);
    if beyond_state(frame.tau,frame.step,state)
        error('ipm:DatasetSeedFrame','Packed frame %d exceeds its scalar history.',k);
    end
    extras = rmfield(payload,standard);
    frameMetadata(k) = struct('tau',frame.tau,'step',frame.step,'extras',extras);
    [descriptor,journal,check] = write_frame(stream,frame,index,'seed_archive_native_copy');
    [stream.frames,mergedCheck] = merge_descriptor(stream.frames,descriptor);
    stream.replayChecks = merge_replay_check(stream.replayChecks,check);
    stream.replayChecks = merge_replay_check(stream.replayChecks,mergedCheck);
    stream.observations = merge_observation(stream.observations,frame.tau,frame.step, ...
        journal.fullRateNorm,journal.coreRateNorm);
end
[~,order] = sortrows([[stream.frames.tau]',[stream.frames.step]'],[1,2]);
stream.frames = stream.frames(order);
stream.seedProvenance.frameMetadata = frameMetadata;
stream.lastFrameRemesh = archive.frames(end).payload.remeshCount;
stream.nextFrameTau = next_frame_tau(stream.frames(end).tau,stream.frameEvery);
% Hash again before publication to ensure the read-only source stayed fixed.
if ~strcmp(file_hash(sourceFile),sourceHash)
    error('ipm:DatasetSeedSourceChanged','The packed source changed while seeding: %s.',sourceFile);
end
publish_catalog(stream,state,log,'seeded_from_verified_archive');
end

function stream = open_stream(state,log,policy,restartInput)
if nargin < 3 || isempty(policy), policy = struct(); end
if nargin < 4, restartInput = []; end
if ~isstruct(policy) || ~isscalar(policy)
    error('ipm:DatasetPolicy','Dataset policy must be a scalar structure.');
end
unexpected = setdiff(fieldnames(policy),{'enabled','file','frameEvery'});
if ~isempty(unexpected)
    error('ipm:DatasetPolicy','Unknown dataset policy field: %s.',unexpected{1});
end
if ~isfield(policy,'enabled'), policy.enabled = state.config.output.saveResults; end
if ~(islogical(policy.enabled) && isscalar(policy.enabled))
    error('ipm:DatasetPolicy','Dataset enabled must be a logical scalar.');
end
stream = struct('enabled',policy.enabled);
if ~stream.enabled, return; end
defaultFile = ~isfield(policy,'file') || isempty(policy.file) || ...
    (isstring(policy.file) && isscalar(policy.file) && strlength(policy.file)==0);
if defaultFile
    resultFile = state.config.output.resultFile;
    if isfield(state.runMetadata,'resultFile') && ...
            strlength(string(state.runMetadata.resultFile)) > 0
        resultFile = state.runMetadata.resultFile;
    end
    policy.file = fullfile(fileparts(char(resultFile)),'data.mat');
end
if ~isfield(policy,'frameEvery')
    policy.frameEvery = 0.5;
    if isfield(state.config.output,'checkpoint') && ...
            isfinite(state.config.output.checkpoint.every)
        policy.frameEvery = state.config.output.checkpoint.every;
    end
end
validateattributes(policy.frameEvery,{'double'}, ...
    {'scalar','real','finite','positive'},mfilename,'policy.frameEvery');
if ~(ischar(policy.file) && isrow(policy.file) && ~isempty(policy.file)) && ...
        ~(isstring(policy.file) && isscalar(policy.file) && strlength(policy.file)>0)
    error('ipm:DatasetPolicy','Dataset file must be a nonempty path.');
end
file = absolute_path(char(policy.file));
if defaultFile && isfile(file)
    existing = load(file,'data');
    sameCase = isfield(existing,'data') && isstruct(existing.data) && ...
        isscalar(existing.data) && isfield(existing.data,'meta') && ...
        isstruct(existing.data.meta) && isscalar(existing.data.meta) && ...
        isfield(existing.data.meta,'caseId') && ...
        strcmp(char(existing.data.meta.caseId),char(state.runMetadata.caseId));
    if ~sameCase
        safeCaseId = regexprep(char(state.runMetadata.caseId),'[^A-Za-z0-9_.-]','_');
        file = fullfile(fileparts(file),['data_',safeCaseId,'.mat']);
    end
    clear existing
end
[folder,stem,extension] = fileparts(file);
if isempty(stem) || ~strcmpi(extension,'.mat')
    error('ipm:DatasetPolicy','Dataset file must name a .mat file.');
end
if ~isfolder(folder)
    [ok,message] = mkdir(folder);
    if ~ok, error('ipm:DatasetDirectory','Cannot create %s: %s',folder,message); end
end
[writerLock,writerCleanup] = acquire_writer_lock(file);
frameFolder = fullfile(folder,[stem,'_frames']);
if ~isfolder(frameFolder)
    [ok,message] = mkdir(frameFolder);
    if ~ok, error('ipm:DatasetDirectory','Cannot create %s: %s',frameFolder,message); end
end
stream = struct('enabled',true,'file',file,'baseFolder',folder, ...
    'frameFolder',frameFolder,'frameEvery',policy.frameEvery, ...
    'caseId',char(state.runMetadata.caseId), ...
    'numericsSha256',numerics_hash(), ...
    'fingerprint','', ...
    'writerLock',writerLock,'writerCleanup',writerCleanup, ...
    'frames',empty_descriptors(),'observations',empty_observations(), ...
    'replayChecks',empty_replay_checks(), ...
    'seedProvenance',struct(),'references',struct(), ...
    'seriesExtras',struct('tau',[],'step',[],'values',struct()), ...
    'diagnosticExtras',struct(), ...
    'analysis',struct(),'lastFrameRemesh',NaN,'nextFrameTau',0);
stream.fingerprint = numerical_fingerprint(state.config,stream.numericsSha256);

% A local catalogue is authoritative. A new continuation directory may
% inherit an earlier catalogue without changing any source files.
if isfile(file)
    stream = import_catalog(stream,file,state,log,false);
else
    source = restart_catalog(restartInput,stream.caseId,stream.fingerprint,state);
    if ~isempty(source) && ~strcmp(source,file)
        stream = import_catalog(stream,source,state,log,true);
    end
end
stream = recover_journal(stream,state,log);
if ~isempty(stream.frames)
    [~,order] = sortrows([[stream.frames.tau]',[stream.frames.step]'],[1,2]);
    stream.frames = stream.frames(order);
    latest = read_frame(stream.frames(end));
    if isfield(latest,'remeshCount'), stream.lastFrameRemesh = latest.remeshCount; end
    stream.nextFrameTau = next_frame_tau(stream.frames(end).tau,stream.frameEvery);
end
end

function stream = record_state(stream,state,log,force,reason,analysis)
if ~stream.enabled, return; end
assert_writer_lock(stream);
if nargin < 5 || isempty(reason), reason = 'running'; end
if nargin < 4 || isempty(force), force = false; end
if nargin < 6, analysis = []; end
if ~(islogical(force) && isscalar(force))
    error('ipm:DatasetRecord','Dataset force must be a logical scalar.');
end
if ~(ischar(reason) && isrow(reason)) && ...
        ~(isstring(reason) && isscalar(reason))
    error('ipm:DatasetRecord','Dataset reason must be text.');
end
reason = char(reason);
if ~strcmp(stream.caseId,char(state.runMetadata.caseId))
    error('ipm:DatasetIdentity','Dataset writer belongs to another case.');
end
[tau,step,historyIndex] = accepted_key(state,log);
validate_cache(state);
tolerance = clock_tolerance(tau);
if ~isempty(stream.observations.tau) && ...
        (tau < stream.observations.tau(end)-tolerance || ...
         step < stream.observations.step(end))
    error('ipm:DatasetBranch','Cannot append an earlier state to this dataset.');
end
same = find([stream.frames.tau]'==tau & [stream.frames.step]'==step,1);
saveFrame = ~isempty(same) || isempty(stream.frames) || force || ...
    tau+tolerance>=stream.nextFrameTau || ...
    state.ops.remeshCount~=stream.lastFrameRemesh;
if saveFrame
    frame = ipm.output.dataFrame(state);
    assert_frame_key(frame,tau,step);
    if ~isempty(same)
        previous = read_frame(stream.frames(same));
        check = compare_frames(previous,frame);
        stream.replayChecks = merge_replay_check(stream.replayChecks,check);
        fullRateNorm = previous.fullRateNorm;
        coreRateNorm = previous.coreRateNorm;
    else
        [descriptor,journal,check] = write_frame(stream,frame,historyIndex,reason);
        stream.frames(end+1,1) = descriptor;
        stream.replayChecks = merge_replay_check(stream.replayChecks,check);
        fullRateNorm = journal.fullRateNorm;
        coreRateNorm = journal.coreRateNorm;
    end
    stream.lastFrameRemesh = state.ops.remeshCount;
    stream.nextFrameTau = next_frame_tau(tau,stream.frameEvery);
else
    rate = state.rhsCache.rhoRate*state.ops.Dx';
    fullRateNorm = max(abs(rate),[],'all');
    coreX = state.ops.x>=-2 & state.ops.x<=2;
    coreY = state.ops.y>=0 & state.ops.y<=2;
    coreRateNorm = NaN;
    if any(coreX) && any(coreY)
        coreRateNorm = max(abs(rate(coreY,coreX)),[],'all');
    end
end
stream.observations = merge_observation(stream.observations, ...
    tau,step,fullRateNorm,coreRateNorm);
if ~isempty(analysis)
    if ~isstruct(analysis) || ~isscalar(analysis)
        error('ipm:DatasetAnalysis','Dataset analysis updates must be a scalar structure.');
    end
    names = fieldnames(analysis);
    for k = 1:numel(names), stream.analysis.(names{k}) = analysis.(names{k}); end
end
publish_catalog(stream,state,log,reason);
end

function publish_catalog(stream,state,log,reason)
tau = state.scale.canonicalTime;
step = state.step;
data = ipm.output.dataSeries(state,log,stream.observations);
if ~isfield(data,'meta') || ~isstruct(data.meta) || ~isscalar(data.meta)
    error('ipm:DatasetSeries','dataSeries must return a scalar data.meta structure.');
end
data.frames = public_descriptors(stream.frames,stream.baseFolder);
data.fieldTau = reshape([stream.frames.tau],[],1);
data.fieldStep = reshape([stream.frames.step],[],1);
data.analysis = stream.analysis;
if ~isempty(fieldnames(stream.references)), data.references = stream.references; end
names = fieldnames(stream.diagnosticExtras);
for k = 1:numel(names)
    if ~isfield(data.diagnostics,names{k})
        data.diagnostics.(names{k}) = stream.diagnosticExtras.(names{k});
    end
end
names = fieldnames(stream.seriesExtras.values);
if ~isempty(names)
    [found,rows] = ismember([stream.seriesExtras.tau,stream.seriesExtras.step], ...
        [data.tau,data.step],'rows');
    for k = 1:numel(names)
        name = names{k};
        values = NaN(numel(data.tau),1);
        source = stream.seriesExtras.values.(name);
        values(rows(found)) = source(found);
        data.(name) = values;
    end
end
data.meta.caseId = stream.caseId;
data.meta.datasetFile = stream.file;
data.meta.baseFolder = stream.baseFolder;
data.meta.storage = 'native_frame_journal';
data.meta.status = reason;
data.meta.fingerprint = stream.fingerprint;
data.meta.numericsSha256 = stream.numericsSha256;
data.meta.replayChecks = stream.replayChecks;
data.meta.fingerprintDefinition = ...
    ['SHA256 canonical config JSON excluding output and time.finalTime/maxSteps/physicalFinalTime, ' ...
    'combined with the numerical MATLAB source tree SHA256'];
data.meta.frameEvery = stream.frameEvery;
data.meta.publishedThroughTau = tau;
data.meta.publishedThroughStep = step;
data.meta.frameCount = numel(stream.frames);
data.meta.nativeGridPreserved = true;
data.meta.spatialInterpolationUsed = false;
data.meta.timeInterpolationUsed = false;
if ~isempty(fieldnames(stream.seedProvenance))
    data.meta.seedProvenance = stream.seedProvenance;
    data.meta.seedArchiveSha256 = stream.seedProvenance.archiveSha256;
    data.meta.seedArchiveFile = stream.seedProvenance.archiveFile;
    if isfield(stream.seedProvenance,'sourceMetadata')
        for name = {'widthDefinitions','monitorFieldAliases'}
            if isfield(stream.seedProvenance.sourceMetadata,name{1})
                data.meta.(name{1}) = stream.seedProvenance.sourceMetadata.(name{1});
            end
        end
    end
end
if ~isfield(data.meta,'definitions'), data.meta.definitions = struct(); end
data.meta.definitions.frameSampling = ...
    'initial accepted state, frameEvery threshold, remesh changes, and forced accepted terminal states';
data.meta.definitions.rho_xt = ...
    'partial_tau(R_X) at fixed native X,Y; accepted density RHS times Dx transpose';
data.meta.definitions.frameJournal = ...
    'immutable MAT files containing frame and journal; catalogue publication is atomic';
assert_writer_lock(stream);
save_catalog(stream.file,data);
end

function stream = import_catalog(stream,file,state,log,inherited)
loaded = load(file,'data');
if ~isfield(loaded,'data') || ~isstruct(loaded.data) || ~isscalar(loaded.data)
    error('ipm:DatasetCatalog','Dataset catalogue %s has no scalar data structure.',file);
end
data = loaded.data;
if ~isfield(data,'meta') || ...
        ~all(isfield(data.meta,{'caseId','fingerprint','storage'})) || ...
        ~strcmp(char(data.meta.caseId),stream.caseId) || ...
        ~strcmp(char(data.meta.fingerprint),stream.fingerprint) || ...
        ~strcmp(char(data.meta.storage),'native_frame_journal')
    error('ipm:DatasetIdentity', ...
        'Refusing to combine %s: caseId, numerical fingerprint, or native storage contract differs.',file);
end
if ~isfield(data,'frames') || ~isstruct(data.frames)
    error('ipm:DatasetCatalog','Dataset %s has no frame descriptor array.',file);
end
assert_history_prefix(data,log,state,inherited,file);
sourceFolder = fileparts(file);
frames = data.frames(:);
for k = 1:numel(frames)
    if inherited && all(isfield(frames(k),{'tau','step'})) && ...
            beyond_state(frames(k).tau,frames(k).step,state)
        continue;
    end
    descriptor = normalize_descriptor(frames(k),sourceFolder);
    if ~key_in_history(descriptor.tau,descriptor.step,log)
        if inherited && beyond_state(descriptor.tau,descriptor.step,state), continue; end
        error('ipm:DatasetBranch','Frame step %d tau %.17g is outside this accepted history.', ...
            descriptor.step,descriptor.tau);
    end
    if beyond_state(descriptor.tau,descriptor.step,state)
        if inherited, continue; end
        error('ipm:DatasetBranch','Existing catalogue %s is ahead of the restart state.',file);
    end
    actualHash = file_hash(descriptor.file);
    if ~strcmp(actualHash,descriptor.sha256)
        error('ipm:DatasetIntegrity','Frame SHA256 mismatch: %s.',descriptor.file);
    end
    descriptor.historyIndex = history_index(descriptor.tau,descriptor.step,log);
    [stream.frames,check] = merge_descriptor(stream.frames,descriptor);
    stream.replayChecks = merge_replay_check(stream.replayChecks,check);
end
if all(isfield(data,{'tau','step','fullRateNorm','coreRateNorm'}))
    count = numel(data.tau);
    if any([numel(data.step),numel(data.fullRateNorm),numel(data.coreRateNorm)]~=count)
        error('ipm:DatasetCatalog','Dataset scalar observation vectors are misaligned.');
    end
    for k = 1:count
        if key_in_history(data.tau(k),data.step(k),log) && ...
                ~beyond_state(data.tau(k),data.step(k),state)
            if isfinite(data.fullRateNorm(k)) && ...
                    (isfinite(data.coreRateNorm(k)) || isnan(data.coreRateNorm(k)))
                stream.observations = merge_observation(stream.observations, ...
                    data.tau(k),data.step(k),data.fullRateNorm(k),data.coreRateNorm(k));
            end
        elseif ~inherited
            error('ipm:DatasetBranch','Existing catalogue observations do not match the restart history.');
        end
    end
end
if isfield(data.meta,'replayChecks') && isstruct(data.meta.replayChecks)
    for k = 1:numel(data.meta.replayChecks)
        check = data.meta.replayChecks(k);
        if all(isfield(check,{'tau','step','kind','fields'})) && ...
                key_in_history(check.tau,check.step,log) && ...
                ~beyond_state(check.tau,check.step,state)
            stream.replayChecks = merge_replay_check(stream.replayChecks,check);
        end
    end
end
if isfield(data,'analysis') && (~inherited || ...
        (~any(data.tau>state.scale.canonicalTime+clock_tolerance(state.scale.canonicalTime)) && ...
         ~any(data.step>state.step)))
    stream.analysis = data.analysis;
end
if isfield(data,'references') && isstruct(data.references) && isscalar(data.references)
    stream.references = data.references;
end
if isfield(data.meta,'seedProvenance') && isstruct(data.meta.seedProvenance) && ...
        isscalar(data.meta.seedProvenance)
    stream.seedProvenance = data.meta.seedProvenance;
end
stream = import_extras(stream,data);
end

function stream = import_extras(stream,data)
for name = {'wx90','wy90','WX90'}
    key = name{1};
    if isfield(data,key)
        values = data.(key);
        if ~isnumeric(values) || ~isreal(values) || numel(values)~=numel(data.tau)
            error('ipm:DatasetExtras','Archive measurement %s is not aligned with scalar history.',key);
        end
        stream.seriesExtras.values.(key) = values(:);
    end
end
if ~isempty(fieldnames(stream.seriesExtras.values))
    stream.seriesExtras.tau = data.tau(:);
    stream.seriesExtras.step = data.step(:);
end
if isfield(data,'diagnostics') && isstruct(data.diagnostics)
    names = setdiff(fieldnames(data.diagnostics), ...
        {'common','gauge','mesh','anisotropic','wallCore'});
    for k = 1:numel(names)
        if strcmp(names{k},'monitor')
            stream.diagnosticExtras.archivedMonitor=struct('tau',data.tau, ...
                'step',data.step,'values',data.diagnostics.monitor, ...
                'sampling','fixed archive history, not extended by native output');
        else
            stream.diagnosticExtras.(names{k}) = data.diagnostics.(names{k});
        end
    end
end
end

function stream = recover_journal(stream,state,log)
files = dir(fullfile(stream.frameFolder,'frame_*.mat'));
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);
    loaded = load(file,'journal');
    if ~isfield(loaded,'journal') || ~isstruct(loaded.journal) || ...
            ~isscalar(loaded.journal)
        error('ipm:DatasetJournal','Invalid immutable frame journal: %s.',file);
    end
    journal = loaded.journal;
    required = {'caseId','fingerprint','tau','t','step','historyIndex', ...
        'nx','ny','reason','fullRateNorm','coreRateNorm'};
    if ~all(isfield(journal,required))
        error('ipm:DatasetJournal','Incomplete frame journal: %s.',file);
    end
    if ~strcmp(char(journal.caseId),stream.caseId) || ...
            ~strcmp(char(journal.fingerprint),stream.fingerprint)
        error('ipm:DatasetIdentity','Frame journal belongs to another case or numerical fingerprint: %s.',file);
    end
    if beyond_state(journal.tau,journal.step,state)
        error('ipm:DatasetBranch', ...
            'Local frame journal %s is ahead of the restart state; use a new output directory.',file);
    end
    index = history_index(journal.tau,journal.step,log);
    descriptor = struct('tau',journal.tau,'t',journal.t,'step',journal.step, ...
        'historyIndex',index,'nx',journal.nx,'ny',journal.ny, ...
        'file',absolute_path(file),'sha256',file_hash(file),'reason',char(journal.reason));
    present = find([stream.frames.tau]'==journal.tau & ...
        [stream.frames.step]'==journal.step,1);
    if isempty(present)
        recovered = read_frame(descriptor);
        if ~isequaln(recovered.fullRateNorm,journal.fullRateNorm) || ...
                ~isequaln(recovered.coreRateNorm,journal.coreRateNorm) || ...
                (isfield(journal,'remeshCount') && ...
                 recovered.remeshCount~=journal.remeshCount)
            error('ipm:DatasetJournal','Frame payload and journal disagree: %s.',file);
        end
    end
    [stream.frames,check] = merge_descriptor(stream.frames,descriptor);
    stream.replayChecks = merge_replay_check(stream.replayChecks,check);
    stream.observations = merge_observation(stream.observations, ...
        journal.tau,journal.step,journal.fullRateNorm,journal.coreRateNorm);
end
end

function [descriptor,journal,check] = write_frame(stream,frame,historyIndex,reason)
check = [];
file = fullfile(stream.frameFolder,sprintf('frame_step%012d_tau%s.mat', ...
    frame.step,num2hex(double(frame.tau))));
journal = struct('schemaVersion',1,'caseId',stream.caseId, ...
    'fingerprint',stream.fingerprint,'tau',frame.tau,'t',frame.t, ...
    'step',frame.step,'historyIndex',historyIndex, ...
    'nx',numel(frame.x),'ny',numel(frame.y),'reason',reason, ...
    'remeshCount',frame.remeshCount,'fullRateNorm',frame.fullRateNorm, ...
    'coreRateNorm',frame.coreRateNorm);
if isfile(file)
    existing = load(file,'frame','journal');
    if ~all(isfield(existing,{'frame','journal'})) || ...
            ~strcmp(char(existing.journal.caseId),stream.caseId) || ...
            ~strcmp(char(existing.journal.fingerprint),stream.fingerprint)
        error('ipm:DatasetFrameConflict','Refusing to replace immutable frame %s.',file);
    end
    check = compare_frames(existing.frame,frame);
    journal = existing.journal;
else
    temporaryFile = [tempname(stream.frameFolder),'.mat'];
    cleanup = onCleanup(@()delete_temporary(temporaryFile));
    save(temporaryFile,'frame','journal','-v7.3');
    if isfile(file)
        error('ipm:DatasetConcurrentWriter','Another writer created immutable frame %s.',file);
    end
    [moved,message] = movefile(temporaryFile,file);
    if ~moved
        error('ipm:DatasetFrameMove','Could not install %s: %s',file,message);
    end
end
descriptor = struct('tau',frame.tau,'t',frame.t,'step',frame.step, ...
    'historyIndex',historyIndex,'nx',numel(frame.x),'ny',numel(frame.y), ...
    'file',absolute_path(file),'sha256',file_hash(file),'reason',char(journal.reason));
end

function save_catalog(file,data)
temporaryFile = [tempname(fileparts(file)),'.mat'];
cleanup = onCleanup(@()delete_temporary(temporaryFile));
save(temporaryFile,'data','-v7.3');
[moved,message] = movefile(temporaryFile,file,'f');
if ~moved, error('ipm:DatasetCatalogMove','Could not publish %s: %s',file,message); end
end

function validate_cache(state)
if ~isfield(state,'rhsCache') || ~isstruct(state.rhsCache)
    error('ipm:DatasetCache','Accepted-state RHS cache is missing.');
end
cache = state.rhsCache;
required = {'rho','scale','x','y','remeshCount','rhoRate'};
if ~all(isfield(cache,required)) || ...
        ~isequaln(cache.rho,state.rho) || ~isequaln(cache.scale,state.scale) || ...
        ~isequal(cache.x,state.ops.x) || ~isequal(cache.y,state.ops.y) || ...
        cache.remeshCount~=state.ops.remeshCount || ...
        ~isequal(size(cache.rhoRate),size(state.rho)) || ...
        any(~isfinite(cache.rhoRate),'all')
    error('ipm:DatasetCache','RHS cache does not match the accepted density, scale, or native mesh.');
end
end

function [tau,step,index] = accepted_key(state,log)
tau = state.scale.canonicalTime;
step = state.step;
index = history_index(tau,step,log);
if index~=numel(log.history.common.canonicalTau)
    error('ipm:DatasetRecord','Dataset recording requires the terminal recorded accepted state.');
end
end

function assert_frame_key(frame,tau,step)
if ~isstruct(frame) || ~isscalar(frame) || ...
        ~all(isfield(frame,{'tau','t','step','x','y','rho','w','rho_xt', ...
        'remeshCount','fullRateNorm','coreRateNorm'})) || ...
        frame.tau~=tau || frame.step~=step || ...
        ~isequal(size(frame.rho),[numel(frame.y),numel(frame.x)]) || ...
        ~isequal(size(frame.w),size(frame.rho)) || ...
        ~isequal(size(frame.rho_xt),size(frame.rho)) || ...
        any(~isfinite(frame.rho),'all') || any(~isfinite(frame.rho_xt),'all')
    error('ipm:DatasetFrame','dataFrame returned an invalid accepted native frame.');
end
end

function [descriptors,check] = merge_descriptor(descriptors,next)
check = [];
same = find([descriptors.tau]'==next.tau & [descriptors.step]'==next.step,1);
if isempty(same)
    descriptors(end+1,1) = next;
elseif ~strcmp(descriptors(same).sha256,next.sha256)
    previousFrame = read_frame(descriptors(same));
    nextFrame = read_frame(next);
    check = compare_frames(previousFrame,nextFrame);
end
end

function check = compare_frames(saved,replayed)
% A same-grid checkpoint rebuilds the accepted RHS rather than storing the
% transient cache. Only these four RHS-derived values may differ by the
% measured reconstruction roundoff. Every other payload field is exact.
if ~isstruct(saved) || ~isscalar(saved) || ~isstruct(replayed) || ...
        ~isscalar(replayed) || ~isempty(setxor(fieldnames(saved),fieldnames(replayed)))
    error('ipm:DatasetFrameConflict','Immutable frame contracts differ.');
end
allowed = {'rho_t','rho_xt','fullRateNorm','coreRateNorm'};
differences = struct();
names = fieldnames(saved);
for k = 1:numel(names)
    name = names{k};
    a = saved.(name); b = replayed.(name);
    if isequaln(a,b), continue; end
    if ~ismember(name,allowed) || ~isnumeric(a) || ~isnumeric(b) || ...
            ~strcmp(class(a),class(b)) || ~isequal(size(a),size(b)) || ...
            any(~isfinite(a),'all') || any(~isfinite(b),'all')
        error('ipm:DatasetFrameConflict', ...
            'Immutable frame step %d tau %.17g differs in exact field %s.',saved.step,saved.tau,name);
    end
    maximumDifference = max(abs(a(:)-b(:)));
    tolerance = 128*eps(max([1;abs(a(:));abs(b(:))]));
    if maximumDifference>tolerance
        error('ipm:DatasetFrameConflict', ...
            ['Immutable frame step %d tau %.17g differs in %s by %.17g ' ...
            '(RHS reconstruction roundoff limit %.17g).'], ...
            saved.step,saved.tau,name,maximumDifference,tolerance);
    end
    differences.(name) = struct('maxDifference',maximumDifference,'tolerance',tolerance);
end
check = [];
if ~isempty(fieldnames(differences))
    check = struct('step',saved.step,'tau',saved.tau, ...
        'kind','derived_rhs_roundoff_only','fields',differences);
end
end

function checks = merge_replay_check(checks,check)
if isempty(check), return; end
if ~strcmp(char(check.kind),'derived_rhs_roundoff_only')
    error('ipm:DatasetReplayAudit','Unknown immutable frame reconstruction audit kind.');
end
index = find([checks.tau]'==check.tau & [checks.step]'==check.step,1);
if isempty(index)
    checks(end+1,1) = check;
else
    names = fieldnames(check.fields);
    for k = 1:numel(names)
        name = names{k};
        if isfield(checks(index).fields,name)
            previous = checks(index).fields.(name);
            value = check.fields.(name);
            value.maxDifference = max(previous.maxDifference,value.maxDifference);
            value.tolerance = max(previous.tolerance,value.tolerance);
            checks(index).fields.(name) = value;
        else
            checks(index).fields.(name) = check.fields.(name);
        end
    end
end
end

function frame = read_frame(descriptor)
if ~strcmp(file_hash(descriptor.file),descriptor.sha256)
    error('ipm:DatasetIntegrity','Frame SHA256 mismatch: %s.',descriptor.file);
end
loaded = load(descriptor.file,'frame');
if ~isfield(loaded,'frame') || ~isstruct(loaded.frame) || ~isscalar(loaded.frame)
    error('ipm:DatasetFrame','Frame variable is missing from %s.',descriptor.file);
end
frame = loaded.frame;
assert_frame_key(frame,descriptor.tau,descriptor.step);
if numel(frame.x)~=descriptor.nx || numel(frame.y)~=descriptor.ny || frame.t~=descriptor.t
    error('ipm:DatasetFrame','Frame payload differs from its descriptor: %s.',descriptor.file);
end
end

function descriptor = normalize_descriptor(input,folder)
required = fieldnames(empty_descriptor());
if ~all(isfield(input,required))
    error('ipm:DatasetCatalog','Frame descriptor is incomplete.');
end
descriptor = empty_descriptor();
for k = 1:numel(required), descriptor.(required{k}) = input.(required{k}); end
descriptor.file = char(descriptor.file);
if ~is_absolute(descriptor.file), descriptor.file = fullfile(folder,descriptor.file); end
descriptor.file = absolute_path(descriptor.file);
descriptor.sha256 = char(descriptor.sha256);
descriptor.reason = char(descriptor.reason);
if ~isfile(descriptor.file)
    error('ipm:DatasetFrameMissing','Native frame file is missing: %s.',descriptor.file);
end
for name = {'tau','t','step','historyIndex','nx','ny'}
    validateattributes(descriptor.(name{1}),{'numeric'}, ...
        {'scalar','real','finite','nonnegative'},mfilename,name{1});
end
end

function observations = merge_observation(observations,tau,step,fullNorm,coreNorm)
validateattributes([tau,step,fullNorm],{'numeric'}, ...
    {'real','finite','nonnegative'},mfilename,'observation');
if ~(isnumeric(coreNorm) && isreal(coreNorm) && isscalar(coreNorm) && ...
        (isnan(coreNorm) || (isfinite(coreNorm) && coreNorm>=0)))
    error('ipm:DatasetObservation','Core residual must be nonnegative finite or NaN if the core is empty.');
end
index = find(observations.tau==tau & observations.step==step,1);
if isempty(index)
    observations.tau(end+1,1) = tau;
    observations.step(end+1,1) = step;
    observations.fullRateNorm(end+1,1) = fullNorm;
    observations.coreRateNorm(end+1,1) = coreNorm;
    [~,order] = sortrows([observations.tau,observations.step],[1,2]);
    names = fieldnames(observations);
    for k = 1:numel(names), observations.(names{k}) = observations.(names{k})(order); end
else
    actual = [fullNorm,coreNorm];
    previous = [observations.fullRateNorm(index),observations.coreRateNorm(index)];
    finite = isfinite(actual) & isfinite(previous);
    scale = max([ones(size(actual));abs(actual);abs(previous)],[],1);
    if any(xor(isnan(actual),isnan(previous))) || ...
            any(abs(actual(finite)-previous(finite))>128*eps(scale(finite)))
        error('ipm:DatasetObservationConflict', ...
            'Accepted residual observations changed at step %d tau %.17g.',step,tau);
    end
end
end

function assert_history_prefix(data,log,state,inherited,file)
if ~all(isfield(data,{'tau','step'}))
    error('ipm:DatasetCatalog','Dataset %s lacks accepted-history keys.',file);
end
previous = struct();
if isfield(data,'diagnostics') && isstruct(data.diagnostics) && ...
        isfield(data.diagnostics,'common')
    previous = data.diagnostics.common;
elseif isfield(data,'history') && isstruct(data.history) && isfield(data.history,'common')
    previous = data.history.common;
end
current = log.history.common;
mapped = { ...
    'tau','canonicalTau'; 't','physicalTime'; 'step','acceptedStep'; ...
    'cl','canonicalCL'; 'co','canonicalCOmega'; 'cr','canonicalCR'; ...
    'Cx','C_x'; 'Cy','C_y'; 'Comega','C_omega'; 'Xshift','X_shift'; ...
    'M','physicalRhoXInf'};
for k = 1:size(mapped,1)
    if isfield(data,mapped{k,1}), previous.(mapped{k,2}) = data.(mapped{k,1}); end
end
names = intersect(fieldnames(previous),fieldnames(current));
count = numel(data.tau);
if numel(data.step)~=count
    error('ipm:DatasetCatalog','Dataset %s has misaligned accepted-history keys.',file);
end
for k = 1:count
    tau = data.tau(k); step = data.step(k);
    if inherited && beyond_state(tau,step,state), continue; end
    row = history_index(tau,step,log);
    for j = 1:numel(names)
        name = names{j};
        if numel(previous.(name))~=count || ...
                numel(current.(name))~=numel(current.canonicalTau)
            error('ipm:DatasetCatalog','Dataset accepted history has misaligned field %s.',name);
        end
        if ~isequaln(previous.(name)(k),current.(name)(row))
            error('ipm:DatasetBranch', ...
                'Dataset %s differs from checkpoint history at step %d tau %.17g (%s).',file,step,tau,name);
        end
    end
end
end

function index = history_index(tau,step,log)
if ~key_in_history(tau,step,log)
    error('ipm:DatasetHistory','No accepted history record matches step %d tau %.17g.',step,tau);
end
common = log.history.common;
index = find(common.canonicalTau(:)==tau & common.acceptedStep(:)==step,1,'last');
end

function found = key_in_history(tau,step,log)
if ~isfield(log,'history') || ~isfield(log.history,'common') || ...
        ~all(isfield(log.history.common,{'canonicalTau','acceptedStep'}))
    error('ipm:DatasetHistory','Accepted history clocks and step identifiers are required.');
end
common = log.history.common;
found = any(common.canonicalTau(:)==tau & common.acceptedStep(:)==step);
end

function result = beyond_state(tau,step,state)
result = tau>state.scale.canonicalTime+clock_tolerance(state.scale.canonicalTime) || step>state.step;
end

function source = restart_catalog(restartInput,caseId,fingerprint,state)
source = '';
hint = '';
if isstruct(restartInput) && isscalar(restartInput)
    if isfield(restartInput,'payload') && isfield(restartInput.payload,'state') && ...
            isfield(restartInput.payload.state,'runMetadata') && ...
            isfield(restartInput.payload.state.runMetadata,'datasetFile')
        hint = char(restartInput.payload.state.runMetadata.datasetFile);
    end
    if isfield(restartInput,'storageFile')
        restartInput = restartInput.storageFile;
    else
        restartInput = '';
    end
end
candidates = cell(0,1);
if ~isempty(hint) && isfile(hint)
    candidates{end+1,1} = absolute_path(hint);
end
if (ischar(restartInput) || (isstring(restartInput) && isscalar(restartInput))) && ...
        strlength(string(restartInput))>0
    folder = fileparts(absolute_path(char(restartInput)));
    % Segmented runs place CPs at run/segment/checkpoints; include the run
    % catalogue above the segment without recursively scanning other runs.
    for k = 1:4
        files = dir(fullfile(folder,'*.mat'));
        for j = 1:numel(files)
            candidates{end+1,1} = fullfile(files(j).folder,files(j).name); %#ok<AGROW>
        end
        parent = fileparts(folder);
        if strcmp(parent,folder), break; end
        folder = parent;
    end
end
candidates = unique(candidates,'stable');
best = [-Inf,-Inf,-Inf,-Inf];
for k = 1:numel(candidates)
    candidate = candidates{k};
    try
        variables = whos('-file',candidate);
    catch
        continue; % An unrelated legacy/report MAT may still be writing.
    end
    index = find(strcmp({variables.name},'data'),1);
    % A checkpoint/result may be very large. Inspect the MAT directory first
    % and load only the data variable from actual dataset candidates.
    if isempty(index) || ~strcmp(variables(index).class,'struct') || ...
            prod(variables(index).size)~=1
        continue;
    end
    try
        loaded = load(candidate,'data');
    catch
        continue;
    end
    data = loaded.data;
    if ~isfield(data,'meta') || ~isstruct(data.meta) || ~isscalar(data.meta) || ...
            ~all(isfield(data.meta,{'caseId','storage','fingerprint'})) || ...
            ~strcmp(char(data.meta.caseId),caseId) || ...
            ~strcmp(char(data.meta.storage),'native_frame_journal')
        continue;
    end
    containsSource = false;
    if all(isfield(data,{'tau','step'})) && numel(data.tau)==numel(data.step)
        containsSource = any(data.tau(:)==state.scale.canonicalTime & data.step(:)==state.step);
    end
    eligibleFrames = 0;
    if isfield(data,'frames') && isstruct(data.frames) && ...
            all(isfield(data.frames,{'tau','step'}))
        eligibleFrames = nnz([data.frames.tau]<=state.scale.canonicalTime+ ...
            clock_tolerance(state.scale.canonicalTime) & [data.frames.step]<=state.step);
    end
    info = dir(candidate);
    score = [strcmp(char(data.meta.fingerprint),fingerprint), ...
        double(containsSource),eligibleFrames,info.datenum];
    firstDifference = find(score~=best,1);
    if ~isempty(firstDifference) && score(firstDifference)>best(firstDifference)
        best = score;
        source = candidate;
    end
    clear loaded data
end
end

function fingerprint = numerical_fingerprint(config,numericsSha256)
if isfield(config,'output'), config = rmfield(config,'output'); end
terminal = {'finalTime','maxSteps','physicalFinalTime'};
if isfield(config,'time')
    present = terminal(isfield(config.time,terminal));
    if ~isempty(present), config.time = rmfield(config.time,present); end
end
value = struct('config',config,'numericsSha256',numericsSha256);
encoded = jsonencode(canonical_value(value),'ConvertInfAndNaN',false);
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(typecast(unicode2native(encoded,'UTF-8'),'int8'));
fingerprint = digest_hex(digest);
end

function hash = numerics_hash()
% Only numerical sources participate. Output code and solve orchestration
% may evolve without changing the numerical identity of a checkpoint.
flowFile = which('ipm.evolve.flow');
if isempty(flowFile), error('ipm:DatasetNumerics','Cannot identify the active numerical runtime.'); end
packageRoot = fileparts(fileparts(flowFile));
parts = {'+evolve','+field','+mesh','+remesh','+config','+diagnostics'};
relative = cell(0,1);
for k = 1:numel(parts)
    files = dir(fullfile(packageRoot,parts{k},'**','*.m'));
    for j = 1:numel(files)
        absolute = fullfile(files(j).folder,files(j).name);
        relative{end+1,1} = strrep(absolute(numel(packageRoot)+2:end),filesep,'/'); %#ok<AGROW>
    end
end
relative = sort(relative);
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
for k = 1:numel(relative)
    digest.update(typecast(unicode2native(relative{k},'UTF-8'),'int8'));
    digest.update(int8(0));
    contentHash = file_hash(fullfile(packageRoot,strrep(relative{k},'/',filesep)));
    digest.update(typecast(uint8(contentHash),'int8'));
    digest.update(int8(0));
end
hash = digest_hex(digest);
end

function descriptors = public_descriptors(descriptors,folder)
prefix = [folder,filesep];
for k = 1:numel(descriptors)
    if startsWith(descriptors(k).file,prefix)
        descriptors(k).file = descriptors(k).file(numel(prefix)+1:end);
    end
end
end

function output = canonical_value(input)
if isa(input,'function_handle')
    output = struct('functionHandle',func2str(input));
elseif isstruct(input)
    if ~isscalar(input)
        output = arrayfun(@canonical_value,input,'UniformOutput',false);
        return;
    end
    output = struct(); names = sort(fieldnames(input));
    for k = 1:numel(names), output.(names{k}) = canonical_value(input.(names{k})); end
elseif iscell(input)
    output = cell(size(input));
    for k = 1:numel(input), output{k} = canonical_value(input{k}); end
else
    output = input;
end
end

function hash = file_hash(file)
[fid,message] = fopen(file,'rb');
if fid<0, error('ipm:DatasetRead','Cannot read %s: %s',file,message); end
cleanup = onCleanup(@()fclose(fid));
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
while true
    bytes = fread(fid,1024*1024,'*uint8');
    if isempty(bytes), break; end
    digest.update(typecast(bytes,'int8'));
end
[message,errorNumber] = ferror(fid);
if errorNumber~=0 && ~feof(fid)
    error('ipm:DatasetRead','Cannot hash %s: %s',file,message);
end
hash = digest_hex(digest);
end

function hash = digest_hex(digest)
bytes = typecast(digest.digest(),'uint8');
hash = lower(reshape(dec2hex(bytes,2).',1,[]));
end

function path = absolute_path(path)
if ~is_absolute(path), path = fullfile(pwd,path); end
file = javaObject('java.io.File',path);
path = char(file.getCanonicalPath());
end

function yes = is_absolute(path)
yes = startsWith(path,filesep) || ~isempty(regexp(path,'^[A-Za-z]:[\\/]','once'));
end

function next = next_frame_tau(tau,every)
next = (floor((tau+clock_tolerance(tau))/every)+1)*every;
if next<=tau, next = tau+every; end
end

function tolerance = clock_tolerance(tau)
tolerance = 64*eps(max(1,abs(tau)));
end

function descriptors = empty_descriptors()
descriptors = repmat(empty_descriptor(),0,1);
end

function descriptor = empty_descriptor()
descriptor = struct('tau',0,'t',0,'step',0,'historyIndex',0,'nx',0,'ny',0, ...
    'file','','sha256','','reason','');
end

function observations = empty_observations()
observations = struct('tau',zeros(0,1),'step',zeros(0,1), ...
    'fullRateNorm',zeros(0,1),'coreRateNorm',zeros(0,1));
end

function checks = empty_replay_checks()
checks = repmat(struct('step',0,'tau',0,'kind','derived_rhs_roundoff_only', ...
    'fields',struct()),0,1);
end

function [lock,cleanup] = acquire_writer_lock(file)
% Keep the lock inode permanently. Removing a lock file can let two writers
% lock different inodes with the same path. The OS releases this advisory
% lock when the process exits, including a crash.
same_process_writer(file);
lockFile = [];
channel = [];
try
    lockFile = javaObject('java.io.RandomAccessFile',[file,'.lock'],'rw');
    channel = lockFile.getChannel();
catch exception
    release_writer_lock([],channel,lockFile);
    error('ipm:DatasetLock','Cannot open dataset writer lock for %s: %s',file,exception.message);
end
try
    lock = channel.tryLock();
catch exception
    % Java also throws OverlappingFileLockException for a second writer in
    % this MATLAB process rather than returning a null lock.
    release_writer_lock([],channel,lockFile);
    error('ipm:DatasetBusy','Dataset already has a writer or its exclusive lock is unavailable: %s (%s)', ...
        file,exception.message);
end
if isempty(lock)
    release_writer_lock([],channel,lockFile);
    error('ipm:DatasetBusy','Dataset already has an active writer: %s.',file);
end
cleanup = onCleanup(@()release_writer_lock(lock,channel,lockFile));
same_process_writer(file,lock);
end

function same_process_writer(file,lock)
% Avoid opening/closing a second descriptor to a file already locked by
% this writer. On some POSIX systems closing any descriptor can release
% another channel's process-wide advisory lock on the same inode.
persistent references
if isempty(references)
    references = containers.Map('KeyType','char','ValueType','any');
end
if nargin>1
    references(file) = javaObject('java.lang.ref.WeakReference',lock);
    return;
end
if isKey(references,file)
    reference = references(file);
    previous = reference.get();
    if ~isempty(previous) && previous.isValid()
        error('ipm:DatasetBusy','Dataset already has an active writer in this MATLAB process: %s.',file);
    end
    remove(references,file);
end
end

function assert_writer_lock(stream)
if ~isfield(stream,'writerLock') || isempty(stream.writerLock) || ...
        ~stream.writerLock.isValid()
    error('ipm:DatasetBusy','Dataset writer no longer holds its exclusive lock.');
end
end

function release_writer_lock(lock,channel,lockFile)
% Cleanup must release every resource even if a preceding close fails.
if ~isempty(lock)
    try
        if lock.isValid(), lock.release(); end
    catch
    end
end
if ~isempty(channel)
    try
        channel.close();
    catch
    end
end
if ~isempty(lockFile)
    try
        lockFile.close();
    catch
    end
end
end

function delete_temporary(file)
if isfile(file), delete(file); end
end
