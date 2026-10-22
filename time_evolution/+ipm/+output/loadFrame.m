function frame = loadFrame(data,index)
%IPM.OUTPUT.LOADFRAME Read exactly one immutable native frame, or packed frame.
validateattributes(index,{'numeric'},{'scalar','integer','>=',1,'<=',numel(data.frames)});
entry = data.frames(index);
if isfield(entry,'payload') && ~isempty(entry.payload)
    frame = entry.payload;
else
    file = char(entry.file);
    absolute = startsWith(file,filesep) || ~isempty(regexp(file,'^[A-Za-z]:[\\/]','once'));
    if ~absolute && isfield(data,'meta') && isfield(data.meta,'baseFolder')
        file = fullfile(data.meta.baseFolder,file);
    end
    if ~isfile(file)
        error('ipm:DatasetFrameMissing', ...
            'Native frame is missing: %s. Move the dataset with its frame directory or use packData.',file);
    end
    if isfield(entry,'sha256') && ~isempty(entry.sha256)
        assert(strcmp(hash_file(file),entry.sha256),'ipm:DatasetFrameIntegrity', ...
            'Native frame content changed: %s',file);
    end
    stored = load(file,'frame');
    assert(isfield(stored,'frame'),'ipm:DatasetFrameSchema','No frame structure in %s.',file);
    frame = stored.frame;
end
assert(frame.tau==entry.tau && frame.step==entry.step, ...
    'ipm:DatasetFrameIndex','Frame time/step does not match its index.');
assert(isequal(size(frame.rho_xt),[numel(frame.y),numel(frame.x)]), ...
    'ipm:DatasetFrameShape','Native grid and residual shape do not match.');
end

function hash = hash_file(file)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(file,'rb'); assert(fid>0,'Cannot read %s',file);
guard = onCleanup(@()fclose(fid));
while true
    bytes = fread(fid,1048576,'*uint8');
    if isempty(bytes), break; end
    md.update(bytes);
end
hash = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
