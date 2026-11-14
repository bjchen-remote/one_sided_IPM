function value=identity(action,root,initial)
%IDENTITY Portable content identity; machine locations are never identities.
if nargin<2||isempty(root),root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));end
if nargin<3,initial='data/initial.json';end
switch action
    case 'root',value=root;
    case 'sha',value=sha(root);
    case 'json',write_json(root,initial);value=[];
    case 'utc'
        value=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
    case 'source'
        root=char(java.io.File(root).getCanonicalPath());
        listing=dir(fullfile(root,'src','**','*'));paths={};
        for k=1:numel(listing)
            [~,~,extension]=fileparts(listing(k).name);
            if ~listing(k).isdir&&ismember(extension,{'.m','.c','.cpp',['.',mexext]})
                paths{end+1}=fullfile(listing(k).folder,listing(k).name); %#ok<AGROW>
            end
        end
        paths{end+1}=fullfile(root,'data','initial.json');
        paths{end+1}=fullfile(root,'main.m');
        paths{end+1}=fullfile(root,'data','mesh.json');paths=sort(unique(paths));
        rows=struct('path',{},'sha256',{});
        for k=1:numel(paths)
            relative=project_relative(root,paths{k});
            rows(end+1)=struct('path',relative,'sha256',sha(paths{k})); %#ok<AGROW>
        end
        value=struct('format','ipm_cr_source_v2','files',rows);
    case 'authority'
        [manifestFile,path]=authority_path(root,initial);meta=jsondecode(fileread(path));
        assert(isstruct(meta)&&isscalar(meta)&&all(isfield(meta,{'format','file','sha256'}))&& ...
            strcmp(meta.format,'ipm_cr_initial_v1'), ...
            'ipm_cr:AuthorityManifest','Use a complete initial-data manifest.');
        [relative,input]=authority_path(root,meta.file);
        actual=sha(input);
        assert(((ischar(meta.sha256)&&isrow(meta.sha256))||(isstring(meta.sha256)&&isscalar(meta.sha256)))&& ...
            ~isempty(regexp(char(meta.sha256),'^[0-9a-fA-F]{64}$','once')), ...
            'ipm_cr:AuthoritySha','The initial MAT SHA must contain 64 hexadecimal digits.');
        assert(strcmpi(actual,meta.sha256),'ipm_cr:AuthoritySha','Initial input SHA differs from its manifest.');
        value=struct('format','ipm_cr_authority_v2','file',relative,'sha256',actual, ...
            'manifestFile',manifestFile,'manifestSha256',sha(path));
    otherwise,error('ipm_cr:IdentityAction','Unknown identity action %s.',action);
end
end

function [relative,path]=authority_path(root,value)
assert((ischar(value)&&isrow(value))||(isstring(value)&&isscalar(value)), ...
    'ipm_cr:AuthorityPath','Use a project-relative initial input path.');
relative=char(value);parts=strsplit(relative,'/','CollapseDelimiters',false);
assert(~isempty(relative)&&~startsWith(relative,'/')&&~contains(relative,'\')&& ...
    ~contains(relative,':')&&~any(ismember(parts,{'','..','.'})), ...
    'ipm_cr:AuthorityPath','Initial inputs must use canonical project-relative paths.');
base=char(java.io.File(root).getCanonicalPath());
path=char(java.io.File(fullfile(root,parts{:})).getCanonicalPath());
assert(startsWith(path,[base,filesep]),'ipm_cr:AuthorityPath','An initial input resolves outside the project.');
end

function relative=project_relative(root,path)
baseFile=java.io.File(root);inputFile=java.io.File(path);
base=char(baseFile.getCanonicalPath());actual=char(inputFile.getCanonicalPath());prefix=[base,filesep];
assert(startsWith(actual,prefix),'ipm_cr:SourcePath','An input resolves outside the project.');
relative=strrep(actual(numel(prefix)+1:end),filesep,'/');
end

function value=sha(path)
%SHA256 Streaming file identity, without reading a large checkpoint at once.
d=java.security.MessageDigest.getInstance('SHA-256');
f=fopen(path,'rb');assert(f>=0,'ipm_c:HashFile','Validation failed.');cleanup=onCleanup(@()fclose(f)); %#ok<NASGU>
while true
    bytes=fread(f,1048576,'*uint8');if isempty(bytes),break;end;d.update(bytes);
end
value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2)',1,[]));
end

function write_json(path,value)
%JSON Replace a UTF-8 JSON record atomically in its containing directory.
temporary=[char(path),'.writing'];fid=fopen(temporary,'wb');
assert(fid>=0,'ipm_cr:JsonOpen','Cannot open JSON staging file.');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,unicode2native(jsonencode(value,'PrettyPrint',true),'UTF-8'),'uint8');clear cleanup
[ok,message]=movefile(temporary,path,'f');
assert(ok,'ipm_cr:JsonRename','%s',message);
end
