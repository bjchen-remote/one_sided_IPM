function value=identity(action,root)
%IDENTITY Portable content identity; machine locations are never identities.
if nargin<2,root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));end
switch action
    case 'root',value=root;
    case 'utc'
        value=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
    case 'source'
        listing=dir(fullfile(root,'src','**','*'));paths={};
        for k=1:numel(listing)
            [~,~,extension]=fileparts(listing(k).name);
            if ~listing(k).isdir&&ismember(extension,{'.m','.c','.cpp',['.',mexext]})
                paths{end+1}=fullfile(listing(k).folder,listing(k).name); %#ok<AGROW>
            end
        end
        paths{end+1}=fullfile(root,'data','initial.json');paths=sort(unique(paths));
        rows=struct('path',{},'sha256',{});
        for k=1:numel(paths)
            relative=project_relative(root,paths{k});
            rows(end+1)=struct('path',relative,'sha256',ipm_cr.io.sha(paths{k})); %#ok<AGROW>
        end
        value=struct('format','ipm_cr_source_v2','files',rows);
    case 'authority'
        path=fullfile(root,'data','initial.json');meta=jsondecode(fileread(path));
        assert(strcmp(meta.format,'ipm_cr_initial_v1')&&isfield(meta,'file')&&isfield(meta,'sha256'), ...
            'ipm_cr:AuthorityManifest','Use a complete initial-data manifest.');
        relative=char(meta.file);parts=strsplit(relative,'/');
        assert(~isempty(relative)&&~startsWith(relative,'/')&& ...
            ~contains(relative,'\')&&~any(ismember(parts,{'','..','.'})), ...
            'ipm_cr:AuthorityPath','Initial input must be project-relative.');
        input=fullfile(root,parts{:});project_relative(root,input);
        actual=ipm_cr.io.sha(input);
        assert(strcmp(actual,meta.sha256),'ipm_cr:AuthoritySha','Initial input SHA differs from its manifest.');
        value=struct('format','ipm_cr_authority_v2','file',relative,'sha256',actual);
    otherwise,error('ipm_cr:IdentityAction','Unknown identity action %s.',action);
end
end

function relative=project_relative(root,path)
baseFile=java.io.File(root);inputFile=java.io.File(path);
base=char(baseFile.getCanonicalPath());actual=char(inputFile.getCanonicalPath());prefix=[base,filesep];
assert(startsWith(actual,prefix),'ipm_cr:SourcePath','An input resolves outside the project.');
relative=strrep(actual(numel(prefix)+1:end),filesep,'/');
end
