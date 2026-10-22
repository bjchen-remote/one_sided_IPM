function report=build_ipm_cr_native(outputDir,backend)
%BUILD_IPM_CR_NATIVE Build the OPTIONAL smooth-WENO/linear-optimal prototype.
% Does not select a transport backend or change default numerical options.
% One binary supports optional reconstructionMode=0 (linear optimal) and
% reconstructionMode=1 (smooth WENO-Z, also the six-argument default).
% The linear mode combines face reconstructions into direct fifth-upwind
% coefficients while retaining the same common-reference ghost differences.
% Uses the interleaved-complex MEX API and strict floating-point flags; the
% common-reference differences in this kernel must not be reassociated.
% backend='auto' prefers the installed macOS CLT directly, bypassing only
% mex's compiler discovery; backend='mex' retains the original mex route.
% backend='clt' explicitly selects CLT. No license or global setup is changed.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'time_evolution'));
packagePath=fullfile(root,'test','+ipm_cr','+explicit');
source=fullfile(packagePath,'wenoPairSmoothMex.cpp');
assert(isfile(source),'ipm_cr:NativeBuildSource','Missing native source.');
if nargin<1||isempty(outputDir)
    results=fullfile(root,'test','results');
    if ~isfolder(results),mkdir(results);end
    [~,token]=fileparts(tempname(results));
    outputDir=fullfile(results,['native_build_',token]);
end
outputDir=char(outputDir);
assert(startsWith(outputDir,filesep),'ipm_cr:NativeBuildOutput','Use an absolute output directory.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'native_build.json')), ...
    'ipm_cr:NativeBuildPreserve','Choose a fresh build report directory.');
if nargin<2||isempty(backend),backend='auto';end
backend=lower(char(backend));
assert(any(strcmp(backend,{'auto','mex','clt'})),'ipm_cr:NativeBuildBackend', ...
    'backend must be auto, mex or clt.');
selectedBackend=backend;
clt='/Library/Developer/CommandLineTools/usr/bin/clang++';
if strcmp(backend,'auto')
    if ismac&&isfile(clt),selectedBackend='clt';else,selectedBackend='mex';end
end
args={'-R2018a','-outdir',packagePath,'-output','wenoPairSmoothMex',source};
if ispc
    args{end+1}='COMPFLAGS=$COMPFLAGS /fp:strict';
else
    args{end+1}='CXXFLAGS=$CXXFLAGS -fno-fast-math -ffp-contract=off';
end
report=struct('kind','optional_smooth_weno_native_build_v3','status','running', ...
    'source',source,'sourceSha256',ipm_c.output.sha256(source), ...
    'mexArguments',{args},'matlabVersion',version,'platform',computer, ...
    'requestedBackend',backend,'selectedBackend',selectedBackend, ...
    'mexExtension',mexext,'fastMathEnabled',false, ...
    'supportedReconstructionModes',[0,1],'defaultReconstructionMode',1, ...
    'linearDerivativeImplementation','direct common-reference fifth-upwind coefficients on shared seven-point support', ...
    'transportDefaultChanged',false,'sourceUnchangedDuringBuild',false);
timer=tic;
try
    clear ipm_cr.explicit.wenoPairSmoothMex
    if strcmp(selectedBackend,'clt')
        [report,buildStatus]=build_clt(report,clt,outputDir,source);
        assert(buildStatus==0,'ipm_cr:NativeBuildCLT', ...
            'Direct CLT build failed (status %d); see %s.',buildStatus,report.buildLog);
    else
        compiler=mex.getCompilerConfigurations('C++','Selected');
        assert(~isempty(compiler),'ipm_cr:NativeBuildCompiler','Select a C++ MEX compiler first.');
        report.compiler=struct('name',compiler.Name,'manufacturer',compiler.Manufacturer, ...
            'version',compiler.Version,'location',compiler.Location);
        buildLog=evalc('mex(args{:});');
        report.buildLog=fullfile(outputDir,'mex.log');
        write_text(report.buildLog,buildLog);
    end
    report.sourceUnchangedDuringBuild=strcmp(report.sourceSha256,ipm_c.output.sha256(source));
    assert(report.sourceUnchangedDuringBuild,'ipm_cr:NativeBuildSourceChanged', ...
        'Native source changed during compilation; this build is not source-pinned.');
    binary=fullfile(packagePath,['wenoPairSmoothMex.',mexext]);
    if strcmp(selectedBackend,'clt')
        assert(strcmp(report.versionSourceSha256,ipm_c.output.sha256(report.versionSource))&& ...
            strcmp(report.exportMapSha256,ipm_c.output.sha256(report.exportMap)), ...
            'ipm_cr:NativeBuildVersionChanged','MATLAB API-version source changed during compilation.');
        [ok,message]=movefile(report.stagedBinary,binary,'f');
        assert(ok,'ipm_cr:NativeBuildPublish','%s',message);
    end
    assert(isfile(binary),'ipm_cr:NativeBuildMissing','MEX did not produce the expected binary.');
    report.binary=binary;report.binarySha256=ipm_c.output.sha256(binary);
    report.status='completed';report.elapsedSeconds=toc(timer);
    write_report(outputDir,report);
catch problem
    report.status='failed';report.elapsedSeconds=toc(timer);
    report.error=struct('identifier',problem.identifier,'message',problem.message,'stack',problem.stack);
    write_report(outputDir,report);rethrow(problem);
end
end

function [report,status]=build_clt(report,compiler,folder,source)
assert(ismac&&isfile(compiler),'ipm_cr:NativeBuildCLTUnavailable', ...
    'Direct CLT builds require the installed macOS Command Line Tools clang++.');
sdk='/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk';
assert(isfolder(sdk),'ipm_cr:NativeBuildSDK','Missing installed CLT macOS SDK.');
releaseRoot=matlabroot;includePath=fullfile(releaseRoot,'extern','include');
versionSource=fullfile(releaseRoot,'extern','version','cpp_mexapi_version.cpp');
platform=computer('arch');
switch platform
    case 'maca64',architecture='arm64';deployment='13.3';
    case 'maci64',architecture='x86_64';deployment='10.15';
    otherwise,error('ipm_cr:NativeBuildArchitecture','Unsupported macOS platform %s.',platform);
end
libraryPath=fullfile(releaseRoot,'bin',platform);
exportMap=fullfile(releaseRoot,'extern','lib',platform,'c_exportsmexfileversion.map');
assert(isfile(versionSource)&&isfile(exportMap)&&isfile(fullfile(includePath,'mex.h')), ...
    'ipm_cr:NativeBuildMatlabSDK','Missing MATLAB headers, API-version source or export map.');
for library={'libmx.dylib','libmex.dylib','libmat.dylib'}
    assert(isfile(fullfile(libraryPath,library{1})),'ipm_cr:NativeBuildLibrary', ...
        'Missing MATLAB library %s.',library{1});
end
report.matlabRoot=releaseRoot;report.sdk=sdk;report.architecture=architecture;
report.minimumMacOSVersion=deployment;report.compilerExecutable=compiler;
[compilerStatus,report.compilerVersion]=system([shell_quote(compiler),' --version 2>&1']);
assert(compilerStatus==0,'ipm_cr:NativeBuildCompiler','Cannot execute installed CLT clang++.');
report.versionSource=versionSource;report.versionSourceSha256=ipm_c.output.sha256(versionSource);
report.exportMap=exportMap;report.exportMapSha256=ipm_c.output.sha256(exportMap);
report.apiContract=struct('defaultRelease','R2018a','interleavedComplex',true, ...
    'macroDefinition','matrix.h derives MX_HAS_INTERLEAVED_COMPLEX from the R2018a target API; it is not redundantly overridden.', ...
    'versionEntryPoint','mexfilerequiredapiversion from the installed cpp_mexapi_version.cpp');
report.stagedBinary=fullfile(folder,['wenoPairSmoothMex.',mexext]);
commandArgs={compiler,'-arch',architecture,['-mmacosx-version-min=',deployment], ...
    '-isysroot',sdk,'-std=c++14','-stdlib=libc++','-O2','-DNDEBUG','-fexceptions', ...
    '-fno-common','-fwrapv','-fno-fast-math','-ffp-contract=off', ...
    '-DMATLAB_MEX_FILE','-DMATLAB_DEFAULT_RELEASE=R2018a', ...
    '-DMATLAB_MEXCMD_RELEASE=R2018a','-DUSE_MEX_CMD','-DMW_NEEDS_VERSION_H', ...
    ['-I',includePath],source,versionSource,'-bundle','-Wl,-twolevel_namespace', ...
    ['-Wl,-exported_symbols_list,',exportMap],['-L',libraryPath], ...
    '-lmx','-lmex','-lmat',['-Wl,-rpath,',libraryPath],'-framework','Cocoa', ...
    '-o',report.stagedBinary};
report.directCltArguments=commandArgs;
report.buildCommand=strjoin(cellfun(@shell_quote,commandArgs,'UniformOutput',false),' ');
write_text(fullfile(folder,'build_command.txt'),[report.buildCommand,newline]);
report.buildLog=fullfile(folder,'clang.log');
[status,buildLog]=system([report.buildCommand,' 2>&1']);
write_text(report.buildLog,buildLog);report.compilerExitStatus=status;
if status==0
    [symbolStatus,symbols]=system(['/usr/bin/nm -gU ',shell_quote(report.stagedBinary),' 2>&1']);
    report.exportedSymbols=symbols;report.symbolInspectionExitStatus=symbolStatus;
    report.requiredEntryPointsVerified=symbolStatus==0&&contains(symbols,'_mexFunction')&& ...
        contains(symbols,'_mexfilerequiredapiversion');
    if ~report.requiredEntryPointsVerified,status=101;return;end
    [linkStatus,report.linkedLibraries]=system(['/usr/bin/otool -L ',shell_quote(report.stagedBinary),' 2>&1']);
    report.linkInspectionExitStatus=linkStatus;
    if linkStatus~=0,status=102;end
end
end

function result=shell_quote(value)
single=char(39);escapedSingle=[single,'"',single,'"',single];
result=[single,strrep(char(value),single,escapedSingle),single];
end

function write_text(path,value)
fid=fopen(path,'w');assert(fid>=0,'ipm_cr:NativeBuildLog','Cannot write %s.',path);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,value,'char');
end

function write_report(folder,report)
save(fullfile(folder,'native_build.mat'),'report');
fid=fopen(fullfile(folder,'native_build.json'),'w');
assert(fid>=0,'ipm_cr:NativeBuildReport','Cannot save MEX build provenance.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');
end
