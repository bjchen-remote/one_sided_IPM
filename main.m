function [summary,state,ctx]=main(input,density)
%MAIN Edit the settings below and press Run in the MATLAB editor.
% [c,rho0]=main('config') returns settings without building or solving.
% main(c) uses c.initial; main(c,rho0) accepts a density expression @(x,y).
% Relative output folders here are relative to this main.m, on every OS.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'),fullfile(root,'tools'),'-begin');
assert(strcmp(char(java.io.File(ipm_cr.io.identity('root')).getCanonicalPath()), ...
    char(java.io.File(root).getCanonicalPath())), ...
    'ipm_cr:MainPath','Another ipm_cr package is shadowing this project.');

% ===== User settings: edit this block, then press Run =====
c=ipm_cr.config('case');
c.run.tauEnd=200;                 % absolute final rescaled time tau
c.run.wallSeconds=600;            % budget for this call, in seconds
c.run.outputDir=fullfile(root,'result');
c.run.resume=true;                % continue this folder's accepted state
c.domain.rmax=1e10;               % common transport and Poisson radius
c.initial='data/initial_small_source.json';    % used when rho0=[] or when calling main(c)

% Write any vectorized initial density here; no Gaussian parameter fit.
rho0=[];                         % original registered Gaussian datum
% rho0=@(x,y) (1+tanh(x)).*exp(-y.^2); % optional vectorized alternative
% Set rho0=[] to use c.initial unchanged.

% Nonuniform thin staircase. Use 'coarse', 'medium', or 'fine'.
c.method.geometry.preset='medium';
c.method.remesh=struct('enabled',true);
c.method.cfl=.4;                  % current-rate step formula; stages are observed
c.stop.physicalTime=[];           % optional original-physical-time threshold
c.save.data=true;                 % native fields, geometry, scales and SHA
c.save.history=true;
c.save.label='main';
% Change outputDir for a changed mesh/initial problem; old data are preserved.
% ===== End of user settings =====

if nargin==0
    input=c;density=rho0;
elseif ischar(input)||(isstring(input)&&isscalar(input))
    assert(nargin==1&&strcmp(input,'config'),'ipm_cr:MainAction', ...
        'Use main, main(''config''), or main(case[,rho0]).');
    summary=c;state=rho0;ctx=[];return
elseif nargin<2
    density=[];
end
assert(isstruct(input)&&isscalar(input)&&isfield(input,'format')&& ...
    strcmp(input.format,'ipm_cr_case_v1'), ...
    'ipm_cr:MainCase','Use a case from main(''config'') or ipm_cr.config(''case'').');
if ~isfield(input,'run'),input.run=struct();end
assert(isstruct(input.run)&&isscalar(input.run),'ipm_cr:MainCase','case.run must be a scalar structure.');
if ~isfield(input.run,'outputDir'),input.run.outputDir=fullfile(root,'result');end
if isfield(input,'run')&&isfield(input.run,'outputDir')
    folder=input.run.outputDir;
    assert((ischar(folder)&&isrow(folder))||(isstring(folder)&&isscalar(folder)), ...
        'ipm_cr:MainOutput','Output directory must be nonempty text.');
    if ~java.io.File(char(folder)).isAbsolute()
        input.run.outputDir=fullfile(root,char(folder));
    end
end
% Validate controls before building or generating any input files.
ipm_cr.config(input,'run');
samplingBinary=fullfile(root,'src','+ipm_cr','+mesh',['sampleMex.',mexext]);
if ~isfile(samplingBinary)
    fprintf('Building the local %s sample kernel...\n',mexext);
    build_ipm_cr_native();
    rehash;
end
selected=which('ipm_cr.mesh.sampleMex');
assert(~isempty(selected)&&strcmp(char(java.io.File(selected).getCanonicalPath()), ...
    char(java.io.File(samplingBinary).getCanonicalPath())), ...
    'ipm_cr:MainKernel','The current project platform kernel is not on the MATLAB path.');
if ~isempty(density)
    input.initial=ipm_cr.model.initial.materialize(input,density);
end
[summary,state,ctx]=ipm_cr.run(input);
if ~isempty(state)
    fprintf('Saved accepted tau=%.12g, step=%d, status=%s\n%s\n', ...
        state.scale.canonicalTime,state.step,summary.status,summary.outputDir);
else
    fprintf('status=%s\n%s\n',summary.status,summary.outputDir);
end
end

