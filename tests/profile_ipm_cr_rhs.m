function report=profile_ipm_cr_rhs(outputDir,repetitions)
%PROFILE_IPM_CR_RHS Measure the complete coupled RHS, excluding startup.
if nargin<2,repetitions=40;end
assert(~isfolder(outputDir),'ipm_cr:ProfilePreserve','Use a fresh output directory.');
mkdir(outputDir);
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
[state,g,wallOps,options,ctx]=ipm_cr.evolve.init();
guess=[state.flow.c_l,state.flow.c_r];
run_rhs();
profile clear;profile on;
timer=tic;
for j=1:repetitions
    flow=run_rhs();
end
elapsed=toc(timer);profile off;data=profile('info');
[~,order]=sort([data.FunctionTable.TotalTime],'descend');
rows=struct([]);
for j=order(1:min(35,numel(order)))
    f=data.FunctionTable(j);
    row=struct('name',f.FunctionName,'file',f.FileName, ...
        'calls',f.NumCalls,'totalSeconds',f.TotalTime);
    if isempty(rows),rows=row;else,rows(end+1)=row;end
end
report=struct('repetitions',repetitions,'profiledSeconds',elapsed, ...
    'profiledMeanRhsSeconds',elapsed/repetitions,'functions',rows, ...
    'gaugeResidual',max(abs(flow.gaugeResidual)), ...
    'poisson',flow.poisson,'matlabVersion',version);
save(fullfile(outputDir,'profile.mat'),'data','report','-v7');
fid=fopen(fullfile(outputDir,'profile.json'),'w');cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');
disp(struct2table(rows(1:min(15,numel(rows)))));
    function flow=run_rhs()
        if isfield(state,'residual')
            flow=ipm_cr.evolve.residual(state.residual,state.scale,g,wallOps,options,guess,true,ctx);
        else
            flow=ipm_cr.evolve.rhs(state.rho,g,wallOps,options,guess,[],true,ctx);
        end
    end
end
