function report=test_release(outputDir)
%TEST_RELEASE Small native run, accepted-only stops, recovery and data export.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'),fullfile(root,'tools'));
if nargin<1,outputDir=fullfile(root,'result','test_release');end
if ~isfolder(outputDir),mkdir(outputDir);end
c=ipm_cr.config('case');c.run.outputDir=fullfile(outputDir,'resume');
c.run.tauEnd=.003;c.run.wallSeconds=120;c.run.verbose=false;c.run.maxSteps=1;
c.method=struct('geometry',struct('coreSpacing',.06,'grid', ...
    struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
    'solverOptions',struct('coreDegree',[32,32],'angularDegree',40,'polarDegree',12), ...
    'maxDt',.001);
c.stop.steady.enabled=true;c.stop.steady.sampleEveryTau=.001;
[opts,method]=ipm_cr.config(c,'run');[s0,ctx0]=ipm_cr.init(method);
assert(all(s0.eta.core==0,'all')&&all(s0.eta.polar==0,'all'),'ipm_cr:ReleaseEta0');
[~,~,refs]=ipm_cr.model.reference(ipm_cr.problem(method.initial));
assert(strcmp(opts.outputDir,c.run.outputDir)&&isequaln(s0.originalAuthority,refs));
assert(isequaln(ctx0.opts,method));
[a,s,ctx]=ipm_cr.run(c);
assert(strcmp(a.status,'maximum_steps')&&s.step==1&&s.scale.canonicalTime==.001);
assert(a.maximumAcceptedStageCfl<=method.cfl*(1+64*eps)&&~a.steadyCandidate);
[cp,pointer]=ipm_cr.io.checkpoint('load',c.run.outputDir);
assert(cp.state.step==s.step&&pointer.sequence==s.step&&~isempty(cp.stopState));
verify_package(a.dataPackage,s,ctx);
c.run.maxSteps=1000;[b,s,ctx]=ipm_cr.run(c);
assert(strcmp(b.status,'target_reached')&&b.resumed&&s.scale.canonicalTime==c.run.tauEnd);
assert(~b.steadyCandidate&&~b.steadyStateCertified&&~b.trajectoryAccuracyCertified);
verify_package(b.dataPackage,s,ctx);
assert(isfile(fullfile(c.run.outputDir,'history.csv'))&&isfile(fullfile(c.run.outputDir,'case.json')));
% Physical-time thresholds stop at the first accepted endpoint after crossing.
c.run.outputDir=fullfile(outputDir,'physical');c.stop.steady.enabled=false;
c.stop.physicalTime=realmin;[r,s]=ipm_cr.run(c);
assert(strcmp(r.status,'physical_time')&&s.step==1&&r.physicalTimeStopOvershoot>=0);
% A tiny invocation budget still commits the original accepted initial state.
c.run.outputDir=fullfile(outputDir,'budget');c.stop.physicalTime=[];c.run.wallSeconds=realmin;
[r,s]=ipm_cr.run(c);assert(strcmp(r.status,'wall_budget')&&s.step==0);
c.run.outputDir=fullfile(outputDir,'stopped');c.run.wallSeconds=120;
mkdir(c.run.outputDir);fid=fopen(fullfile(c.run.outputDir,'STOP'),'w');fclose(fid);
[r,s]=ipm_cr.run(c);assert(strcmp(r.status,'stopped')&&s.step==0);
% Deliberately loose thresholds test only the positive stop/resume control
% path. They are not physical steady-state evidence or production criteria.
c.run.outputDir=fullfile(outputDir,'candidate_control');
c.stop.steady.enabled=true;c.stop.steady.residualTol=1e9;c.stop.steady.changeTol=1e9;
c.stop.steady.windowTau=1e-4;c.stop.steady.sampleEveryTau=1e-4;c.stop.steady.windows=1;
[r,s]=ipm_cr.run(c);assert(strcmp(r.status,'steady_candidate')&&s.step==1);
[again,s2]=ipm_cr.run(c);
assert(strcmp(again.status,'steady_candidate')&&s2.step==s.step&&s2.scale.canonicalTime==s.scale.canonicalTime);
% An old accepted state cannot be relabeled with a later manifest identity.
[~,name]=fileparts(tempname);initial=make_initial(['release_',name],struct('widthX',.9));
custom=rmfield(method,'background');custom.initial=initial;[state,own]=ipm_cr.init(custom);
assert(all(state.eta.core==0,'all')&&all(state.eta.polar==0,'all'));
file=fullfile(root,initial);text=fileread(file);restore=onCleanup(@()restore_text(file,text));
fid=fopen(file,'w');fprintf(fid,'%s\n',text);fclose(fid);
expect(@()ipm_cr.pack(state,own,fullfile(outputDir,'changed_initial')),'ipm_cr:PackIdentity');
clear restore;
report=struct('kind','ipm_cr_release_integration_v1','passed',true, ...
    'defaultEtaExactlyZero',true,'physicalCflPreserved',true,'stopStateResumable',true, ...
    'nativePackageShaAndFields',true,'physicalTimeAcceptedBoundary',true, ...
    'budgetAndStopAcceptedOnly',true,'candidateResumeControlOnly',true, ...
    'exportRejectsChangedInitial',true,'steadyStateCertified',false);
ipm_cr.io.json(fullfile(outputDir,'release.json'),report);
end

function verify_package(info,s,ctx)
file=fullfile(info.folder,info.fieldFile);
assert(strcmp(ipm_cr.io.sha(file),info.fieldSHA256));v=load(file,'data');d=v.data;
assert(d.step==s.step&&isequaln(d.scale,s.scale)&&isequaln(d.method,ctx.opts));
for patch={'core','polar'}
    key=patch{1};assert(isequaln(d.fields.(key).W,s.W.(key))&& ...
        isequaln(d.fields.(key).F,s.rhs.F.(key))&&isequaln(d.fields.(key).Q,s.rhs.Q.(key))&& ...
        isequaln(d.mesh.(key).active,ctx.mesh.(key).active));
end
assert(strcmp(ipm_cr.io.sha(fullfile(info.folder,info.metadataFile)),info.metadataSHA256));
m=jsondecode(fileread(fullfile(info.folder,info.metadataFile)));
assert(m.gauge.I==s.rhs.measurement.I&&m.gauge.P==s.rhs.measurement.P&& ...
    m.gauge.C==s.rhs.measurement.C&&~m.steadyStateCertified);
end

function restore_text(file,text)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',text);
end

function expect(fn,id)
try,fn();error('ipm_cr:ExpectedFailure','Expected %s.',id);
catch problem,assert(strcmp(problem.identifier,id),'ipm_cr:WrongFailure', ...
        'Expected %s, got %s.',id,problem.identifier);end
end
