function report=test_steady(outputDir)
%TEST_STEADY Synthetic state-machine tests; no PDE solve or native module.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','steady_test');end
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'steady.json')),'ipm_cr:SteadyTestOutput','Use a fresh output directory.');
timer=tic;
report=struct('kind','ipm_cr_synthetic_steady_logic_v1','passed',false,'status','running', ...
    'synthetic',true,'numericalSolverRun',false,'nativeRequired',false,'checks',struct(), ...
    'steadyStateCertified',false,'trajectoryAccuracyCertified',false, ...
    'scope','Synthetic arrays exercise candidate bookkeeping, norms and rejection. They are not a physical steady solution or numerical accuracy evidence.');
try
    [base,ctx]=initial();p=struct('windowTau',1);
    s=base;s.rhs.F.core(1,1)=2e-7;s.rhs.Q.polar(3)=8e-7;
    [~,r]=ipm_cr.diag.steady(s,ctx,p);
    check(r.regions.global.sampleCount==9&&r.regions.near.sampleCount==8&& ...
        r.regions.body.sampleCount==8&&r.regions.global.duplicateActiveOverlapRetained,'overlap copies');
    check(r.regions.global.F.denominator==2&&r.regions.global.Q.denominator==4&& ...
        r.regions.global.F.mixedInfinity==1e-7&&r.regions.global.Q.mixedInfinity==2e-7,'F/Q normalization');
    check(r.regions.global.Q.maximumLocation.X==201&& ...
        strcmp(r.regions.global.Q.maximumLocation.patch,'polar'),'native peak location');
    report.checks.nativeNormsAndOverlap=true;
    [t,r]=ipm_cr.diag.steady(base,ctx,p);
    check(~r.candidate&&~r.window.performed,'initial observation');
    for j=1:3
        s=at(base,j,j);[t,r]=ipm_cr.diag.steady(s,ctx,p,t);
        check(r.window.performed&&r.window.passed&&r.consecutive==j,'nonoverlapping window');
        check(r.candidate==(j==3),'three-window candidate');
    end
    saved=t;
    for j=1:5
        [t,repeat]=ipm_cr.diag.steady(s,ctx,p,t);
        check(repeat.duplicate&&repeat.consecutive==3&&isequaln(t,saved),'duplicate idempotence');
    end
    report.checks.nonoverlappingWindowsAndDuplicateIdempotence=true;
    s=at(base,4,3.25);[t,r]=ipm_cr.diag.steady(s,ctx,p,t);
    check(~r.candidate&&~r.window.performed&&r.consecutive==3,'candidate is not sticky');
    s=at(base,5,3.5);s.rhs.F.polar(3)=1e-3;
    [t,r]=ipm_cr.diag.steady(s,ctx,p,t);
    check(~r.currentResidualPassed&&~r.candidate&&r.consecutive==0&& ...
        r.regions.near.F.mixedInfinity==0&&r.regions.body.F.mixedInfinity==0&& ...
        r.regions.global.F.mixedInfinity>1e-6,'unmasked far F');
    s=at(base,6,4);[t,r]=ipm_cr.diag.steady(s,ctx,p,t);
    check(r.currentResidualPassed&&r.window.performed&&~r.window.passed&& ...
        ~r.window.residualPassed&&r.consecutive==0,'intermediate residual failure retained');
    report.checks.nonstickyCandidateAndWholeWindowResiduals=true;
    s=base;s.rhs.Q.polar(3)=1e-3;[~,r]=ipm_cr.diag.steady(s,ctx,p);
    check(~r.currentResidualPassed&&r.regions.global.F.mixedInfinity==0&& ...
        r.regions.near.Q.mixedInfinity==0&&r.regions.body.Q.mixedInfinity==0&& ...
        r.regions.global.Q.mixedInfinity>1e-6,'unmasked far Q');
    report.checks.fullFAndQBothRequired=true;
    p2=struct('windowTau',1,'windows',1,'minTau',3);
    [u,~]=ipm_cr.diag.steady(base,ctx,p2);
    [u,r]=ipm_cr.diag.steady(at(base,1,1),ctx,p2,u);
    check(r.window.passed&&~r.candidate,'minimum tau');
    [~,r]=ipm_cr.diag.steady(at(base,2,3),ctx,p2,u);
    check(r.candidate&&r.consecutive==2,'minimum tau reached');
    [u,~]=ipm_cr.diag.steady(base,ctx,p);
    [~,r]=ipm_cr.diag.steady(at(base,1,10),ctx,p,u);
    check(r.completedWindows==1&&r.consecutive==1&&~r.candidate,'one long interval is one window');
    report.checks.minimumTauAndNoInventedWindows=true;
    p3=struct('windowTau',1,'windows',1,'residualTol',1e-3,'changeTol',1e-8);
    for field={'W','q','F','Q','u1','u2'}
        f=field{1};[u,~]=ipm_cr.diag.steady(base,ctx,p3);s=at(base,1,1);
        if strcmp(f,'W'),s.W.core(1,1)=s.W.core(1,1)+1e-4;
        elseif ismember(f,{'q','F','Q'}),s.rhs.(f).core(1,1)=s.rhs.(f).core(1,1)+1e-4;
        else,s.rhs.physical.core.(f)(1,1)=s.rhs.physical.core.(f)(1,1)+1e-4;end
        [~,r]=ipm_cr.diag.steady(s,ctx,p3,u);
        check(r.currentResidualPassed&&~r.window.passed&& ...
            r.window.changes.(f).mixedInfinity>p3.changeTol,['window field ',f]);
    end
    [u,~]=ipm_cr.diag.steady(base,ctx,p3);s=at(base,1,1);s.rhs.Uref=s.rhs.Uref+1e-4;
    [~,r]=ipm_cr.diag.steady(s,ctx,p3,u);
    check(r.window.changes.u1.absoluteInfinity==0&&r.window.changes.u1Centered.absoluteInfinity>0&& ...
        ~r.window.passed,'physical and Uref-subtracted velocity both retained');
    for field={'c_l','c_omega','c_r','delta','Uref'}
        f=field{1};[u,~]=ipm_cr.diag.steady(base,ctx,p3);s=at(base,1,1);s.rhs.(f)=s.rhs.(f)+1e-4;
        [~,r]=ipm_cr.diag.steady(s,ctx,p3,u);
        check(~r.window.passed&&r.window.changes.rates.mixedInfinity>p3.changeTol,['rate ',f]);
    end
    report.checks.allRequestedFieldsAndRatesCompared=true;
    a=base;a.W.core(1,1)=2.1;[u,~]=ipm_cr.diag.steady(a,ctx,p3);
    b=at(base,1,1);b.W.core(1,2)=2.1;
    [~,r]=ipm_cr.diag.steady(b,ctx,p3,u);
    check(~r.window.passed&&r.window.changes.W.absoluteInfinity>.09,'no shift or amplitude fit');
    report.checks.noAlignmentOrFit=true;
    [u,~]=ipm_cr.diag.steady(base,ctx,p);
    changed=p;changed.windowTau=2;
    expect(@()ipm_cr.diag.steady(base,ctx,changed,u),'ipm_cr:SteadyPolicyChanged');
    other=ctx;other.mesh.core.X(1,1)=-2;
    expect(@()ipm_cr.diag.steady(base,other,p,u),'ipm_cr:SteadySetupChanged');
    bad=base;bad.W.core(1,1)=3;
    expect(@()ipm_cr.diag.steady(bad,ctx,p,u),'ipm_cr:SteadySequence');
    [u,~]=ipm_cr.diag.steady(at(base,3,3),ctx,p);
    expect(@()ipm_cr.diag.steady(at(base,4,2),ctx,p,u),'ipm_cr:SteadySequence');
    report.checks.policySetupAndSequenceRejected=true;
    bad=base;bad.rhs.Q.polar(3)=NaN;
    expect(@()ipm_cr.diag.steady(bad,ctx,p),'ipm_cr:SteadyFinite');
    bad=base;bad.rhs.physical.core.u1(1,1)=Inf;
    expect(@()ipm_cr.diag.steady(bad,ctx,p),'ipm_cr:SteadyFinite');
    bad=base;bad.rhs.c_l=NaN;
    expect(@()ipm_cr.diag.steady(bad,ctx,p),'ipm_cr:SteadyFinite');
    bad=base;bad.scale.canonicalTime=Inf;
    expect(@()ipm_cr.diag.steady(bad,ctx,p),'ipm_cr:SteadyFinite');
    expect(@()ipm_cr.diag.steady(base,ctx,struct('windows',0)),'ipm_cr:SteadyPolicy');
    check(~r.steadyStateCertified&&~r.trajectoryAccuracyCertified&&~r.gridConvergenceCertified&& ...
        ~r.candidateIsAccuracyCertificate,'candidate is not a certificate');
    report.checks.nonfiniteRejectedAndNoCertificate=true;
    report.passed=true;report.status='passed_synthetic_candidate_logic';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'steady.json'),report);rethrow(problem);
end
report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'steady.json'),report);
fprintf('CR STEADY: synthetic native norms/windows/rejections passed; no PDE qualification (%.3fs)\n',report.seconds);
end

function [s,ctx]=initial()
g=struct('core',struct('X',[-1,0,1;-1,0,1],'Y',[0,0,0;1,1,1],'active',true(2,3)), ...
    'polar',struct('X',[1,3,201],'Y',[1,0,0],'active',true(1,3)));
field=@(v)struct('core',v*ones(size(g.core.X)),'polar',v*ones(size(g.polar.X)));
physical=struct('core',struct('u1',3*ones(2,3),'u2',-.5*ones(2,3)), ...
    'polar',struct('u1',3*ones(1,3),'u2',-.5*ones(1,3)));
s=struct('W',field(2),'scale',struct('canonicalTime',0),'step',0, ...
    'rhs',struct('q',field(4),'F',field(0),'Q',field(0),'physical',physical, ...
    'c_l',0,'c_omega',0,'c_r',0,'delta',0,'Uref',1));
ctx=struct('mesh',g);
end

function s=at(s,step,tau)
s.step=step;s.scale.canonicalTime=tau;
end

function check(condition,label)
assert(condition,'ipm_cr:SteadyTest','Synthetic candidate test failed: %s.',label);
end

function expect(callback,identifier)
try
    callback();error('ipm_cr:SteadyTestExpected','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:SteadyTestRejection', ...
        'Expected %s, received %s.',identifier,problem.identifier);
end
end
