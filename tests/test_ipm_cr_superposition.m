function report=test_ipm_cr_superposition(outputDir)
%TEST_IPM_CR_SUPERPOSITION Direct-system tests of cached linear responses.
% Uses the original Gaussian density's actual native differentiated source,
% nonzero OLD and NEW inner traces, and nonzero exterior right hand sides.
% Decay/Robin and Dirichlet systems are checked independently. This tests
% discrete linear equivalence and original componentwise backward errors;
% it makes no continuum spatial-accuracy or speed claim.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir)
    parent=fullfile(root,'results');if ~isfolder(parent),mkdir(parent);end
    [~,token]=fileparts(tempname(parent));
    outputDir=fullfile(root,'results',['superposition_',token]);
end
outputDir=char(outputDir);
assert(startsWith(outputDir,filesep),'ipm_cr:SuperpositionOutput','Use an absolute output directory.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'superposition.json')), ...
    'ipm_cr:SuperpositionPreserve','Choose a fresh test output directory.');
timer=tic;cfg=ipm_cr.problem();physics=cfg.physics;
density=@(X,Y)ipm_cr.model.initial( ...
    struct('X',X,'Y',Y,'ny',size(X,1),'nx',size(X,2)),physics);
origin=density(0,0);equivalenceTolerance=5e-12;
report=struct('kind','native_source_cached_trace_superposition_direct_system_test_v1', ...
    'status','running','passed',false,'equivalenceTolerance',equivalenceTolerance, ...
    'productionAccuracyGatesRelaxed',false,'timeAdvanced',0, ...
    'sourceDefinition','Original Gaussian density, centered at its analytic origin, differentiated on the native Cartesian and log-polar spectral grids.', ...
    'oldInnerTraceIsNonzero',true,'newInnerTraceIsNonzero',true, ...
    'componentwiseOriginalMatrixResidualsComputed',true, ...
    'continuumAccuracyCertified',false,'cases',struct([]));
try
    for boundary={'decay','dirichlet'}
        mode=boundary{1};s=ipm_cr.psolver(struct('coreDegree',[32,32], ...
            'angularDegree',40,'polarDegree',12,'rOuter',16, ...
            'outerBoundary',mode,'requireAccuracy',true));
        c=s.Core;p=s.Polar;
        rc=density(c.X,c.Y)-origin;rp=density(p.X,p.Y)-origin;
        qc=rc*c.Dx';
        qp=((p.Dt*rp).*cos(p.theta)-(rp*p.Dtheta').*sin(p.theta))./p.r;
        oldInner=.08*sin(p.theta)+.02*sin(2*p.theta);
        newInner=-.03*sin(p.theta)+.014*sin(3*p.theta)+.005*sin(5*p.theta);
        outer=.015*sin(p.theta)+.007*sin(3*p.theta);
        oldInner([1,end])=0;newInner([1,end])=0;outer([1,end])=0;
        base=ipm_cr.poisson.polar.solve(p,qp,oldInner,outer);
        cached=ipm_cr.poisson.polar.withInnerTrace(p,base,newInner);
        direct=ipm_cr.poisson.polar.solve(p,qp,newInner,outer);
        unchanged=ipm_cr.poisson.polar.withInnerTrace(p,base,oldInner);
        cachedZero=ipm_cr.poisson.polar.withInnerTrace(p,cached,zeros(1,p.nt));
        directZero=ipm_cr.poisson.polar.solve(p,qp,zeros(1,p.nt),outer);
        updated=polar_comparison(p,cached,direct);
        reset=polar_comparison(p,cachedZero,directZero);
        sameTraceIsIdentity=isequaln(base,unchanged);

        % Reuse the full homogeneous core fields already cached when the
        % Schur matrix was built. Compare with the original direct solve at
        % the IDENTICAL source and artificial boundary values.
        p0=ipm_cr.poisson.polar.solve(p,qp,[],outer);
        boundaryValue=ipm_cr.poisson.polar.evaluate(p,p0,c.traceX,c.traceY);
        c0=ipm_cr.poisson.core.solve(c,qc,boundaryValue.psi);
        lambda=newInner(p.interior)';
        coreCached=struct('psi',c0.psi+ ...
            reshape(s.Coupling.coreResponseMatrix*lambda,size(c.X)));
        coreDirect=ipm_cr.poisson.core.solve(c,qc, ...
            boundaryValue.psi+s.Coupling.boundaryMap*lambda);
        core=struct('psi',difference(coreCached.psi,coreDirect.psi), ...
            'cachedBackwardError',core_backward(c,qc,coreCached), ...
            'directBackwardError',core_backward(c,qc,coreDirect));
        [~,coupled]=s.solve(struct('core',qc,'polar',qp),outer);
        entry=struct('outerBoundary',mode,'nativeSourceInfinity',max(abs(qp),[],'all'), ...
            'sourceShape',size(qp),'oldInnerInfinity',max(abs(oldInner)), ...
            'newInnerInfinity',max(abs(newInner)),'outerRhsInfinity',max(abs(outer)), ...
            'changedNonzeroInner',updated,'changedBackToZeroInner',reset, ...
            'sameInnerTraceIsIdentity',sameTraceIsIdentity,'cachedCoreResponse',core, ...
            'strictCoupledSolverAudit',coupled, ...
            'harmonicResidualInfinity',p.harmonicResidualInfinity, ...
            'harmonicRefinementUpdates',p.harmonicRefinementUpdates);
        entry.fastTraceLinearGatePass=updated.cachedBackwardError.maximum<=s.Config.linearTolerance;
        entry.fastResetLinearGatePass=reset.cachedBackwardError.maximum<=s.Config.linearTolerance;
        entry.fastCoreLinearGatePass=core.cachedBackwardError.maximum<=s.Config.linearTolerance;
        % A helper-level cancellation failure must remain visible. The
        % production Solver owns the documented direct fallback and full
        % strict audit; this unit test does not relax either requirement.
        entry.helperWouldRequireDirectFallback=~(entry.fastTraceLinearGatePass&& ...
            entry.fastResetLinearGatePass&&entry.fastCoreLinearGatePass);
        entry.passed=sameTraceIsIdentity&&updated.rhsExactlyIdentical&&reset.rhsExactlyIdentical&& ...
            updated.modal.mixedRelativeInfinity<=equivalenceTolerance&& ...
            updated.psi.mixedRelativeInfinity<=equivalenceTolerance&& ...
            reset.modal.mixedRelativeInfinity<=equivalenceTolerance&& ...
            reset.psi.mixedRelativeInfinity<=equivalenceTolerance&& ...
            core.psi.mixedRelativeInfinity<=equivalenceTolerance&& ...
            updated.directBackwardError.maximum<=s.Config.linearTolerance&& ...
            reset.directBackwardError.maximum<=s.Config.linearTolerance&& ...
            core.directBackwardError.maximum<=s.Config.linearTolerance&&coupled.accepted;
        if isempty(report.cases),report.cases=entry;else,report.cases(end+1)=entry;end %#ok<AGROW>
        report.elapsedSeconds=toc(timer);write_report(outputDir,report);
        assert(entry.passed,'ipm_cr:SuperpositionEquivalence', ...
            '%s cached responses differ from the original direct systems or their strict audits failed.',mode);
        fprintf('CR SUPERPOSITION %s modal %.3e psi %.3e core %.3e fastBackward %.3e fallback=%d\n', ...
            mode,updated.modal.absoluteInfinity,updated.psi.absoluteInfinity,core.psi.absoluteInfinity, ...
            updated.cachedBackwardError.maximum,entry.helperWouldRequireDirectFallback);
    end
    report.status='completed';report.passed=all([report.cases.passed]);
    report.elapsedSeconds=toc(timer);write_report(outputDir,report);
catch problem
    report.status='failed';report.passed=false;report.elapsedSeconds=toc(timer);
    report.error=struct('identifier',problem.identifier,'message',problem.message,'stack',problem.stack);
    write_report(outputDir,report);rethrow(problem);
end
end

function result=polar_comparison(p,cached,direct)
result=struct('modal',difference(cached.modal,direct.modal), ...
    'psi',difference(cached.psi,direct.psi),'rhs',difference(cached.rhs,direct.rhs), ...
    'rhsExactlyIdentical',isequaln(cached.rhs,direct.rhs), ...
    'cachedBackwardError',polar_backward(p,cached), ...
    'directBackwardError',polar_backward(p,direct));
end

function result=polar_backward(p,solution)
componentwise=zeros(size(solution.modal));residualInfinity=zeros(1,numel(p.mu));
for k=1:numel(p.mu)
    A=p.operators{k};u=solution.modal(:,k);b=solution.rhs(:,k);
    residual=A*u-b;
    componentwise(:,k)=abs(residual)./max(abs(A)*abs(u)+abs(b),realmin);
    residualInfinity(k)=norm(residual,inf);
end
[maximum,index]=max(componentwise(:));[row,mode]=ind2sub(size(componentwise),index);
result=struct('maximum',maximum,'perModeMaximum',max(componentwise,[],1), ...
    'perModeResidualInfinity',residualInfinity,'componentwise',componentwise, ...
    'maximumRow',row,'maximumMode',mode, ...
    'definition','abs(A_original*u-rhs)./max(abs(A_original)*abs(u)+abs(rhs),realmin)');
end

function result=core_backward(c,source,solution)
u=solution.psi;residual=-c.Dyy*u-u*c.Dxx'-source;
denominator=abs(c.Dyy)*abs(u)+abs(u)*abs(c.Dxx')+abs(source);
componentwise=abs(residual(c.iy,c.ix))./max(denominator(c.iy,c.ix),realmin);
result=struct('maximum',max(componentwise,[],'all'),'componentwise',componentwise, ...
    'definition','Original Cartesian collocation interior componentwise Poisson backward error');
end

function result=difference(a,b)
error=max(abs(a-b),[],'all');scale=max(1,max(abs(b),[],'all'));
result=struct('absoluteInfinity',error,'mixedRelativeInfinity',error/scale, ...
    'directInfinityScale',max(abs(b),[],'all'));
end

function write_report(folder,report)
temporary=fullfile(folder,'superposition.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'superposition.mat'),'f');
assert(ok,'ipm_cr:SuperpositionReport','%s',message);
temporary=fullfile(folder,'superposition.json.writing');fid=fopen(temporary,'w');
assert(fid>=0,'ipm_cr:SuperpositionReport');cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'superposition.json'),'f');
assert(ok,'ipm_cr:SuperpositionReport','%s',message);
end
