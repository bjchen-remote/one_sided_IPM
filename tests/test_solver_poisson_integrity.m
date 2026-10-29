function report=test_solver_poisson_integrity(outputDir)
%TEST_SOLVER_POISSON_INTEGRITY Bounded Poisson diagnostics/regression tests.
% No time evolution, no continuum/tail certification, no private reporters.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
timer=tic;
cfg=struct('coreDegree',[16,16],'angularDegree',20,'polarDegree',10, ...
    'rOuter',4,'auditPoints',13,'requireAccuracy',false, ...
    'overlapValueTolerance',1,'overlapFluxTolerance',1);
s=ipm_cr.psolver(cfg);c=s.Core;p=s.Polar;
source=struct('core',source_value(c.X,c.Y),'polar',source_value(p.X,p.Y));
[solution,original]=s.solve(source);
queries=s.AuditQueries;
audit=struct('core',source_value(queries.core.x,queries.core.y), ...
    'polar',source_value(queries.polar.x,queries.polar.y),'solverId',s.Id);
withAudit=source;withAudit.nativeReconstructionAudit=audit;
[audited,info]=s.solve(withAudit);
assert(isequal(solution.core.psi,audited.core.psi)&&isequal(solution.polar.psi,audited.polar.psi), ...
    'ipm_cr:NewPoissonIntegrity','Diagnostics changed the numerical solution.');
assert(info.accepted==original.accepted&&info.originalAccuracyGatesPreserved&& ...
    info.nativeReconstructionDefect.available&&~info.strongResidualAvailable&& ...
    ~info.nativeReconstructionDefect.independentContinuumCertificate, ...
    'ipm_cr:NewPoissonIntegrity','Native reconstruction must remain a separately labelled diagnostic.');
vcore=ipm_cr.poisson.core.evaluate(c,solution.core,queries.core.x,queries.core.y);
vpolar=ipm_cr.poisson.polar.evaluate(p,solution.polar,queries.polar.x,queries.polar.y);
expected=max(max(abs(-vcore.laplacian-audit.core),[],'all'), ...
    max(abs(-vpolar.laplacian-audit.polar),[],'all'));
% Tensor regrouping and nominal polar coordinates can change roundoff.
% This tolerance is 5e-11 times the common Laplacian/source amplitude,
% rather than a bitwise requirement or a spatial-accuracy assertion.
roundoffFactor=5e-11;
referenceScale=max([1;abs(vcore.laplacian(:));abs(vpolar.laplacian(:)); ...
    abs(audit.core(:));abs(audit.polar(:))]);
roundoffTolerance=roundoffFactor*referenceScale;
auditPlan=ipm_cr.poisson.auditPlan(c,p,s.Config);
[tensorCore,tensorPolar]=ipm_cr.poisson.nativeLaplacian(c,p,solution,auditPlan);
tensorCoreDifference=max(abs(tensorCore-vcore.laplacian),[],'all');
tensorPolarDifference=max(abs(tensorPolar-vpolar.laplacian),[],'all');
assert(max([tensorCoreDifference,tensorPolarDifference, ...
    abs(info.nativeReconstructionDefect.absoluteMaximum-expected)])<=roundoffTolerance, ...
    'ipm_cr:NewPoissonIntegrity','Tensor and generic fixed-query audits differ beyond the declared rounding tolerance.');
genericPanelDefect=zeros(p.nPanels,1);
for panel=1:p.nPanels
    selected=queries.polarPanel==panel;
    genericPanelDefect(panel)=max(abs(-vpolar.laplacian(selected)-audit.polar(selected)));
end
assert(max(abs(info.nativeReconstructionDefect.polarPanelAbsoluteMaximum-genericPanelDefect))<=roundoffTolerance, ...
    'ipm_cr:NewPoissonIntegrity','Per-panel tensor defects differ from generic evaluation.');

% Cached and uncached interpolation must evaluate the same polynomials.
x=[-.9;-.31;.17;.91];y=[.11;.37;.49;.77];
cm=ipm_cr.poisson.core.prepare(c,x,y);
cached=ipm_cr.poisson.core.evaluate(c,solution.core,x,y,[],cm);
plain=ipm_cr.poisson.core.evaluate(c,solution.core,x,y);
assert(isequaln(cached,plain),'ipm_cr:NewPoissonIntegrity','Core cached evaluation differs.');
psi=ipm_cr.poisson.core.evaluate(c,solution.core,x,y,{'psi'},cm);
assert(isequal(fieldnames(psi),{'psi'})&&isequal(psi.psi,plain.psi), ...
    'ipm_cr:NewPoissonIntegrity','Core psi-only evaluation changed its value or returned unused derivatives.');
x=[.8,1.4,2.7;-.9,-1.1,-2.3];y=[.13,.4,.8;.21,.6,1.2];
pm=ipm_cr.poisson.polar.prepare(p,x,y);
cached=ipm_cr.poisson.polar.evaluate(p,solution.polar,x,y,[],pm);
plain=ipm_cr.poisson.polar.evaluate(p,solution.polar,x,y);
assert(isequaln(cached,plain),'ipm_cr:NewPoissonIntegrity','Polar cached evaluation differs.');
psi=ipm_cr.poisson.polar.evaluate(p,solution.polar,x,y,{'psi'},pm);
assert(isequal(fieldnames(psi),{'psi'})&&isequal(psi.psi,plain.psi), ...
    'ipm_cr:NewPoissonIntegrity','Polar psi-only evaluation changed its value or returned unused derivatives.');

% Compare sparse assembly with the former dense panel arithmetic, entrywise.
[Dt,Dtt,flux]=dense_radial(p);
assert(issparse(p.Dt)&&issparse(p.Dtt)&&issparse(p.flux)&& ...
    isequal(full(p.Dt),Dt)&&isequal(full(p.Dtt),Dtt)&&isequal(full(p.flux),flux), ...
    'ipm_cr:NewPoissonIntegrity','Sparse panel assembly changed matrix values.');
for k=1:numel(p.mu)
    A=-Dtt+p.mu(k)*eye(p.nr);A(p.interfaceRows,:)=flux(p.interfaceRows,:);
    A(1,:)=0;A(1,1)=1;A(end,:)=Dt(end,:);A(end,end)=A(end,end)+sqrt(p.mu(k));
    assert(isequal(full(p.operators{k}),A), ...
        'ipm_cr:NewPoissonIntegrity','A sparse modal operator differs from its dense reference.');
end

% Dyadic samples preserve an exact common-translation check in double precision.
a=struct('u1',[1;2],'u2',[.5;-.5]);
b=struct('u1',[1.125;2],'u2',[.5;-.375]);
m=ipm_cr.poisson.overlapVelocityMetrics(a,b,[0,0],.25);
a.u1=a.u1+2^20;b.u1=b.u1+2^20;
shifted=ipm_cr.poisson.overlapVelocityMetrics(a,b,[2^20,0],.25);
assert(m.absoluteMismatch==shifted.absoluteMismatch&& ...
    m.centeredLocalScaledMismatch==shifted.centeredLocalScaledMismatch&& ...
    shifted.legacyScaledMismatch<m.legacyScaledMismatch, ...
    'ipm_cr:NewPoissonIntegrity','The local centered diagnostic must not be hidden by common translation.');

% Every new gate supplements the legacy gate, and missing data cannot pass it.
strictCfg=cfg;strictCfg.requireNativeDefect=true;strictCfg.nativeDefectTolerance=1e-12;
strict=ipm_cr.psolver(strictCfg);
bad=withAudit;bad.nativeReconstructionAudit.solverId=strict.Id;
bad.nativeReconstructionAudit.core=bad.nativeReconstructionAudit.core+1;
[~,badInfo]=strict.solve(bad);[~,missingInfo]=strict.solve(source);
assert(badInfo.legacyAccuracyAccepted&&~badInfo.nativeDefectGateAccepted&&~badInfo.accepted&& ...
    ~missingInfo.nativeDefectGateAccepted&&~missingInfo.accepted, ...
    'ipm_cr:NewPoissonIntegrity','The native gate must reject a large defect and unavailable data.');
strictCfg=cfg;strictCfg.overlapAbsoluteVelocityTolerance=realmin;
strictCfg.overlapCenteredVelocityTolerance=realmin;
strict=ipm_cr.psolver(strictCfg);[~,velocityInfo]=strict.solve(source);
assert(velocityInfo.legacyAccuracyAccepted&&~velocityInfo.additionalVelocityGatesAccepted&& ...
    ~velocityInfo.accepted,'ipm_cr:NewPoissonIntegrity','Additional velocity gates were not enforced.');

% Exercise the production native maps and analytic offset, including the
% nonseparable return branch whose includePsi=false used to be ignored.
g=ipm_cr.mesh.build();adapterCfg=cfg;adapterCfg.rOuter=g.polar.r(end);
adapterCfg.coreHalfWidth=1.1;adapterCfg.coreHeight=1.1;
cache=ipm_cr.poisson.plan(g,adapterCfg);
native=struct('core',source_value(g.core.X,g.core.Y), ...
    'polar',source_value(g.polar.X,g.polar.Y));
analytic=@(x,y).125*source_value(x,y);
[flow,adapter]=ipm_cr.poisson.apply(native,g,cache,analytic,false);
assert(~cache.coreSeparable&&isempty(flow.core.psi)&&isempty(flow.polar.psi)&& ...
    ~adapter.returnIncludesPsi&&adapter.nativeReconstructionDefect.available&& ...
    adapter.nativeReconstructionDefectSourceMapsCached&& ...
    ~adapter.nativeReconstructionDefectIsIndependentContinuumResidual, ...
    'ipm_cr:NewPoissonIntegrity','Adapter diagnostics or includePsi=false were not respected.');
q=cache.solver.AuditQueries;
qc=reshape(cache.auditSourceCore.core*native.core(:)+cache.auditSourceCore.polar*native.polar(:),size(q.core.x));
qp=reshape(cache.auditSourcePolar.core*native.core(:)+cache.auditSourcePolar.polar*native.polar(:),size(q.polar.x));
qc=qc+analytic(q.core.x,q.core.y);qp=qp+analytic(q.polar.x,q.polar.y);
vc=ipm_cr.poisson.core.evaluate(cache.solver.Core,flow.solution.core,q.core.x,q.core.y,{'laplacian'});
vp=ipm_cr.poisson.polar.evaluate(cache.solver.Polar,flow.solution.polar,q.polar.x,q.polar.y,{'laplacian'});
adapterExpected=max(max(abs(-vc.laplacian-qc),[],'all'),max(abs(-vp.laplacian-qp),[],'all'));
adapterScale=max([1;abs(vc.laplacian(:));abs(vp.laplacian(:));abs(qc(:));abs(qp(:))]);
adapterTolerance=roundoffFactor*adapterScale;
assert(abs(adapter.nativeReconstructionDefect.absoluteMaximum-adapterExpected)<=adapterTolerance, ...
    'ipm_cr:NewPoissonIntegrity','Native-map or analytic-offset audit values were not included.');
coreExpected=max(abs(-vc.laplacian-qc),[],'all');
polarExpected=max(abs(-vp.laplacian-qp),[],'all');
assert(abs(adapter.nativeReconstructionDefect.coreAbsoluteMaximum-coreExpected)<=adapterTolerance&& ...
    abs(adapter.nativeReconstructionDefect.polarAbsoluteMaximum-polarExpected)<=adapterTolerance, ...
    'ipm_cr:NewPoissonIntegrity','Core/polar analytic-offset defect components differ from generic queries.');
for panel=1:cache.solver.Polar.nPanels
    selected=q.polarPanel==panel;
    expectedPanel=max(abs(-vp.laplacian(selected)-qp(selected)));
    assert(abs(adapter.nativeReconstructionDefect.polarPanelAbsoluteMaximum(panel)-expectedPanel)<=adapterTolerance, ...
        'ipm_cr:NewPoissonIntegrity','A native reconstructed analytic-offset panel defect differs.');
end
report=struct('passed',true,'kind','poisson_integrity_and_cached_audit_v1','seconds',toc(timer), ...
    'numericalSolutionUnchangedByDiagnostics',true,'defaultAcceptanceUnchanged',true, ...
    'cachedEvaluationIdentical',true,'sparseOperatorsIdenticalToDenseReference',true, ...
    'newGatesSupplementLegacy',true,'nonseparableIncludePsiHonored',true, ...
    'tensorAuditComparison',struct('relativeAmplitudeTolerance',roundoffFactor, ...
        'absoluteTolerance',roundoffTolerance,'coreMaximumDifference',tensorCoreDifference, ...
        'polarMaximumDifference',tensorPolarDifference,'adapterAbsoluteTolerance',adapterTolerance, ...
        'bitwiseEquivalenceRequired',false), ...
    'nativeDefect',adapter.nativeReconstructionDefect,'timeAdvanced',0, ...
    'continuumAccuracyCertified',false,'tailAccuracyCertified',false);
if nargin>0&&~isempty(outputDir)
    if ~isfolder(outputDir),mkdir(outputDir);end
    path=fullfile(outputDir,'poisson_integrity.json');
    assert(~isfile(path),'ipm_cr:NewPoissonIntegrity','Choose a fresh report directory.');
    fid=fopen(path,'w');assert(fid>=0,'ipm_cr:NewPoissonIntegrity','Cannot open report.');
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');
end
fprintf('Poisson integrity checks passed in %.3fs; native reconstruction defect %.3g.\n', ...
    report.seconds,report.nativeDefect.scaledMaximum);
end

function q=source_value(x,y)
q=y.*(8-4*(x.^2+y.^2)).*exp(-x.^2-y.^2);
end

function [Dt,Dtt,flux]=dense_radial(p)
Dt=zeros(p.nr);Dtt=Dt;flux=Dt;
for j=1:p.nPanels
    indices=p.panelIndices{j};D=p.panelD{j};
    Dt(indices,indices)=Dt(indices,indices)+D;
    Dtt(indices,indices)=Dtt(indices,indices)+D*D;
end
Dt(p.interfaceRows,:)=.5*Dt(p.interfaceRows,:);
Dtt(p.interfaceRows,:)=.5*Dtt(p.interfaceRows,:);
for j=1:p.nPanels-1
    row=j*p.degree+1;il=p.panelIndices{j};ir=p.panelIndices{j+1};
    flux(row,il)=flux(row,il)+p.panelD{j}(end,:);
    flux(row,ir)=flux(row,ir)-p.panelD{j+1}(1,:);
end
end
