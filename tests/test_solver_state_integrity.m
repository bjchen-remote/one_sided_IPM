function report=test_solver_state_integrity(outputDir)
%TEST_SOLVER_STATE_INTEGRITY Stage clocks, fresh endpoints and q time meaning.
% Uses two single short production steps and an independent tangent/source
% check. It does not qualify a long trajectory or a continuum error bound.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','solver_state_integrity');end
if ~isfolder(outputDir),mkdir(outputDir);end
timer=tic;
method=struct('geometry',struct('coreSpacing',.06, ...
    'grid',struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
    'solverOptions',struct('coreDegree',[48,48],'angularDegree',64, ...
    'polarDegree',16,'rOuter',1e10),'maxDt',1e-4);
[initial,ctx]=ipm_cr.init(method);g=ctx.mesh;
assert(ctx.opts.background.tailCorrection,'ipm_cr:StateIntegrityTail', ...
    'The clock regression requires the explicitly time-dependent tail background.');
clock=clock_trials(g);steps=struct([]);
for name={'rk4','ssprk3'}
    integrator=name{1};
    [s,info]=ipm_cr.step(initial,ctx,1e-4, ...
        struct('timeIntegrator',integrator,'maxDt',1e-4));
    [fresh,W]=ipm_cr.rhs(s,ctx,true);
    assert(isequaln(W,s.W)&&isequaln(fresh.q,s.rhs.q), ...
        'ipm_cr:StateIntegrityEndpoint', ...
        'The accepted scale, W and q must describe the same evaluated endpoint.');
    ferror=relative_field(fresh.F,s.rhs.F);qerror=relative_field(fresh.Q,s.rhs.Q);
    assert(ferror<2e-10&&qerror<2e-9, ...
        'ipm_cr:StateIntegrityFreshRhs', ...
        'Fresh endpoint evaluation differs beyond gauge-root roundoff.');
    assert(s.scale.canonicalTime==initial.scale.canonicalTime+info.dt&& ...
        info.realizedMaximumStageCfl<=info.configuredCfl*(1+64*eps)&& ...
        strcmp(info.timeIntegrator,integrator), ...
        'ipm_cr:StateIntegrityStep','A short accepted step violated its time or all-stage CFL contract.');
    row=struct('method',integrator,'dt',info.dt,'retries',info.retries, ...
        'clockHigh',s.scale.physicalTime,'clockLow',s.scale.physicalTimeLow, ...
        'freshWAndQSourceBitwiseIdentical',true,'freshFRelativeDifference',ferror, ...
        'freshQTimeDerivativeRelativeDifference',qerror, ...
        'maximumStageCfl',info.realizedMaximumStageCfl);
    if isempty(steps),steps=row;else,steps(end+1)=row;end %#ok<AGROW>
end
representation=source_tangent(initial,ctx);
assembly=complete_assembly(initial,ctx);
report=struct('kind','ipm_cr_state_integrity_v1','passed',true,'clock',clock, ...
    'shortSteps',steps,'representation',representation,'completeAssembly',assembly, ...
    'trajectoryAccuracyCertified',false,'continuumErrorCertified',false,'seconds',toc(timer));
ipm_cr.io.json(fullfile(outputDir,'solver_state_integrity.json'),report);
fprintf('CR STATE: compensated stages/endpoints, q time derivative and active time cancellation passed (%.3fs)\n',report.seconds);
end

function report=clock_trials(g)
% The increment is smaller than an ulp of HIGH. Recovering it from rounded
% z(3)-HIGH would return zero, so this directly distinguishes the old path.
previous=struct('logC_l',0,'logC_omega',0,'physicalTime',2^53, ...
    'physicalTimeLow',.125,'X_shift',0,'canonicalTime',.125,'densityOrigin',0);
z=[0;0;previous.physicalTime;0;.375;0];
unit=ipm_cr.time.stageScale(previous,z,.25);
assert(unit.physicalTime==2^53&&unit.physicalTimeLow==.375&& ...
    unit.canonicalTime==.375,'ipm_cr:StateIntegrityClock', ...
    'A physical increment below the high-clock ulp must be retained in the low component.');
assert(z(3)-previous.physicalTime==0,'ipm_cr:StateIntegrityClockFixture');
zero=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
base=struct('scale',previous,'step',0,'W',zero,'rhs',fake_rhs(previous,zero));
rows=struct([]);captured=[];
callbacks=struct('stage',@stage,'finish',@finish,'notePeak',@ignore);
for name={'rk4','ssprk3'}
    captured=[];
    [s,r]=ipm_cr.time.trial(base,.25,name{1},g,struct('cfl',.3),callbacks);
    if strcmp(name{1},'rk4')
        b=ipm_cr.time.rk4();b=b.b;increment=.25*(b(1)+b(2)+b(3)+b(4));
        expected=[.125,.125,.25,increment];
    else
        increment=.25*(1/6+1/6+2/3);expected=[.25,.125,increment];
    end
    assert(isequal(captured,expected)&&r.physicalClockIncrement==increment&& ...
        isequaln(s.scale,s.rhs.evaluatedScale)&&s.scale.physicalTimeLow==.375, ...
        'ipm_cr:StateIntegrityClockTableau', ...
        'Intermediate and endpoint RHS clocks must use the tableau increments committed by finish.');
    row=struct('method',name{1},'stageIncrements',captured,'weightedClockIncrement',increment, ...
        'endpointHigh',s.scale.physicalTime,'endpointLow',s.scale.physicalTimeLow, ...
        'evaluatedAndCommittedScaleBitwiseIdentical',true);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
report=struct('high',previous.physicalTime,'initialLow',previous.physicalTimeLow, ...
    'increment',.25,'roundedHighDifference',0,'records',rows);

    function [f,W]=stage(W,z,~,~,increment)
        scale=ipm_cr.time.stageScale(previous,z,increment);
        captured(end+1)=increment;f=fake_rhs(scale,zero);
    end
    function s=finish(old,W,rho,z,f,h,increment)
        s=old;s.W=rho;s.rhs=f;s.scale=ipm_cr.time.stageScale(old.scale,z,increment);
        assert(s.scale.canonicalTime==old.scale.canonicalTime+h, ...
            'ipm_cr:StateIntegrityCanonicalEndpoint');
        s.step=old.step+1; %#ok<INUSD>
    end
end

function f=fake_rhs(scale,zero)
f=struct('F',zero,'c_l',0,'c_omega',0,'c_r',0,'transportOrigin',0, ...
    'rate',0,'cflPosition',struct(),'gaugeResidual',zeros(3,1),'evaluatedScale',scale);
end
function ignore(varargin) %#ok<INUSD>
end

function report=source_tangent(s,ctx)
% Differentiate the defining q representation, rather than the assembled F.
% The eta tangent is linear, while the analytic background is evaluated at
% four perturbed scales. Freezing the instantaneous rates is sufficient for
% the first derivative and avoids another Poisson/gauge solve in this check.
g=ctx.mesh;f=s.rhs;model=ctx.opts.background;scale=s.scale;
speed=exp(scale.logC_omega-scale.logC_l);
k=[f.c_l;f.c_omega;speed;f.c_l*scale.X_shift+f.c_r;1; ...
    f.transportOrigin+f.c_omega*scale.densityOrigin];
hs=[1e-3,5e-4,2.5e-4];errors=zeros(size(hs));
for level=1:numel(hs)
    h=hs(level);minus2=source(-2*h);minus=source(-h);plus=source(h);plus2=source(2*h);
    derivative=struct();
    for patch={'core','polar'}
        name=patch{1};derivative.(name)=(minus2.(name)-8*minus.(name)+ ...
            8*plus.(name)-plus2.(name))/(12*h);
    end
    errors(level)=relative_field(derivative,f.Q);
end
assert(errors(end)<2e-7,'ipm_cr:StateIntegrityQTangent', ...
    'Q must match the time derivative of q=D_h eta+H_X.');
assert(isfield(f,'representationAudit')&&f.representationAudit.available&& ...
    ~f.representationAudit.continuumErrorCertified&& ...
    strcmp(f.representationAudit.QDefinition,'time derivative of q=D_h eta+H_X'), ...
    'ipm_cr:StateIntegrityRepresentation','The Q representation must be stated and independently audited.');
d=ipm_cr.mesh.diff(f.F,g);defect=zeros(1,2);i=0;
for patch={'core','polar'}
    i=i+1;name=patch{1};active=g.(name).active;
    defect(i)=max(abs(f.Q.(name)(active)-d.(name).x(active)));
    assert(defect(i)==f.representationAudit.patches.(name).absoluteInfinity, ...
        'ipm_cr:StateIntegrityDefect','The reported discrete representation defect must equal its named quantity.');
end
assert(max(defect)==f.representationAudit.maximumAbsoluteDefect, ...
    'ipm_cr:StateIntegrityDefect');
report=struct('finiteDifferenceSteps',hs,'sourceTimeDerivativeRelativeErrors',errors, ...
    'representationDefectByPatch',defect,'auditMatchesMeasuredDefectExactly',true, ...
    'QIsDefinedAsSourceTimeDerivative',true,'continuumErrorCertified',false);

    function q=source(h)
        moved=s;moved.eta=struct('core',s.eta.core+h*f.residualRhs.core, ...
            'polar',s.eta.polar+h*f.residualRhs.polar);
        moved.scale=scale;moved.scale.logC_l=scale.logC_l+h*k(1);
        moved.scale.logC_omega=scale.logC_omega+h*k(2);
        moved.scale.physicalTime=scale.physicalTime+h*k(3);
        moved.scale.X_shift=scale.X_shift+h*k(4);
        moved.scale.canonicalTime=scale.canonicalTime+h;
        moved.scale.densityOrigin=scale.densityOrigin+h*k(6);
        eta=ipm_cr.mesh.fill(moved.eta,g);deta=ipm_cr.mesh.diff(eta,g);
        numerical=ipm_cr.mesh.fill(struct('core',deta.core.x,'polar',deta.polar.x),g);
        q=struct();
        for patch={'core','polar'}
            name=patch{1};p=g.(name);
            q.(name)=numerical.(name)+ipm_cr.model.source(model,moved.scale,p.X,p.Y);
        end
    end
end

function report=complete_assembly(s,ctx)
% On active nodes the full transport equation has no explicit timeH term:
% F=transport(eta)-u.grad(H)+cw*eta-T0 + gauge/background terms.
g=ctx.mesh;f=s.rhs;eta=ipm_cr.mesh.fill(s.eta,g);bg=struct();forcing=struct();
vx=struct();vy=struct();
for patch={'core','polar'}
    name=patch{1};p=g.(name);v=f.physical.(name);
    bg.(name)=ipm_cr.model.background(ctx.opts.background,s.scale,p.X,p.Y);
    forcing.(name)=-v.u1.*bg.(name).x-v.u2.*bg.(name).y;
    vx.(name)=(v.u1-f.Uref)+f.c_l*p.X+f.delta;
    vy.(name)=v.u2+f.c_l*p.Y;
end
prepared=ipm_cr.op.prepare(eta,g,struct('transport',ctx.opts.transport));
transport=prepared(vx,vy,struct('cfl',false));
T=ipm_cr.mesh.fill(struct('core',transport.B.core+forcing.core, ...
    'polar',transport.B.polar+forcing.polar),g);
T0=T.core(1,g.core.originIndex);maximum=0;timeMagnitude=0;
for patch={'core','polar'}
    name=patch{1};p=g.(name);b=bg.(name);active=p.active;
    base=T.(name)+f.c_omega*eta.(name)-T0;
    gaugeX=f.c_l*p.X+f.c_r;gaugeY=f.c_l*p.Y;
    expected=base+(-gaugeX.*b.x-gaugeY.*b.y+f.c_omega*b.H+f.c_r*b.originX);
    if strcmp(name,'core'),expected(1,g.core.originIndex)=0;end
    error=max(abs(expected(active)-f.F.(name)(active)));
    maximum=max(maximum,error);timeMagnitude=max(timeMagnitude,max(abs(b.timeH(active))));
    assert(isequaln(expected(active),f.F.(name)(active)), ...
        'ipm_cr:StateIntegrityCompleteAssembly', ...
        'Active complete F must be assembled without subtracting and then readding explicit timeH.');
end
assert(timeMagnitude>0,'ipm_cr:StateIntegrityTimeFixture', ...
    'The regression must actually exercise a nonzero explicit background time term.');
report=struct('activeFormulaBitwiseIdentical',true,'maximumDifference',maximum, ...
    'explicitBackgroundTimeTermMaximum',timeMagnitude, ...
    'eliminationRoundoffDiagnostic',f.representationAudit.eliminatedTimeCancellationDifference);
end

function error=relative_field(a,b)
error=0;
for patch={'core','polar'}
    name=patch{1};error=max(error,max(abs(a.(name)-b.(name)),[],'all')/ ...
        max(1,max(abs(b.(name)),[],'all')));
end
end
