function report=test_direct_cfl()
%TEST_DIRECT_CFL Direct step arithmetic, no retry, SSPRK3 and one origin.
% Manufactured autonomous density plus a manufactured scale/clock ODE;
% no physical Poisson, gauge or PDE trajectory is advanced.
[state,control]=fixture(20,false);control.reset();
[advanced,reason]=ipm_c.evolve.advance(state);
expected=.4/(20+.7);
assert(isempty(reason)&&abs(advanced.timeStep.dt-expected)<=4*eps(expected)&& ...
    advanced.timeStep.rate==20+.7&&advanced.timeStep.transportRate==20&& ...
    advanced.timeStep.reactionRate==.7&& ...
    abs(advanced.timeStep.realizedCfl-.4)<=8*eps(.4), ...
    'ipm_c:DirectCflArithmetic','CFL must include full transport and abs(c_omega), including negative reaction.');
assert(state.flow.c_l~=0&&abs(2*state.flow.c_l+state.flow.c_omega)~=.7, ...
    'ipm_c:DirectCflReactionFixture','The fixture must distinguish c_omega from the conservative reaction.');
assert(control.count()==3&&advanced.timeStep.rhsEvaluations==3&& ...
    advanced.timeStep.initialStageCacheHit&&advanced.timeStep.retryCount==0&& ...
    advanced.timeStep.attemptedSteps==1&&advanced.timeStep.implicitStageSolves==0&& ...
    ~advanced.timeStep.remeshApplied&&advanced.ops.remeshCount==state.ops.remeshCount, ...
    'ipm_c:DirectCflSingleStep','A cached SSPRK3 step must evaluate only its two later stages and fresh endpoint.');
assert(all(advanced.rho(1,:)==advanced.rho(1,1)),'ipm_c:DirectCflOrigin');
noCache=rmfield(state,'rhsCache');control.reset();
[uncached,~]=ipm_c.evolve.advance(noCache);
assert(control.count()==4&&~uncached.timeStep.initialStageCacheHit&& ...
    uncached.timeStep.rhsEvaluations==4&&isequaln(advanced.rho,uncached.rho)&& ...
    isequaln(advanced.scale,uncached.scale),'ipm_c:DirectCflCache');

canonical=state;canonical.config.time.finalTime=.013;
[short,~]=ipm_c.evolve.advance(canonical);
assert(short.timeStep.dt==.013&&strcmp(short.timeStep.activeLimiter,'canonical_final_time'));
physical=state;physical.config.time.physicalFinalTime=.012;
[short,~]=ipm_c.evolve.advance(physical);
assert(short.timeStep.dt==.012&&strcmp(short.timeStep.activeLimiter,'physical_final_time'));
[zeroRate,~]=fixture(0,false,0);
zeroRate.config.time.maxDt=.01;
[short,~]=ipm_c.evolve.advance(zeroRate);
assert(short.timeStep.dt==.01&&isinf(short.timeStep.dtCfl)&& ...
    strcmp(short.timeStep.activeLimiter,'max_dt'));

saved=state;control.setFailure(true);control.reset();caught=false;
try
    [saved,~]=ipm_c.evolve.advance(state);
catch problem
    assert(strcmp(problem.identifier,'ipm_c:DirectCflManufacturedFailure'));caught=true;
end
assert(caught&&control.count()==1&&isequaln(saved,state), ...
    'ipm_c:DirectCflNoRetry','A failed stage must propagate once and leave the accepted state unchanged.');
control.setFailure(false);

enabled=state;enabled.config.remesh=struct('adaptiveRemesh',true);
control.reset();control.setRemeshRejected(false);
[observed,~]=ipm_c.evolve.advance(enabled);
assert(control.remeshCount()==1&&control.count()==3&& ...
    ~observed.timeStep.remeshApplied&&strcmp(observed.timeStep.remesh.status,'not_triggered')&& ...
    isequaln(observed.rho,advanced.rho)&&isequaln(observed.scale,advanced.scale)&& ...
    isequaln(observed.rhsCache,advanced.rhsCache), ...
    'ipm_c:DirectCflRemeshObservation','An enabled observer must run once without recomputing an unchanged RHS.');
control.reset();control.setRemeshRejected(true);
[rejectedMesh,~]=ipm_c.evolve.advance(enabled);
assert(control.remeshCount()==1&&control.count()==3&& ...
    strcmp(rejectedMesh.timeStep.remesh.status,'candidate_rejected')&& ...
    ~rejectedMesh.timeStep.remeshApplied&&rejectedMesh.step==advanced.step&& ...
    isequaln(rejectedMesh.rho,advanced.rho)&&isequaln(rejectedMesh.flow,advanced.flow)&& ...
    isequaln(rejectedMesh.scale,advanced.scale)&&isequaln(rejectedMesh.rhsCache,advanced.rhsCache), ...
    'ipm_c:DirectCflRemeshRollback','A rejected mesh must retain the accepted explicit step, clocks and RHS cache.');
disabled=enabled;disabled.config.remesh.adaptiveRemesh=false;control.reset();
[~,~]=ipm_c.evolve.advance(disabled);
assert(control.remeshCount()==0&&control.count()==3,'ipm_c:DirectCflDisabledRemesh');
control.setRemeshRejected(false);

counts=[4,8,16];T=.4;densityErrors=zeros(1,3);clockErrors=densityErrors;rows=cell(1,3);
for level=1:numel(counts)
    [s,c]=fixture(0,true);initial=s.rho;s.config.time.finalTime=T;
    s.config.time.maxDt=T/counts(level);c.reset();
    for step=1:counts(level)
        [s,reason]=ipm_c.evolve.advance(s);assert(isempty(reason));
        assert(all(s.rho(1,:)==s.rho(1,1))&&s.timeStep.initialStageCacheHit);
    end
    densityErrors(level)=max(abs(s.rho-exp(-.7*T)*initial),[],'all');
    exactClock=exp(.7^2/4)*sqrt(pi)/2*(erf(T+.7/2)-erf(.7/2));
    clockErrors(level)=abs(s.scale.physicalTime-exactClock);
    assert(c.count()==3*counts(level)&&s.step==counts(level));
    rows{level}=struct('substeps',counts(level),'densityError',densityErrors(level), ...
        'physicalClockError',clockErrors(level),'rhsEvaluations',c.count());
end
densityOrders=log2(densityErrors(1:end-1)./densityErrors(2:end));
clockOrders=log2(clockErrors(1:end-1)./clockErrors(2:end));
assert(all(densityOrders>2.7)&&all(densityOrders<3.3)&& ...
    all(clockOrders>2.7)&&all(clockOrders<3.3), ...
    'ipm_c:DirectCflThirdOrder','Shared SSPRK3 must retain third-order density and physical-clock convergence.');
report=struct('kind','single_direct_CFL_SSPRK3_polar_state_manufacture_v1','passed',true, ...
    'negativeComOmegaIncluded',true,'conservativeReactionNotUsed',true, ...
    'oneFailedNativeCallWithoutRetry',true,'acceptedStatePreservedOnFailure',true, ...
    'uniqueOrigin',true,'cachedAndUncachedEquivalent',true,'records',{rows}, ...
    'densityOrders',densityOrders,'physicalClockOrders',clockOrders, ...
    'physicalPdeAdvanced',false,'implicitCallbackUsed',false, ...
    'enabledRemeshObservedOnce',true,'rejectedRemeshPreservesAcceptedStep',true);
end

function [state,control]=fixture(transportRate,orderMode,baseCw)
if nargin<3,baseCw=-.7;end
nr=4;nt=7;count=0;fail=false;remeshCount=0;rejectRemesh=false;
rho=[repmat(.3,1,nt);reshape(linspace(.1,.9,(nr-1)*nt),nr-1,nt)];
ops=struct('dynamicScaleGeometry','isotropic','remeshCount',0);
ops.rescaling=struct('enabled',true,'lengthGauge','wall_gradient_integral', ...
    'cOmegaGauge','wall_omega_global_spline_maximum','gaugeRateGuess',[0,0]);
ops.spatialGeometryKey=struct('r',linspace(0,1,nr)', ...
    'theta',linspace(0,pi,nt),'operatorVersion',1,'remeshVersion',0);
ops.spatialOperators=struct('evaluate',@evaluate,'derivativeX',@(R)R, ...
    'pack',@(R)pack(R,nr,nt),'unpack',@(x)unpack(x,nr,nt), ...
    'preconditioner',@forbidden,'remesh',@observe_remesh);
scale=ipm.evolve.initialScale(ops);[F,flow]=evaluate(rho,ops,scale);
time=struct('timeIntegrator','ssprk3','cfl',.4,'maxDt',.1,'minDt',1e-12, ...
    'finalTime',1,'physicalFinalTime',Inf);
state=struct('rho',rho,'flow',flow,'scale',scale,'ops',ops,'config',struct('time',time), ...
    'step',0,'normalizedTime',0,'runMetadata',struct());
state.rhsCache=ipm.evolve.makeRhsCache(rho,F,flow,ops,scale);
control=struct('count',@get_count,'reset',@reset,'setFailure',@set_failure, ...
    'remeshCount',@get_remesh_count,'setRemeshRejected',@set_remesh_rejected);
    function [F,flow]=evaluate(value,~,s)
        count=count+1;
        if fail,error('ipm_c:DirectCflManufacturedFailure','Manufactured single-stage failure.');end
        assert(all(value(1,:)==value(1,1)),'ipm_c:DirectCflOrigin');
        F=-.7*value;cw=baseCw;cl=.3;
        if orderMode,cw=-.7-2*s.canonicalTime;cl=0;end
        zero=zeros(size(value));
        flow=struct('u1',zero,'u2',zero,'psi',zero,'source',value, ...
            'transportU1',zero,'transportU2',zero,'poissonResidual',0,'timeSpeed',1, ...
            'c_l',cl,'c_x',cl,'c_y',cl,'c_omega',cw,'c_r',0, ...
            'canonicalCL',cl,'canonicalCX',cl,'canonicalCY',cl,'canonicalCOmega',cw, ...
            'canonicalCR',0,'conservativeSource',2*cl+cw,'canonicalConservativeSource',2*cl+cw, ...
            'fullQ',F,'generatorCflRate',transportRate);
    end
    function value=get_count(),value=count;end
    function value=get_remesh_count(),value=remeshCount;end
    function reset(),count=0;remeshCount=0;end
    function set_failure(value),fail=value;end
    function set_remesh_rejected(value),rejectRemesh=value;end
    function [accepted,info]=observe_remesh(accepted)
        remeshCount=remeshCount+1;
        info=struct('applied',false,'requested',rejectRemesh,'pending',rejectRemesh, ...
            'status','not_triggered','manufacturedFixture',true);
        if rejectRemesh,info.status='candidate_rejected';end
        accepted.ops.rescaling.manufacturedRemeshObserver=remeshCount;
    end
end
function x=pack(R,nr,nt)
assert(isequal(size(R),[nr,nt])&&all(R(1,:)==R(1,1)));
positive=R(2:end,:).';x=[R(1,1);positive(:)];
end
function R=unpack(x,nr,nt)
assert(isequal(size(x),[1+(nr-1)*nt,1]));
R=[repmat(x(1),1,nt);reshape(x(2:end),nt,nr-1).'];
end
function varargout=forbidden(varargin) %#ok<STOUT,INUSD>
error('ipm_c:DirectCflUnexpectedSolver','Implicit callbacks must never run.');
end
