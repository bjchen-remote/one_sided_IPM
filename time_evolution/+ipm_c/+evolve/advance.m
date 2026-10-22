function [out,stopReason]=advance(state)
%ADVANCE One direct-CFL SSPRK3 step of centered density and origin amplitude.
% This path has no nonlinear solve, step retry or error-estimator branch.
% Automatic remeshing is a separate transaction after the explicit step.
% An exception leaves the caller's accepted state untouched.
out=state;stopReason='';time=state.config.time;
assert(strcmp(time.timeIntegrator,'ssprk3'),'ipm_c:DirectCflMethod', ...
    'The direct polar solver uses the shared explicit SSPRK3 stages.');
transportRate=state.flow.generatorCflRate;
assert(isnumeric(transportRate)&&isreal(transportRate)&&isscalar(transportRate)&& ...
    isfinite(transportRate)&&transportRate>=0, ...
    'ipm_c:DirectCflRate','The complete polar transport rate must be finite and nonnegative.');
assert(isnumeric(state.flow.c_omega)&&isreal(state.flow.c_omega)&& ...
    isscalar(state.flow.c_omega)&&isfinite(state.flow.c_omega), ...
    'ipm_c:DirectCflRate','The physical amplitude reaction rate must be finite.');
reactionRate=abs(state.flow.c_omega);rate=transportRate+reactionRate;
assert(isfinite(rate),'ipm_c:DirectCflRate','The combined transport and reaction rate must be finite.');
assert(isscalar(time.cfl)&&isfinite(time.cfl)&&time.cfl>0&& ...
    isscalar(time.maxDt)&&isfinite(time.maxDt)&&time.maxDt>0&& ...
    isscalar(time.minDt)&&isfinite(time.minDt)&&time.minDt>0, ...
    'ipm_c:DirectCflConfiguration','Use positive finite CFL, maximum and minimum steps.');
clockSpeed=exp(state.scale.logC_omega-state.scale.logC_l);
assert(isfinite(clockSpeed)&&clockSpeed>0,'ipm_c:DirectCflClock', ...
    'The physical clock rate must be finite and positive.');
canonicalRemaining=time.finalTime-state.scale.canonicalTime;
oldPhysicalLow=state.scale.physicalTimeLow;
physicalRemaining=(time.physicalFinalTime-state.scale.physicalTime)-oldPhysicalLow;
if physicalRemaining<=0,stopReason='physical_final_time';return;end
if canonicalRemaining<=0,stopReason='final_time';return;end
if rate>0,dtCfl=time.cfl/rate;else,dtCfl=Inf;end
limits=[time.maxDt,dtCfl,canonicalRemaining,physicalRemaining/clockSpeed];
[dt,index]=min(limits);names={'max_dt','cfl','canonical_final_time','physical_final_time'};
if dt<time.minDt||state.scale.canonicalTime+dt==state.scale.canonicalTime
    stopReason='time_step_below_minimum';return;
end
cache=[];if isfield(state,'rhsCache'),cache=state.rhsCache;end
% Deliberately one call, with no catch/retry. The workflow owns persistence.
[rho,flow,scale,rhsCache,info]=ipm.evolve.stepRk( ...
    state.rho,dt,state.ops,'ssprk3',state.scale,cache);
trial=state;trial.rho=rho;trial.flow=flow;trial.scale=scale;trial.rhsCache=rhsCache;
trial.step=state.step+1;trial.normalizedTime=scale.canonicalTime;
assert(ipm.evolve.isFinite(trial)&&isfield(flow,'fullQ')&& ...
    all(isfinite(flow.fullQ),'all')&&isfinite(flow.generatorCflRate)&&flow.generatorCflRate>=0, ...
    'ipm_c:DirectCflNonfinite','The complete explicit endpoint must remain finite.');
assert((scale.physicalTime>state.scale.physicalTime|| ...
    (scale.physicalTime==state.scale.physicalTime&&scale.physicalTimeLow>oldPhysicalLow))&& ...
    scale.canonicalTime>state.scale.canonicalTime&&all(rho(1,:)==0), ...
    'ipm_c:DirectCflClock','A direct step must advance both clocks.');
if isfinite(time.physicalFinalTime)
    assert((scale.physicalTime-time.physicalFinalTime)+scale.physicalTimeLow<= ...
        64*eps(max(1,time.physicalFinalTime)), ...
        'ipm_c:DirectCflPhysicalHorizon','The explicit endpoint exceeds the physical horizon.');
end
trial.ops.rescaling.gaugeRateGuess=[flow.c_l,flow.c_r];
trial.timeStep=struct('dt',dt,'method','ssprk3','rate',rate, ...
    'transportRate',transportRate,'reactionRate',reactionRate,'dtCfl',dtCfl, ...
    'maxDtLimit',limits(1),'canonicalLimit',limits(3),'physicalLimit',limits(4), ...
    'physicalClockSpeed',clockSpeed,'activeLimiter',names{index}, ...
    'realizedCfl',dt*rate,'terminalCfl',dt*(flow.generatorCflRate+abs(flow.c_omega)), ...
    'realizedCflScope','initial_accepted_full_transport_plus_absolute_comega', ...
    'rhsEvaluations',info.rhsEvaluations,'initialStageCacheHit',info.initialStageCacheHit, ...
    'terminalRhsEvaluated',info.terminalRhsEvaluated,'attemptedSteps',1,'retryCount',0, ...
    'explicitCflLimitApplied',true,'implicitStageSolves',0,'remeshApplied',false);
if ~isfield(trial,'runMetadata'),trial.runMetadata=struct();end
trial.runMetadata.lastExplicitStep=info;
if isfield(trial.config,'remesh')
    policy=trial.config.remesh;
    enabled=~isfield(policy,'adaptiveRemesh')||policy.adaptiveRemesh;
    if enabled
        if isfield(trial.ops,'spatialOperators')&&isfield(trial.ops.spatialOperators,'remesh')
            [trial,remeshInfo]=trial.ops.spatialOperators.remesh(trial);
        else
            [trial,remeshInfo]=ipm_c.remesh.apply(trial);
        end
        % apply owns transfer, fresh candidate audit, commit and rollback.
        % A rejected proposal retains this already accepted physical step.
        trial.timeStep.remeshApplied=remeshInfo.applied;
        trial.timeStep.remesh=remeshInfo;
        trial.timeStep.rhsEvaluationScope='explicit_step_before_remesh';
        trial.timeStep.terminalCfl=dt*(trial.flow.generatorCflRate+abs(trial.flow.c_omega));
    end
end
out=trial;
end
