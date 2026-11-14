function [out,info]=step(state,ctx,finalTau,controls)
%STEP One explicit step; numerical verification belongs to the test suite.
% CFL selects the initial interval. RK4 partitions ordinary continuous peak
% crossings of original spline knots by up to four secant iterates of S_x=0.
% Numerical diagnostics never accept/reject a state or select a fallback.
if nargin<4,controls=struct();end
policy=struct('cfl',ctx.opts.cfl,'maxDt',ctx.opts.maxDt);
assert(isstruct(controls)&&isscalar(controls)&&all(ismember(fieldnames(controls),fieldnames(policy))), ...
    'ipm_cr:StepControls','Use CFL and maxDt step controls.');
for name=fieldnames(controls)',policy.(name{1})=controls.(name{1});end
h=min([policy.maxDt,.95*policy.cfl/state.rhs.rate, ...
    finalTau-state.scale.canonicalTime]);
initialH=h;probeH=h;hasResidual=isfield(state,'eta');calls=0;peaks=struct([]);
callbacks=struct('stage',@stage,'finish',@finish,'notePeak',@note_peak);
started=tic;
[out,records]=trial_at(h);
partition=struct('kind','direct_tableau_interval','fitted',false, ...
    'initialDt',h,'segmentDt',h,'knotX',NaN,'knotIndex',NaN,'direction',0, ...
    'secantIterations',0,'extraRhsEvaluations',0,'knotSlope',NaN,'rawPeakKnotDefect',NaN, ...
    'coordinateIdentityUlps',256,'startClassification',struct(),'endClassification',struct());
[out,records,h,partition]=ipm_cr.time.knot.fit(state,out,records,h,ctx.wall.x,@trial_at,partition);
partition.extraRhsEvaluations=calls-4;
[maximumRate,index]=max(records.rates);
[~,clock]=ipm_cr.time.clock.stage(state.scale,pack(out.scale),records.physicalClockIncrement);
order=4;stages=4;scheme='rk4_fixed4_spline_knot_partition_v4';
info=struct('method',scheme,'timeIntegrator','rk4', ...
    'temporalOrder',order,'spatialOrderCertified',false,'dt',h,'seconds',toc(started), ...
    'initialProposedDt',initialH,'timePartition',partition, ...
    'rate',state.rhs.rate,'maximumStageRate',maximumRate, ...
    'realizedMaximumStageCfl',h*maximumRate,'configuredCfl',policy.cfl, ...
    'cflStepFormulaApplied',true,'cflStepFormula','min(maxDt,0.95*cfl/currentRate,remainingTau)', ...
    'stageRates',records.rates,'stageCfl',h*records.rates, ...
    'stageAbscissae',records.abscissae,'stageCount',stages, ...
    'maximumCflStageIndex',index,'cflPosition',out.rhs.cflPosition, ...
    'maximumStageCflPosition',records.positions{index}, ...
    'stagePeakDiagnostics',records.peaks,'allAttemptPeakDiagnostics',peaks, ...
    'stageGaugeResidualInfinity',records.gaugeResidualInfinity, ...
    'maximumStageGaugeResidual',max(records.gaugeResidualInfinity), ...
    'stagePhysicalClockRates',records.clockRates,'rhsEvaluations',calls, ...
    'initialRhsReused',true,'terminalRhsEvaluated',true,'terminalRhsReusedNextStep',true, ...
    'cflScope','initial interval from current RHS; ordinary spline-knot partition; stage CFL recorded without rejection', ...
    'physicalClockIncrement',records.physicalClockIncrement,'physicalClockUpdate',clock, ...
    'compensatedPhysicalClock',true,'originAmplitudeAdvancedWithSameTableau',true, ...
    'sspMethod',false,'analyticBackgroundVariableChange',hasResidual);

    function [value,record]=trial_at(interval)
        probeH=interval;
        [value,record]=trial(state,interval,ctx.mesh,callbacks);
    end
    function [f,W]=stage(field,z,previous,wantDiagnostics,increment)
        calls=calls+1;
        staged=state;staged.scale=ipm_cr.time.clock.stage(state.scale,z,increment);staged.rhs=previous;
        if hasResidual,staged.eta=field;else,staged.W=field;end
        [f,W]=ipm_cr.rhs(staged,ctx,wantDiagnostics);
        note_peak(f,z(5),(z(5)-state.scale.canonicalTime)/probeH,false,true);
    end
    function s=finish(old,field,W,z,f,step,increment)
        s=old;s.W=W;s.rhs=f;
        s.scale=ipm_cr.time.clock.stage(old.scale,z,increment);
        s.scale.canonicalTime=old.scale.canonicalTime+step;
        if hasResidual,s.eta=field;end
        s.step=old.step+1;
    end
    function note_peak(f,tau,abscissa,reused,returned)
        p=ipm_cr.time.knot.peak(f);p.canonicalTime=tau;p.abscissa=abscissa;
        p.initialRhsReused=reused;p.rhsReturned=returned;
        if isempty(peaks),peaks=p;else,peaks(end+1)=p;end
    end
end

function z=pack(s)
z=[s.logC_l;s.logC_omega;s.physicalTime;s.X_shift;s.canonicalTime;s.densityOrigin];
end

function [s,records]=trial(initial,h,g,callbacks)
%TRIAL One complete explicit trial, with unchanged fixed-tableau arithmetic.
% No acceptance, retry, event search, or geometry ownership occurs here.
initialZ=[initial.scale.logC_l;initial.scale.logC_omega; ...
    initial.scale.physicalTime;initial.scale.X_shift; ...
    initial.scale.canonicalTime;initial.scale.densityOrigin];
hasResidual=isfield(initial,'eta');initialGrid=initial.W;
if hasResidual,initialGrid=initial.eta;end
stage=callbacks.stage;note_peak=callbacks.notePeak;finish=callbacks.finish;
tableau=rk4();
[s,records]=attempt_rk4(initial,h);

    function [s,records]=attempt_rk4(old,h)
        s=[];z=initialZ;f0=old.rhs;k0=scalar_rate(z,f0);
        note_peak(f0,z(5),0,true,true);
        records=record_start(f0,k0(3));r0=grid_rhs(f0);
        q1=add(initialGrid,1,r0,h*tableau.A(2,1));z1=z+h*tableau.A(2,1)*k0;
        f1=stage(q1,z1,f0,false,h*tableau.A(2,1)*k0(3));k1=scalar_rate(z1,f1);
        records=record_add(records,f1,k1(3),tableau.c(2));
        r1=grid_rhs(f1);
        q2=add(initialGrid,1,r1,h*tableau.A(3,2));z2=z+h*tableau.A(3,2)*k1;
        f2=stage(q2,z2,f1,false,h*tableau.A(3,2)*k1(3));k2=scalar_rate(z2,f2);
        records=record_add(records,f2,k2(3),tableau.c(3));
        r2=grid_rhs(f2);
        q3=add(initialGrid,1,r2,h*tableau.A(4,3));z3=z+h*tableau.A(4,3)*k2;
        f3=stage(q3,z3,f2,false,h*tableau.A(4,3)*k2(3));k3=scalar_rate(z3,f3);
        records=record_add(records,f3,k3(3),tableau.c(4));
        b=tableau.b;r3=grid_rhs(f3);
        qNew=struct('core',initialGrid.core+h*(b(1)*r0.core+b(2)*r1.core+ ...
            b(3)*r2.core+b(4)*r3.core), ...
            'polar',initialGrid.polar+h*(b(1)*r0.polar+b(2)*r1.polar+ ...
            b(3)*r2.polar+b(4)*r3.polar));
        zNew=z+h*(b(1)*k0+b(2)*k1+b(3)*k2+b(4)*k3);
        zNew(5)=z(5)+h;
        increment=h*(b(1)*k0(3)+b(2)*k1(3)+b(3)*k2(3)+b(4)*k3(3));
        qNew=ipm_cr.mesh.exchange(qNew,g);
        [fEnd,rhoActual]=stage(qNew,zNew,f3,true,increment);
        records=record_add(records,fEnd,exp(zNew(2)-zNew(1)),1);
        s=finish(old,qNew,rhoActual,zNew,fEnd,h,increment);records.physicalClockIncrement=increment;
    end
    function r=grid_rhs(f)
        if hasResidual,r=f.residualRhs;else,r=f.F;end
    end

end

function r=record_start(f,clockRate)
r=struct('rates',f.rate,'positions',{{f.cflPosition}},'abscissae',0, ...
    'peaks',ipm_cr.time.knot.peak(f), ...
    'gaugeResidualInfinity',max(abs(f.gaugeResidual)),'clockRates',clockRate,'physicalClockIncrement',NaN);
end
function r=record_add(r,f,clockRate,abscissa)
r.rates(end+1)=f.rate;r.positions{end+1}=f.cflPosition;
r.abscissae(end+1)=abscissa;r.gaugeResidualInfinity(end+1)=max(abs(f.gaugeResidual));
r.clockRates(end+1)=clockRate;
r.peaks(end+1)=ipm_cr.time.knot.peak(f);
end
function r=add(a,wa,b,wb)
r=struct('core',wa*a.core+wb*b.core,'polar',wa*a.polar+wb*b.polar);
end
function k=scalar_rate(z,f)
k=[f.c_l;f.c_omega;exp(z(2)-z(1));f.c_l*z(4)+f.c_r;1;f.transportOrigin+f.c_omega*z(6)];

end

function tableau=rk4()
%RK4TABLEAU Classical four-stage, fourth-order EXPLICIT temporal method.
% Production advance and transport tests share these actual coefficients.
% Its SSP coefficient is zero; a physical CFL check alone does not certify
% the stability or spatial accuracy of an overset spatial discretization.
persistent retained
if isempty(retained)
    retained=struct('A',[0,0,0,0;.5,0,0,0;0,.5,0,0;0,0,1,0], ...
        'b',[1/6,1/3,1/3,1/6],'c',[0;.5;.5;1], ...
        'order',4,'stages',4,'sspCoefficient',0,'implicit',false, ...
        'id','classical_explicit_rk4');
end
tableau=retained;
end
