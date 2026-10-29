function [rates,info]=rates(R,r,wallOps,baseEvaluator,Uref,guess,m,observe)
%GAUGE Newton solve of the original I/P/C equations with exact speed Jacobian.
% The caller binds one immutable density and its centered natural-spline
% peak. The maintained solver may supply a cold seed; every accepted root
% is then solved and audited using the SAME supplied observables.
% A nonunique global maximum has no single P derivative. Preserve that
% observable and fail closed rather than silently selecting a smooth branch.
assert(isfield(m,'peak')&&isfield(m.peak,'valueDerivativeUnique')&& ...
    isfield(m.peak,'tieNonunique')&&m.peak.valueDerivativeUnique&&~m.peak.tieNonunique, ...
    'ipm_cr:ExplicitGaugePeak', ...
    'The global wall-gradient maximum is nonunique; its I/P/C rate closure has no unique peak derivative.');
target=[0;0;wallOps.rescaling.peakPositionRestoringRate* ...
    (m.C-wallOps.rescaling.referenceWallGradientCentroid)];
tolerance=16*eps;stagnationLimit=2e-12;
conditionFloor=1e-12;
evaluations=0;iterations=0;condition=NaN;stagnationAccepted=false;
fallbackUsed=false;fallbackInfo=struct();
[z,found,g,v,b]=newton([guess(1);guess(2)+Uref]);
if ~found
    [~,seed,fallbackInfo]=ipm_cr.gauge.seed(R,r,wallOps,baseEvaluator,Uref,guess, ...
        struct('measurementOptions',struct('interpolation','linear'),'rootTolerance',1e-10));
    fallbackUsed=true;
    [z,found,g,v,b]=newton([seed.c_l;seed.delta]);
end
assert(found,'ipm_cr:ExplicitGaugeRoot', ...
    ['The complete original I/P/C equations require a finite, sufficiently conditioned ' ...
    'speed Jacobian (rcond >= %.3g); last rcond %.3g.'],conditionFloor,condition);
cw=-v(2);complete=observe(b+cw*R);residual=complete.rate.normalized+target;
assert(all(isfinite(residual))&&norm(residual,inf)<=1e-8, ...
    'ipm_cr:ExplicitGaugeRoot','The assembled wall RHS failed I/P/C.');
rates=struct('c_l',z(1),'delta',z(2),'c_r',z(2)-Uref,'c_omega',cw,'Uref',Uref);
info=struct('definition','original_IPC_exact_transport_speed_Jacobian_Newton', ...
    'auditResidual',residual,'rootResidual',g,'rootTolerance',tolerance, ...
    'stagnationResidualLimit',stagnationLimit,'boundedStagnationAccepted',stagnationAccepted, ...
    'rootCondition',condition,'rootConditionFloor',conditionFloor, ...
    'iterations',iterations,'evaluations',evaluations, ...
    'branchJacobianUsed',true,'fallbackUsed',fallbackUsed,'fallbackUsedOnlyAsSeed',fallbackUsed, ...
    'fallbackInformation',fallbackInfo);

    function [point,found,g,v,b]=newton(point)
        found=false;
        for iteration=1:8
            iterations=iterations+1;
            [g,v,b,J]=probe(point);merit=norm(g,inf);
            if ~regular_jacobian(J),return;end
            if merit<=tolerance,found=true;return;end
            direction=-J\g;
            limit=4*max(1,norm(point,inf));
            if norm(direction,inf)>limit,direction=direction*(limit/norm(direction,inf));end
            accepted=false;
            for damping=2.^(-(0:7))
                candidate=point+damping*direction;trial=probe(candidate);
                if norm(trial,inf)<merit||norm(trial,inf)<=tolerance
                    point=candidate;accepted=true;break;
                end
            end
            if ~accepted
                % No further decrease in double precision. This preserves
                % the prior root bound, reports the exception explicitly,
                % and does not alter the independent full-RHS audit.
                if merit<=stagnationLimit,found=true;stagnationAccepted=true;end
                return;
            end
        end
        [g,v,b,J]=probe(point);
        found=regular_jacobian(J)&&norm(g,inf)<=tolerance;
    end
    function regular=regular_jacobian(J)
        condition=NaN;
        if all(isfinite(J),'all'),condition=rcond(J);end
        regular=isfinite(condition)&&condition>=conditionFloor;
    end
    function [g,v,b,J]=probe(point)
        evaluations=evaluations+1;
        if nargout>3,[b,partials]=baseEvaluator(point(1),point(2));
        else,b=baseEvaluator(point(1),point(2));end
        observed=observe(b);v=observed.rate.normalized;
        g=[v(1)-v(2);v(3)+target(3)];
        assert(all(isfinite(g)),'ipm_cr:ExplicitGaugeRoot','Nonfinite I/P/C residual.');
        if nargout>3
            cl=observe(partials.c_l);delta=observe(partials.delta);
            a=cl.rate.normalized;d=delta.rate.normalized;
            J=[a(1)-a(2),d(1)-d(2);a(3),d(3)];
        end
    end
end
