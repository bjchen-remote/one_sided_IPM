function [rates,info]=rates(R,r,wallOps,baseEvaluator,Uref,guess,m,observe) %#ok<INUSD>
%RATES Direct Newton solve of the original I/P/C speed equations.
% The chosen global spline maximizer supplies its derivative, also at a tie.
% Residuals and conditioning are observations, not acceptance criteria.
target=[0;0;wallOps.rescaling.peakPositionRestoringRate* ...
    (m.C-wallOps.rescaling.referenceWallGradientCentroid)];
evaluations=0;iterations=0;
z=[guess(1);guess(2)+Uref];
[g,v,b,J]=probe(z);
for iteration=1:8
    iterations=iteration;
    point=z-J\g;
    if isequal(point,z),break;end
    z=point;[g,v,b,J]=probe(z);
end
condition=rcond(J);
cw=-v(2);complete=observe(b+cw*R);residual=complete.rate.normalized+target;
rates=struct('c_l',z(1),'delta',z(2),'c_r',z(2)-Uref,'c_omega',cw,'Uref',Uref);
info=struct('definition','original_IPC_exact_transport_speed_Jacobian_Newton', ...
    'auditResidual',residual,'rootResidual',g,'rootCondition',condition, ...
    'iterations',iterations,'evaluations',evaluations,'branchJacobianUsed',true);

    function [g,v,b,J]=probe(point)
        evaluations=evaluations+1;
        [b,partials]=baseEvaluator(point(1),point(2));
        observed=observe(b);v=observed.rate.normalized;
        g=[v(1)-v(2);v(3)+target(3)];
        cl=observe(partials.c_l);delta=observe(partials.delta);
        a=cl.rate.normalized;d=delta.rate.normalized;
        J=[a(1)-a(2),d(1)-d(2);a(3),d(3)];
    end
end
