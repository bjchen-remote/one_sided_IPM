function tableau=rk4Tableau()
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
