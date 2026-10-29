function [s,records]=trial(initial,h,method,g,options,callbacks)
%TRIAL One complete explicit trial, with unchanged fixed-tableau arithmetic.
% No acceptance, retry, event search, or geometry ownership occurs here.
initialZ=[initial.scale.logC_l;initial.scale.logC_omega; ...
    initial.scale.physicalTime;initial.scale.X_shift; ...
    initial.scale.canonicalTime;initial.scale.densityOrigin];
hasResidual=isfield(initial,'eta');initialGrid=initial.W;
if hasResidual,initialGrid=initial.eta;end
stage=callbacks.stage;note_peak=callbacks.notePeak;finish=callbacks.finish;
tableau=[];
if strcmp(method,'rk4')
    tableau=ipm_cr.time.rk4();[s,records]=attempt_rk4(initial,h);
else
    [s,records]=attempt_ssprk3(initial,h);
end
    function [s,records]=attempt_rk4(old,h)
        s=[];z=initialZ;f0=old.rhs;k0=scalar_rate(z,f0);
        note_peak(f0,z(5),0,true,true);
        records=record_start(f0,k0(3));r0=grid_rhs(f0);
        q1=add(initialGrid,1,r0,h*tableau.A(2,1));z1=z+h*tableau.A(2,1)*k0;
        f1=stage(q1,z1,f0,false,h*tableau.A(2,1)*k0(3));k1=scalar_rate(z1,f1);
        records=record_add(records,f1,k1(3),tableau.c(2));
        if ~inside_cfl(f1.rate,h),return;end
        r1=grid_rhs(f1);
        q2=add(initialGrid,1,r1,h*tableau.A(3,2));z2=z+h*tableau.A(3,2)*k1;
        f2=stage(q2,z2,f1,false,h*tableau.A(3,2)*k1(3));k2=scalar_rate(z2,f2);
        records=record_add(records,f2,k2(3),tableau.c(3));
        if ~inside_cfl(f2.rate,h),return;end
        r2=grid_rhs(f2);
        q3=add(initialGrid,1,r2,h*tableau.A(4,3));z3=z+h*tableau.A(4,3)*k2;
        f3=stage(q3,z3,f2,false,h*tableau.A(4,3)*k2(3));k3=scalar_rate(z3,f3);
        records=record_add(records,f3,k3(3),tableau.c(4));
        if ~inside_cfl(f3.rate,h),return;end
        b=tableau.b;r3=grid_rhs(f3);
        qNew=struct('core',initialGrid.core+h*(b(1)*r0.core+b(2)*r1.core+ ...
            b(3)*r2.core+b(4)*r3.core), ...
            'polar',initialGrid.polar+h*(b(1)*r0.polar+b(2)*r1.polar+ ...
            b(3)*r2.polar+b(4)*r3.polar));
        zNew=z+h*(b(1)*k0+b(2)*k1+b(3)*k2+b(4)*k3);
        zNew(5)=z(5)+h;
        increment=h*(b(1)*k0(3)+b(2)*k1(3)+b(3)*k2(3)+b(4)*k3(3));
        qNew=ipm_cr.mesh.fill(qNew,g);
        [fEnd,rhoActual]=stage(qNew,zNew,f3,true,increment);
        records=record_add(records,fEnd,exp(zNew(2)-zNew(1)),1);
        if ~inside_cfl(fEnd.rate,h),return;end
        s=finish(old,qNew,rhoActual,zNew,fEnd,h,increment);records.physicalClockIncrement=increment;
    end
    function [s,records]=attempt_ssprk3(old,h)
        % Preserve the historical SSPRK3 density/scalar arithmetic.
        s=[];z=initialZ;f0=old.rhs;k0=scalar_rate(z,f0);
        note_peak(f0,z(5),0,true,true);
        records=record_start(f0,k0(3));
        q1=add(initialGrid,1,grid_rhs(f0),h);z1=z+h*k0;
        q1=ipm_cr.mesh.fill(q1,g);
        f1=stage(q1,z1,f0,false,h*k0(3));k1=scalar_rate(z1,f1);
        records=record_add(records,f1,k1(3),1);
        if ~inside_cfl(f1.rate,h),return;end
        q2=add(initialGrid,.75,add(q1,1,grid_rhs(f1),h),.25);z2=.75*z+.25*(z1+h*k1);
        q2=ipm_cr.mesh.fill(q2,g);
        f2=stage(q2,z2,f1,false,h*(k0(3)+k1(3))/4);k2=scalar_rate(z2,f2);
        records=record_add(records,f2,k2(3),.5);
        if ~inside_cfl(f2.rate,h),return;end
        qNew=add(initialGrid,1/3,add(q2,1,grid_rhs(f2),h),2/3);zNew=(z+2*(z2+h*k2))/3;
        zNew(5)=z(5)+h;
        increment=h*(k0(3)/6+k1(3)/6+2*k2(3)/3);
        qNew=ipm_cr.mesh.fill(qNew,g);
        [fEnd,rhoActual]=stage(qNew,zNew,f2,true,increment);
        records=record_add(records,fEnd,exp(zNew(2)-zNew(1)),1);
        if ~inside_cfl(fEnd.rate,h),return;end
        s=finish(old,qNew,rhoActual,zNew,fEnd,h,increment);records.physicalClockIncrement=increment;
    end
    function r=grid_rhs(f)
        if hasResidual,r=f.residualRhs;else,r=f.F;end
    end
    function yes=inside_cfl(rate,h)
        yes=h*rate<=options.cfl*(1+64*eps);
    end
end

function r=record_start(f,clockRate)
r=struct('rates',f.rate,'positions',{{f.cflPosition}},'abscissae',0, ...
    'peaks',ipm_cr.time.peak(f), ...
    'gaugeResidualInfinity',max(abs(f.gaugeResidual)),'clockRates',clockRate,'physicalClockIncrement',NaN);
end
function r=record_add(r,f,clockRate,abscissa)
r.rates(end+1)=f.rate;r.positions{end+1}=f.cflPosition;
r.abscissae(end+1)=abscissa;r.gaugeResidualInfinity(end+1)=max(abs(f.gaugeResidual));
r.clockRates(end+1)=clockRate;
r.peaks(end+1)=ipm_cr.time.peak(f);
end
function r=add(a,wa,b,wb)
r=struct('core',wa*a.core+wb*b.core,'polar',wa*a.polar+wb*b.polar);
end
function k=scalar_rate(z,f)
k=[f.c_l;f.c_omega;exp(z(2)-z(1));f.c_l*z(4)+f.c_r;1;f.transportOrigin+f.c_omega*z(6)];
assert(all(isfinite(k))&&k(3)>0,'ipm_cr:ExplicitScale', ...
    'Every scalar rate must be finite, with a positive physical clock speed.');
end
