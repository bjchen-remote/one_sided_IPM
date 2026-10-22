function [flow,rho]=rhs(rho,g,wallOps,options,guess,background,wantDiagnostics,ctx)
%EVALUATE Complete coupled fifth-order transport RHS and original wall I/P/C semantic gauge.
% All ghost values are exchanged before differentiation/transport. Both
% interiors are independent explicit DOFs; neither overlap is projected.
if nargin<5,guess=[0,0];end
if nargin<6,background=[];end
if nargin<7,wantDiagnostics=true;end
validateattributes(wantDiagnostics,{'logical'},{'scalar'});
if nargin<8||isempty(ctx),ctx=ipm_cr.context(g,wallOps,options);end
split=~isempty(background);
if split
    [residual,sync]=ipm_cr.mesh.fill(background.residual,g);
    rho=struct('core',residual.core+background.core.H,'polar',residual.polar+background.polar.H);
    differentiated=residual;
else
    [rho,sync]=ipm_cr.mesh.fill(rho,g);differentiated=rho;
end
d=ipm_cr.mesh.diff(differentiated,g);
qNative=ipm_cr.mesh.fill(struct('core',d.core.x,'polar',d.polar.x),g);
if split
    q=struct('core',qNative.core+background.core.x,'polar',qNative.polar+background.polar.x);
    [physical,pinfo]=ctx.velocity(qNative,background.source,wantDiagnostics);
    forcing=struct('core',-physical.core.u1.*background.core.x-physical.core.u2.*background.core.y, ...
        'polar',-physical.polar.u1.*background.polar.x-physical.polar.u2.*background.polar.y);
else
    q=qNative;[physical,pinfo]=ctx.velocity(q,[],wantDiagnostics);
end
% The synthetic two-ray carrier represents exactly the ORIGINAL wall
% sample points. It is only an adapter for the unchanged I/P/C routines.
r=ctx.wallR;map=ctx.wallMap;
R=reshape(map.core*rho.core(:)+map.polar*rho.polar(:),map.shape);
R(1,:)=0;
[m,observeRate]=ipm_cr.gauge.wall(R,r,wallOps,struct('interpolation','linear'));
Uref=ipm_cr.mesh.sample(struct('core',physical.core.u1,'polar',physical.polar.u1),g,m.C,0);
prepare=ipm_cr.op.prepare(differentiated,g,struct('transport',options.transport));
[rates,ginfo]=ipm_cr.gauge.rates(R,r,wallOps,@base,Uref,guess,m,observeRate);
[B,rate,position,~,residualTransport,F0]=transport(rates.c_l,rates.delta);
if split
    T0=residualTransport.core(1,g.core.originIndex);
    % The prescribed background may also depend explicitly on physical
    % time. timeH is its CENTERED canonical-time contribution, including
    % dt/dtau. Subtract it only in eta's equation; the full-W gauge base
    % and the full origin-amplitude transport remain unchanged.
    residualRhs=struct('core',residualTransport.core+rates.c_omega*residual.core-T0-background.core.timeH, ...
        'polar',residualTransport.polar+rates.c_omega*residual.polar-T0-background.polar.timeH);
    residualRhs=ipm_cr.mesh.fill(residualRhs,g);
    residualRhs.core(1,g.core.originIndex)=0;
    F=struct();Q=[];
    if wantDiagnostics
        dF=ipm_cr.mesh.diff(residualRhs,g);
        dq=ipm_cr.mesh.fill(struct('core',dF.core.x,'polar',dF.polar.x),g);
        Q=struct();
    end
    for patch={'core','polar'}
        key=patch{1};p=g.(key);b=background.(key);
        gaugeX=rates.c_l*p.X+rates.c_r;gaugeY=rates.c_l*p.Y;
        Ht=-gaugeX.*b.x-gaugeY.*b.y+rates.c_omega*b.H+rates.c_r*b.originX+b.timeH;
        F.(key)=residualRhs.(key)+Ht;
        if wantDiagnostics
            Htx=(rates.c_omega-rates.c_l)*b.x-gaugeX.*b.xx-gaugeY.*b.xy+b.timeX;
            Q.(key)=dq.(key)+Htx;
        end
    end
else
    F=struct('core',B.core+rates.c_omega*rho.core-F0, ...
        'polar',B.polar+rates.c_omega*rho.polar-F0);
    F=ipm_cr.mesh.fill(F,g);
    Q=[];
    if wantDiagnostics
        dF=ipm_cr.mesh.diff(F,g);
        Q=struct('core',dF.core.x,'polar',dF.polar.x);
    end
end
F.core(1,g.core.originIndex)=0;
wallF=reshape(map.core*F.core(:)+map.polar*F.polar(:),map.shape);
wallF(1,:)=0;
if wantDiagnostics
    % Accepted endpoints rebuild the density observables independently.
    complete=ipm_cr.gauge.wall(R,r,wallOps, ...
        struct('interpolation','linear','rateField',wallF));
else
    % The density is immutable within this RHS. Its original I/P/C
    % functionals can be reused for an independently assembled actual F;
    % no observable or field is cached across states.
    complete=observeRate(wallF);
end
target=[0;0;wallOps.rescaling.peakPositionRestoringRate* ...
    (m.C-wallOps.rescaling.referenceWallGradientCentroid)];
actualResidual=complete.rate.normalized+target;
assert(all(isfinite(actualResidual))&&max(abs(actualResidual))<=1e-5, ...
    'ipm_cr:ExplicitGaugeAudit','The actual synchronized complete RHS must satisfy original I/P/C.');
flow=struct('F',F,'Q',Q,'q',q,'c_l',rates.c_l,'c_omega',rates.c_omega, ...
    'c_r',rates.c_r,'delta',rates.delta,'Uref',Uref,'transportOrigin',F0, ...
    'rate',rate+abs(rates.c_omega),'transportRate',rate,'cflPosition',position, ...
    'gaugeResidual',actualResidual,'gaugeSolverResidual',ginfo.auditResidual, ...
    'gaugeSolverInformation',ginfo, ...
    'measurement',m,'poisson',pinfo, ...
    'sourceSynchronization',sync,'physical',physical);
flow.diagnosticsComplete=wantDiagnostics;
if split
    flow.residualRhs=residualRhs;
    flow.densityRepresentation='W=eta+G-G(0); analytic prescribed background with canonical time derivatives';
    flow.analyticBackgroundSource=true;
end

    function [wall,partials]=base(cl,delta)
        if nargout>1,[b,~,~,partials]=transport(cl,delta);
        else,b=transport(cl,delta);end
        if ~split,b=ipm_cr.mesh.fill(b,g);end
        wall=reshape(map.core*b.core(:)+map.polar*b.polar(:),map.shape);
        wall(1,:)=wall(1,1);
        if nargout>1
            for field={'c_l','delta'}
                key=field{1};v=partials.(key);
                if ~split,v=ipm_cr.mesh.fill(v,g);end
                partials.(key)=reshape(map.core*v.core(:)+map.polar*v.polar(:),map.shape);
                partials.(key)(1,:)=partials.(key)(1,1);
            end
        end
    end
    function [b,rate,position,partials,T,originTransport]=transport(cl,delta)
        vx=struct('core',(physical.core.u1-Uref)+cl*g.core.X+delta, ...
            'polar',(physical.polar.u1-Uref)+cl*g.polar.X+delta);
        vy=struct('core',physical.core.u2+cl*g.core.Y,'polar',physical.polar.u2+cl*g.polar.Y);
        partials=[];T=[];
        request=struct('cfl',nargout>1&&nargout~=4,'partials',nargout==4);
        result=prepare(vx,vy,request);
        b=result.B;rate=result.rate;position=result.position;partials=result.partials;
        originTransport=b.core(1,g.core.originIndex);
        if split
            originTransport=originTransport-vx.core(1,g.core.originIndex)*background.core.originX;
            T=ipm_cr.mesh.fill(struct('core',b.core+forcing.core, ...
                'polar',b.polar+forcing.polar),g);
            cr=delta-Uref;
            for patch={'core','polar'}
                key=patch{1};p=g.(key);bg=background.(key);
                b.(key)=T.(key)-(cl*p.X+cr).*bg.x-cl*p.Y.*bg.y;
            end
            if nargout==4
                for field={'c_l','delta'}
                    key=field{1};v=ipm_cr.mesh.fill(partials.(key),g);
                    for patch={'core','polar'}
                        name=patch{1};p=g.(name);bg=background.(name);
                        if strcmp(key,'c_l'),v.(name)=v.(name)-p.X.*bg.x-p.Y.*bg.y;
                        else,v.(name)=v.(name)-bg.x;end
                    end
                    partials.(key)=v;
                end
            end
        end
    end
end
