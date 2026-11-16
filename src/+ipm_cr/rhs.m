function [evaluation,W]=rhs(s,ctx,wantDiagnostics)
%RHS Complete same-state RHS on the context's immutable mathematical setup.
% The optional second output is full centered W reconstructed by this RHS.
% Analytic background timeH/timeX already include the physical clock speed;
% only eta ghosts are interpolated in the background-subtracted branch.
if nargin<3,wantDiagnostics=true;end
assert(isfield(s,'ctxId')&&strcmp(s.ctxId,ctx.id),'ipm_cr:StaleState','State and context differ.');
guess=[s.rhs.c_l,s.rhs.c_r];
g=ctx.mesh;wallOps=ctx.wall;options=ctx.opts;background=[];
if isfield(s,'eta')
    model=options.background;scale=s.scale;
    background=struct('eta',s.eta,'scale',scale, ...
        'core',ipm_cr.model.background(model,scale,g.core.X,g.core.Y,wantDiagnostics), ...
        'polar',ipm_cr.model.background(model,scale,g.polar.X,g.polar.Y,wantDiagnostics));
    W=[];
else
    W=s.W;
end
validateattributes(wantDiagnostics,{'logical'},{'scalar'});
split=~isempty(background);
if split
    [eta,sync]=ipm_cr.mesh.exchange(background.eta,g);
    W=struct('core',eta.core+background.core.H,'polar',eta.polar+background.polar.H);
    differentiated=eta;
else
    [W,sync]=ipm_cr.mesh.exchange(W,g);differentiated=W;
end
d=ipm_cr.space.basis.gradient(differentiated,g,ctx.spline);
qNative=ipm_cr.mesh.exchange(struct('core',d.core.x,'polar',d.polar.x),g);
if split
    q=struct('core',qNative.core+background.core.x,'polar',qNative.polar+background.polar.x);
    [physical,pinfo]=ctx.velocity(q,wantDiagnostics);
    forcing=struct('core',-physical.core.u1.*background.core.x-physical.core.u2.*background.core.y, ...
        'polar',-physical.polar.u1.*background.polar.x-physical.polar.u2.*background.polar.y);
else
    q=qNative;[physical,pinfo]=ctx.velocity(q,wantDiagnostics);
end
% The synthetic two-ray carrier represents exactly the ORIGINAL wall
% sample points. It is only an adapter for the unchanged I/P/C routines.
r=ctx.wallR;map=ctx.wallMap;
R=ipm_cr.mesh.interpolate.sample(W,map);
R(1,:)=0;
[m,observeRate]=ipm_cr.gauge.observe.wall(R,r,wallOps,struct('interpolation','linear'));
Uref=ipm_cr.mesh.interpolate.sample(struct('core',physical.core.u1,'polar',physical.polar.u1),g,m.C,0);
wallPlan=ctx.wallProbe;
[prepare,probe]=ipm_cr.space.transport(differentiated,g,ctx.spline,wallPlan.nodes);
keys={'core','polar'};probeValues=struct();
for patch=keys
    key=patch{1};nodes=wallPlan.nodes.(key);sites=wallPlan.sites.(key);
    p=g.(key);v=physical.(key);
    probeValues.(key)=struct('X',p.X(nodes),'Y',p.Y(nodes), ...
        'u',v.u1(nodes)-Uref,'v',v.u2(nodes));
    if split
        bg=background.(key);
        probeValues.(key).forcing=forcing.(key)(nodes);
        probeValues.(key).siteX=p.X(sites);probeValues.(key).siteY=p.Y(sites);
        probeValues.(key).Hx=bg.x(sites);probeValues.(key).Hy=bg.y(sites);
    end
end
[rates,ginfo]=ipm_cr.gauge.rates(R,r,wallOps,@base,Uref,guess,m,observeRate);
[B,rate,position,residualTransport,F0]=transport(rates.c_l,rates.delta);
if split
    T0=residualTransport.core(1,g.core.originIndex);
    % The prescribed background may also depend explicitly on physical
    % time. timeH is its CENTERED canonical-time contribution, including
    % dt/dtau. Subtract it only in eta's equation; the full-W gauge base
    % and the full origin-amplitude transport remain unchanged.
    residualBase=struct('core',residualTransport.core+rates.c_omega*eta.core-T0, ...
        'polar',residualTransport.polar+rates.c_omega*eta.polar-T0);
    residualRhs=struct('core',residualBase.core-background.core.timeH, ...
        'polar',residualBase.polar-background.polar.timeH);
    residualRhs=ipm_cr.mesh.exchange(residualRhs,g);
    residualRhs.core(1,g.core.originIndex)=0;
    F=struct();Q=[];assemblyDifference=struct();
    if wantDiagnostics
        dF=ipm_cr.space.basis.gradient(residualRhs,g,ctx.spline);
        dq=ipm_cr.mesh.exchange(struct('core',dF.core.x,'polar',dF.polar.x),g);
        Q=struct();
    end
    for patch={'core','polar'}
        key=patch{1};p=g.(key);b=background.(key);
        gaugeX=rates.c_l*p.X+rates.c_r;gaugeY=rates.c_l*p.Y;
        Ht=-gaugeX.*b.x-gaugeY.*b.y+rates.c_omega*b.H+rates.c_r*b.originX;
        % Keep the prescribed ghost representation, where interpolation and
        % analytic background evaluation need not commute. At ACTIVE points
        % cancel the explicit timeH algebraically before floating arithmetic.
        F.(key)=residualRhs.(key)+(Ht+b.timeH);
        stable=residualBase.(key)+Ht;active=p.active;
        if wantDiagnostics
            assemblyDifference.(key)=max(abs(F.(key)(active)-stable(active)));
        end
        F.(key)(active)=stable(active);
        if wantDiagnostics
            Htx=(rates.c_omega-rates.c_l)*b.x-gaugeX.*b.xx-gaugeY.*b.xy+b.timeX;
            Q.(key)=dq.(key)+Htx;
        end
    end
else
    F=struct('core',B.core+rates.c_omega*W.core-F0, ...
        'polar',B.polar+rates.c_omega*W.polar-F0);
    F=ipm_cr.mesh.exchange(F,g);
    Q=[];
    if wantDiagnostics
        dF=ipm_cr.space.basis.gradient(F,g,ctx.spline);
        Q=struct('core',dF.core.x,'polar',dF.polar.x);
        if isfield(g.core,'stair'),Q=ipm_cr.mesh.exchange(Q,g);end
    end
end
F.core(1,g.core.originIndex)=0;
wallF=ipm_cr.mesh.interpolate.sample(F,map);
wallF(1,:)=0;
if wantDiagnostics
    % Accepted endpoints rebuild the density observables independently.
    complete=ipm_cr.gauge.observe.wall(R,r,wallOps, ...
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
evaluation=struct('F',F,'Q',Q,'q',q,'c_l',rates.c_l,'c_omega',rates.c_omega, ...
    'c_r',rates.c_r,'delta',rates.delta,'Uref',Uref,'transportOrigin',F0, ...
    'rate',rate+abs(rates.c_omega),'transportRate',rate,'cflPosition',position, ...
    'gaugeResidual',actualResidual,'gaugeSolverResidual',ginfo.auditResidual, ...
    'gaugeSolverInformation',ginfo, ...
    'measurement',m,'poisson',pinfo, ...
    'sourceSynchronization',sync,'physical',physical);
evaluation.diagnosticsComplete=wantDiagnostics;
if split
    evaluation.residualRhs=residualRhs;
    evaluation.densityRepresentation='W=eta+G-G(0); analytic prescribed background with canonical time derivatives';
    evaluation.analyticBackgroundSource=true;
    evaluation.completeFAssembly='explicit background time terms eliminated on active nodes; original ghost representation retained';
    if wantDiagnostics
        % q=D_h eta+H_X, so q_tau and D_h(W_tau) need not agree exactly.
        % This is a discrete representation defect, not continuum error.
        dComplete=ipm_cr.space.basis.gradient(F,g,ctx.spline);defect=struct();
        for patch={'core','polar'}
            key=patch{1};active=g.(key).active;
            delta=Q.(key)(active)-dComplete.(key).x(active);
            defect.(key)=struct('absoluteInfinity',max(abs(delta)), ...
                'activeSampleCount',nnz(active));
        end
        evaluation.representationAudit=struct('available',true,'QDefinition','time derivative of q=D_h eta+H_X', ...
            'comparison','Q-D_h F on active nodes','patches',defect, ...
            'maximumAbsoluteDefect',max(defect.core.absoluteInfinity,defect.polar.absoluteInfinity), ...
            'eliminatedTimeCancellationDifference',assemblyDifference, ...
            'continuumErrorCertified',false);
    end
end

    function [wall,partials]=base(cl,delta)
        for patch=keys
            key=patch{1};v=probeValues.(key);
            vx.(key)=v.u+cl*v.X+delta;vy.(key)=v.v+cl*v.Y;
        end
        [b,partials]=probe(vx,vy);
        if split
            for patch=keys
                key=patch{1};b.(key)=b.(key)+probeValues.(key).forcing;
            end
        end
        b=ipm_cr.mesh.interpolate.wallExchange(b,wallPlan);
        partials.c_l=ipm_cr.mesh.interpolate.wallExchange(partials.c_l,wallPlan);
        partials.delta=ipm_cr.mesh.interpolate.wallExchange(partials.delta,wallPlan);
        if split
            cr=delta-Uref;
            for patch=keys
                key=patch{1};v=probeValues.(key);
                b.(key)=b.(key)-(cl*v.siteX+cr).*v.Hx-cl*v.siteY.*v.Hy;
                partials.c_l.(key)=partials.c_l.(key)-v.siteX.*v.Hx-v.siteY.*v.Hy;
                partials.delta.(key)=partials.delta.(key)-v.Hx;
            end
        end
        wall=sample_probe(b);
        partials.c_l=sample_probe(partials.c_l);partials.delta=sample_probe(partials.delta);
    end
    function values=sample_probe(field)
        values=ipm_cr.mesh.interpolate.sample(field,wallPlan.map);
        values(1,:)=values(1,1);
    end
    function [b,rate,position,T,originTransport]=transport(cl,delta)
        vx=struct('core',(physical.core.u1-Uref)+cl*g.core.X+delta, ...
            'polar',(physical.polar.u1-Uref)+cl*g.polar.X+delta);
        vy=struct('core',physical.core.u2+cl*g.core.Y,'polar',physical.polar.u2+cl*g.polar.Y);
        T=[];
        result=prepare(vx,vy);
        b=result.B;rate=result.rate;position=result.position;
        originTransport=b.core(1,g.core.originIndex);
        if split
            originTransport=originTransport-vx.core(1,g.core.originIndex)*background.core.originX;
            T=ipm_cr.mesh.exchange(struct('core',b.core+forcing.core, ...
                'polar',b.polar+forcing.polar),g);
            cr=delta-Uref;
            for patch=keys
                key=patch{1};p=g.(key);bg=background.(key);
                b.(key)=T.(key)-(cl*p.X+cr).*bg.x-cl*p.Y.*bg.y;
            end
        end
    end
end
