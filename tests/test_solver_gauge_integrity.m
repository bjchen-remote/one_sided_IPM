function report=test_solver_gauge_integrity(outputDir)
%TEST_SOLVER_GAUGE_INTEGRITY Root regularity, unchanged IPC and split Jacobian.
% Frozen-state tests only: no trajectory or physical clock is advanced.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'results','solver_gauge_integrity');end
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'solver_gauge_integrity.json')), ...
    'ipm_cr:GaugeIntegrityOutput','Use a fresh output directory.');
timer=tic;identity=ipm_cr.io.identity('source');
report=struct('kind','ipm_cr_gauge_integrity_v1','passed',false,'timeAdvanced',0, ...
    'sourceIdentity',identity,'sourceUnchanged',false, ...
    'scope','Finite-axis original IPC; no continuum-order or trajectory certificate.');
try
    [R,r,ops]=wall_fixture();
    report.centeredPeak=centered_peak_test(R,r,ops);
    report.rootGates=root_gate_test(R,r,ops);
    report.finiteAxisQuadrature=quadrature_test();
    report.splitJacobian=split_jacobian_test();
    report.sourceUnchanged=isequaln(identity,ipm_cr.io.identity('source'));
    assert(report.sourceUnchanged,'ipm_cr:GaugeIntegritySource','Numerical source changed during testing.');
    report.passed=true;report.status='passed_solver_gauge_integrity';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'solver_gauge_integrity.json'),report);
    rethrow(problem);
end
report.seconds=toc(timer);ipm_cr.io.json(fullfile(outputDir,'solver_gauge_integrity.json'),report);
fprintf('CR gauge integrity: root gates, one centered peak and split Jacobian passed (%.3fs)\n',report.seconds);
end

function [R,r,ops]=wall_fixture()
r=(0:80)'/40;x=[-flipud(r(2:end));r]';X=[r,-r];
primitive=@(z)erf(z-.137)+.07*z;
R=primitive(X)-primitive(0);R(1,:)=0;
settings=struct('wallIntegralRadius',1,'referenceWallIntegral',17, ...
    'referenceWallGradientCentroid',.012,'peakPositionRestoringRate',.5, ...
    'cOmegaGaugeFloor',1e-12,'pointStrainConditionFloor',1e-12,'pinX',.05);
ops=struct('x',x,'Dx',ipm_cr.mesh.fd(x,1,7),'rescaling',settings);
end

function result=centered_peak_test(R,r,ops)
originalOps=ops;
[raw,rawObserve]=ipm_cr.gauge.measure(R,r,ops,struct('interpolation','linear'));
anchor=max(raw.q);expected=ipm_cr.gauge.peak(raw.q-anchor,ops);
expected.evaluationAnchor=anchor;
expected.centeredMaximumCandidateValue=expected.maximumCandidateValue;
expected.value=expected.value+anchor;expected.gridValue=expected.gridValue+anchor;
expected.maximumCandidateValue=expected.maximumCandidateValue+anchor;
expected.tieToleranceUsesCenteredValues=true;
% Profile the actual entry point, rather than trusting a self-reported count.
profile clear;profile on;
cleanup=onCleanup(@()profile('off'));
[m,observe]=ipm_cr.gauge.wall(R,r,ops,struct('interpolation','linear'));
profile off;details=profile('info');clear cleanup;
peakFile=fullfile(ipm_cr.root(),'src','+ipm_cr','+gauge','peak.m');
table=details.FunctionTable;calls=0;
for j=1:numel(table)
    if strcmp(table(j).FileName,peakFile)&&~contains(table(j).FunctionName,'>')
        calls=calls+table(j).NumCalls;
    end
end
assert(calls==1,'ipm_cr:GaugeIntegrityPeakWork','One wall measurement must compute the global peak once.');
assert(isequaln(m.peak,expected)&&m.P==expected.value&& ...
    isequaln(m.integral,raw.integral)&&isequaln(m.centroid,raw.centroid)&& ...
    isequaln(ops,originalOps),'ipm_cr:GaugeIntegrityObservable', ...
    'The centered peak, finite-axis IPC and immutable references must be preserved.');
X=[r,-r];F=sin(X).*exp(-X.^2);
bound=observe(F);old=rawObserve(F);
expectedRate=old.rate.normalized;
expectedRate(2)=expected.weights*old.rate.qDot(expected.indices)'/expected.value;
independent=ipm_cr.gauge.wall(R,r,ops,struct('interpolation','linear','rateField',F));
assert(isequaln(bound.rate.normalized,expectedRate)&& ...
    isequaln(bound.rate,independent.rate)&&isequaln(bound.peak,m.peak), ...
    'ipm_cr:GaugeIntegrityRate','Bound and independently rebuilt actual-F measurements must agree.');
result=struct('peakCallCount',calls,'legacyCenteredPeakBitwiseIdentical',true, ...
    'integralAndCentroidBitwiseUnchanged',true,'referencesUnchanged',true, ...
    'independentCompleteRateBitwiseIdentical',true,'peakX',m.peak.x,'peakValue',m.P, ...
    'referenceI',m.referenceWallIntegral,'referenceC',m.referenceWallGradientCentroid);
end

function result=root_gate_test(R,r,ops)
[m,observe]=ipm_cr.gauge.wall(R,r,ops);
ops.rescaling.referenceWallGradientCentroid=m.C;
X=[r,-r];templates={X,X.^2,X.^3,X.^4,exp(-X.^2),X.*R};
K=zeros(2,numel(templates));
for j=1:numel(templates)
    seen=observe(templates{j});v=seen.rate.normalized;K(:,j)=[v(1)-v(2);v(3)];
end
assert(rcond(K*K')>1e-8,'ipm_cr:GaugeIntegrityFixture','Independent IPC test directions are required.');
coefficients=K'/(K*K');P1=zeros(size(R));P2=P1;
for j=1:numel(templates)
    P1=P1+coefficients(j,1)*templates{j};P2=P2+coefficients(j,2)*templates{j};
end
kind='regular';
[~,regular]=ipm_cr.gauge.rates(R,r,ops,@base,0,[0,0],m,observe);
assert(regular.rootCondition>=regular.rootConditionFloor&&~regular.fallbackUsed, ...
    'ipm_cr:GaugeIntegrityRegular','A regular exact root should pass without a fallback.');
kind='zero';expect(@()ipm_cr.gauge.rates(R,r,ops,@base,0,[0,0],m,observe),'ipm_cr:ExplicitGaugeRoot');
kind='illconditioned';expect(@()ipm_cr.gauge.rates(R,r,ops,@base,0,[0,0],m,observe),'ipm_cr:ExplicitGaugeRoot');
kind='nonfinite';expect(@()ipm_cr.gauge.rates(R,r,ops,@base,0,[0,0],m,observe),'ipm_cr:ExplicitGaugeRoot');
% Exact integer first differences make q==1 bitwise; this isolates the
% genuinely flat maximum from a high-order stencil's last-bit noise.
tieR=(0:4)';flat=[tieR,-tieR];tieOps=ops;tieOps.x=-4:4;
tieOps.Dx=spdiags([-ones(9,1),ones(9,1)],[0,1],9,9);
tieOps.Dx(9,:)=0;tieOps.Dx(9,8:9)=[-1,1];
tieOps.rescaling.wallIntegralRadius=3;
[tied,tiedObserve]=ipm_cr.gauge.wall(flat,tieR,tieOps);
assert(tied.peak.tieNonunique&&~tied.peak.valueDerivativeUnique, ...
    'ipm_cr:GaugeIntegrityTie','A flat positive gradient must expose its nonunique maximum.');
expect(@()ipm_cr.gauge.rates(flat,tieR,tieOps,@base,0,[0,0],tied,tiedObserve),'ipm_cr:ExplicitGaugePeak');
result=struct('regularCondition',regular.rootCondition,'conditionFloor',regular.rootConditionFloor, ...
    'zeroResidualSingularRejected',true,'zeroResidualIllconditionedRejected',true, ...
    'nonfiniteJacobianRejected',true,'nonuniqueGlobalMaximumRejected',true, ...
    'globalMaximumDefinitionUnchanged',true);

    function [b,partials]=base(cl,delta)
        a=P1;d=P2;
        if strcmp(kind,'zero'),a(:)=0;d(:)=0;
        elseif strcmp(kind,'illconditioned'),d=1e-14*d;end
        b=cl*a+delta*d;
        if strcmp(kind,'nonfinite')
            % Every derivative-field entry is finite, but its wall derivative
            % overflows. A nonfinite observed speed Jacobian must not pass.
            a=(realmax/2)*repmat((-1).^(1:size(R,1))',1,size(R,2));a(1,:)=0;
        end
        partials=struct('c_l',a,'delta',d);
    end
end

function result=quadrature_test()
a=.25;b=.2;exactI=2+2*a/3;
energy=2+2*(b*b+2*a)/3+2*a*a/5;exactC=(4*b/3+4*a*b/5)/energy;
counts=[41,81,161];errors=zeros(numel(counts),2);values=errors;
for j=1:numel(counts)
    x=linspace(-2,2,counts(j));q=1+b*x+a*x.^2;
    settings=struct('wallIntegralRadius',1,'referenceWallIntegral',17, ...
        'referenceWallGradientCentroid',.012,'cOmegaGaugeFloor',1e-12,'pointStrainConditionFloor',1e-12);
    ops=struct('x',x,'rescaling',settings);
    integral=ipm_cr.gauge.integral(q,ops,settings);centroid=ipm_cr.gauge.centroid(q,ops);
    values(j,:)=[integral.value,centroid.value];errors(j,:)=abs(values(j,:)-[exactI,exactC]);
    assert(centroid.quadratureOrder==2&&centroid.referenceValue==settings.referenceWallGradientCentroid, ...
        'ipm_cr:GaugeIntegrityQuadrature','The original finite-axis quadrature and reference must be retained.');
end
orders=log2(errors(1:end-1,:)./errors(2:end,:));
assert(all(errors>0,'all')&&all(orders>1.9&orders<2.1,'all'), ...
    'ipm_cr:GaugeIntegrityContinuumBias','Finite-axis I/C must expose their measurable second-order continuum bias.');
result=struct('wallCounts',counts,'continuousI',exactI,'continuousC',exactC, ...
    'discreteValues',values,'absoluteContinuumBias',errors,'observedOrders',orders, ...
    'definition','Original finite-axis piecewise-linear I and positive trapezoidal C; weights are unchanged.');
end

function result=split_jacobian_test()
method=struct('geometry',struct('coreSpacing',.06, ...
    'grid',struct('nr',193,'ntheta',65,'firstSpacing',.06)), ...
    'solverOptions',struct('coreDegree',[48,48],'angularDegree',64,'polarDegree',16,'rOuter',1e10));
[s,ctx]=ipm_cr.init(method);g=ctx.mesh;originalReferences=s.originalAuthority;
assert(isfield(s,'eta')&&ctx.opts.background.tailCorrection,'ipm_cr:GaugeIntegritySplit');
for patch={'core','polar'}
    key=patch{1};p=g.(key);s.eta.(key)=.002*(p.X+.1*p.X.^2).*exp(-.5*(p.X.^2+p.Y.^2));
end
[f,W]=ipm_cr.rhs(s,ctx,true);eta=ipm_cr.mesh.fill(s.eta,g);
background=struct();forcing=struct();physical=f.physical;
for patch={'core','polar'}
    key=patch{1};p=g.(key);
    background.(key)=ipm_cr.model.background(ctx.opts.background,s.scale,p.X,p.Y,false);
    bg=background.(key);forcing.(key)=-physical.(key).u1.*bg.x-physical.(key).u2.*bg.y;
end
map=ctx.wallMap;r=ctx.wallR;R=sample_wall(W);R(1,:)=0;
[m,observe]=ipm_cr.gauge.wall(R,r,ctx.wall);
prepared=ipm_cr.op.prepare(eta,g,struct('transport',ctx.opts.transport));
Uref=f.Uref;point=[f.c_l;f.delta];
[~,partials]=base(point(1),point(2));J=zeros(2);
for j=1:2
    name='c_l';if j==2,name='delta';end
    value=observe(partials.(name));v=value.rate.normalized;J(:,j)=[v(1)-v(2);v(3)];
end
steps=[1e-4,1e-5,1e-6];errors=zeros(size(steps));wallErrors=zeros(numel(steps),2);
for j=1:numel(steps)
    fd=zeros(2);
    for k=1:2
        h=steps(j)*max(1,abs(point(k)));e=zeros(2,1);e(k)=h;
        plus=base(point(1)+e(1),point(2)+e(2));minus=base(point(1)-e(1),point(2)-e(2));
        vp=observe(plus);vm=observe(minus);
        p=vp.rate.normalized;v=vm.rate.normalized;
        fd(:,k)=([p(1)-p(2);p(3)]-[v(1)-v(2);v(3)])/(2*h);
        name='c_l';if k==2,name='delta';end
        wallErrors(j,k)=norm((plus-minus)/(2*h)-partials.(name),inf)/max(1,norm(partials.(name),inf));
    end
    errors(j)=norm(fd-J,inf)/max(1,norm(J,inf));
end
assert(all(errors(2:end)<2e-7)&&all(wallErrors(2:end,:)<2e-7,'all'), ...
    'ipm_cr:GaugeIntegritySplitJacobian','The production split speed Jacobian must agree with independent central differences.');
[rates,info]=ipm_cr.gauge.rates(R,r,ctx.wall,@base,Uref,[f.c_l,f.c_r],m,observe);
fresh=base(rates.c_l,rates.delta);fresh=fresh+rates.c_omega*R;
independent=ipm_cr.gauge.wall(R,r,ctx.wall,struct('rateField',fresh));
target=[0;0;ctx.wall.rescaling.peakPositionRestoringRate*(m.C-ctx.wall.rescaling.referenceWallGradientCentroid)];
audit=independent.rate.normalized+target;
assert(norm(audit,inf)<=1e-8&&isequaln(s.originalAuthority,originalReferences)&&s.scale.canonicalTime==0, ...
    'ipm_cr:GaugeIntegritySplitAudit','Frozen split rates must pass the independent complete-RHS IPC audit.');
result=struct('residualPresent',true,'tailCorrection',true,'relativeSteps',steps, ...
    'jacobian',J,'condition',rcond(J),'jacobianMixedErrors',errors,'wallPartialMixedErrors',wallErrors, ...
    'actualCompleteIPCRateResidual',audit,'rootInformation',info,'referencesUnchanged',true, ...
    'timeAdvanced',0,'scope','Fresh production split density/source/Poisson; independent split transport callback and central speed differences.');

    function wall=sample_wall(field)
        wall=reshape(map.core*field.core(:)+map.polar*field.polar(:),map.shape);
    end
    function [wall,partials]=base(cl,delta)
        vx=struct('core',(physical.core.u1-Uref)+cl*g.core.X+delta, ...
            'polar',(physical.polar.u1-Uref)+cl*g.polar.X+delta);
        vy=struct('core',physical.core.u2+cl*g.core.Y,'polar',physical.polar.u2+cl*g.polar.Y);
        applied=prepared(vx,vy,struct('cfl',false,'partials',nargout>1));
        T=ipm_cr.mesh.fill(struct('core',applied.B.core+forcing.core,'polar',applied.B.polar+forcing.polar),g);
        cr=delta-Uref;full=struct();
        for patch={'core','polar'}
            key=patch{1};p=g.(key);bg=background.(key);
            full.(key)=T.(key)-(cl*p.X+cr).*bg.x-cl*p.Y.*bg.y;
        end
        wall=sample_wall(full);wall(1,:)=wall(1,1);partials=[];
        if nargout>1
            partials=applied.partials;
            for parameter={'c_l','delta'}
                name=parameter{1};v=ipm_cr.mesh.fill(partials.(name),g);
                for patch={'core','polar'}
                    key=patch{1};p=g.(key);bg=background.(key);
                    if strcmp(name,'c_l'),v.(key)=v.(key)-p.X.*bg.x-p.Y.*bg.y;
                    else,v.(key)=v.(key)-bg.x;end
                end
                value=sample_wall(v);value(1,:)=value(1,1);partials.(name)=value;
            end
        end
    end
end

function expect(callback,identifier)
try
    callback();error('ipm_cr:GaugeIntegrityExpectedFailure','Expected %s.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'ipm_cr:GaugeIntegrityRejection', ...
        'Expected %s; received %s.',identifier,problem.identifier);
end
end
