function report=test_ipm_cr_tail_background(outputDir,options)
%TEST_IPM_CR_TAIL_BACKGROUND Independent dipole-tail representation checks.
% No trajectory is advanced. The Green-integral test concerns the initial
% physical density time derivative, not a certificate for a finite-time IVP.
if nargin<2,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','tail_background');end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options), ...
    {'includeInitial','initialOptions'})),'ipm_cr:TailTestOptions');
includeInitial=true;if isfield(options,'includeInitial'),includeInitial=options.includeInitial;end
initialOptions=struct();if isfield(options,'initialOptions'),initialOptions=options.initialOptions;end
validateattributes(includeInitial,{'logical'},{'scalar'});
outputDir=char(outputDir);assert(startsWith(outputDir,filesep),'ipm_cr:TailTestOutput');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'tail_background.json')),'ipm_cr:TailTestPreserve','Choose a fresh output directory.');
cfg=ipm_cr.problem();p=cfg.physics;
model=struct('kind','gaussian_primitive','amplitude',p.gaussianAmplitude, ...
    'centerX',p.gaussianCenterX,'widthX',p.gaussianWidthX, ...
    'widthY',p.gaussianWidthY,'skewWeight',p.gaussianSkewWeight,'tailCorrection',true);
timer=tic;report=struct('kind','independent_vertical_dipole_tail_background_checks_v1', ...
    'status','running','model',model,'timeAdvanced',0,'sourceInventory',source_inventory(root), ...
    'sourcesUnchanged',false,'allRequestedChecksPassed',false,'initial',struct(), ...
    'limitations',{{'The dipole term retains the exact far-wall X^-3 coefficient; it is an enrichment, not the full initial transport RHS in the interior.', ...
    'The Green-integral radius sweep tests the initial physical time derivative. Full finite-time spatial order remains the separate unmasked IVP test.', ...
    'Centered eta may contain a spatially constant component fixed by the origin convention; this test does not remove or fit a numerical tail.'}});
persist(outputDir,report);
try
    report.momentAndFarTail=far_tail_test(model);
    report.derivatives=derivative_test(model);
    report.scaledArithmetic=scaled_test(model);
    report.backgroundAssembly=assembly_test(model);
    report.kinematic=kinematic_test(model);
    if includeInitial,report.initial=initial_test(initialOptions);end
    report.sourcesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
    assert(report.sourcesUnchanged,'ipm_cr:TailTestSources');
    report.allRequestedChecksPassed=report.momentAndFarTail.passed&&report.derivatives.passed&& ...
        report.scaledArithmetic.passed&&report.backgroundAssembly.passed&&report.kinematic.passed&& ...
        (~includeInitial||report.initial.passed);
    if report.allRequestedChecksPassed,report.status='passed_requested_tail_checks';
    else,report.status='completed_measurements_with_failed_checks';end
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function result=far_tail_test(m)
% Independent half-plane Green differentiation, without subtracting nearly
% equal image kernels. Source quadrature covers the complete Gaussian up to
% explicitly reported exponentially small integration tails.
[xs,wx]=gauss(128,m.centerX-10*m.widthX,m.centerX+13.3*m.widthX);
[ys,wy]=gauss(96,0,8*m.widthY);[SX,SY]=ndgrid(xs,ys);
[~,source]=primitive(m,SX,SY);weights=wx*wy.';weighted=weights.*source;
My=sum(weighted.*SY,'all');exactMy=m.amplitude*m.widthX*sqrt(pi)*m.widthY^2/2;
C=My/pi;exactC=m.amplitude*m.widthX*m.widthY^2/(2*sqrt(pi));
s=scale0();probeX=[2,4];probeY=[.5,.8];d=ipm_cr.model.tail(m,s,probeX,probeY);
rho=primitive(m,probeX,probeY);D=1+probeX.^2+probeY.^2;
inferredC=-d.H.*(m.widthY^2*D.^2)./(4*probeX.*probeY.^2.*rho);
X=2.^(4:10);Y=[.3,.8,1.3];rows=struct([]);allPassed=true;
for j=1:numel(Y)
    y=Y(j);u1=zeros(size(X));u2=u1;
    for k=1:numel(X)
        dx=X(k)-SX;plus=dx.^2+(y+SY).^2;minus=dx.^2+(y-SY).^2;
        u1(k)=-sum(weighted.*SY.*(dx.^2+SY.^2-y^2)./(plus.*minus),'all')/pi;
        u2(k)=-2*y*sum(weighted.*SY.*dx./(plus.*minus),'all')/pi;
    end
    [~,qx,qy]=primitive(m,X,y+zeros(size(X)));
    physical=-u1.*qx-u2.*qy;
    tail=ipm_cr.model.tail(m,s,X,y+zeros(size(X)),false);
    difference=physical-tail.H;
    pRaw=log2(abs(physical(1:end-1)./physical(2:end)));
    pResidual=log2(abs(difference(1:end-1)./difference(2:end)));
    plateau=m.amplitude*m.widthX*sqrt(pi);
    leading=-4*C*plateau*y^2/m.widthY^2*exp(-(y/m.widthY)^2);
    signCorrect=all(physical<0)&&all(tail.H<0)&&leading<0;
    passed=signCorrect&&all(pRaw(end-1:end)>2.95)&&all(pResidual(end-1:end)>3.85);
    rows=append(rows,struct('Y',y,'X',X,'greenVelocity1',u1,'greenVelocity2',u2, ...
        'physicalInitialRhoTimeDerivative',physical,'dipoleVerticalEnrichment',tail.H, ...
        'physicalMinusEnrichment',difference,'rawDecayOrders',pRaw,'remainderDecayOrders',pResidual, ...
        'exactLeadingXCubeCoefficient',leading,'measuredXCubeCoefficients',X.^3.*physical, ...
        'allRadiusSamplesRetained',true,'signCorrect',signCorrect,'passed',passed));
    allPassed=allPassed&&passed;
end
left=ipm_cr.model.tail(m,s,-X,.8+zeros(size(X)),false);
coefficientError=metric(inferredC,exactC+zeros(size(inferredC)),2e-12);
result=struct('definition','C=M_y/pi; M_y=integral_{y>=0} y*rho0_x dx dy', ...
    'quadratureOrder',[128,96],'sourceIntegrationBox',[xs(1),xs(end),ys(1),ys(end)], ...
    'quadratureMoment',My,'analyticMoment',exactMy,'quadratureC',C,'analyticC',exactC, ...
    'coefficientInferredFromProvider',inferredC,'providerCoefficientError',coefficientError, ...
    'momentRelativeError',abs(My-exactMy)/exactMy,'records',rows,'negativeX',-X, ...
    'negativeXEnrichment',left.H,'exteriorMaskApplied',false, ...
    'scope','Uncentered physical initial time derivative. No fitted offset, numerical field alignment or tail deletion.', ...
    'passed',allPassed&&coefficientError.passed&&abs(My-exactMy)<=2e-12*exactMy);
end

function result=derivative_test(m)
s=scale0();X=[-2,-.3,.1,.7,2,4];Y=[0,.1,.8,1.2,.4,2];
d=ipm_cr.model.tail(m,s,X,Y);h=5e-4;step=1e-3;
ax=fivepoint(@(x)vertical_A(m,x,Y),X,h);ay=fivepoint(@(y)vertical_A(m,X,y),Y,h);
axx=second_fivepoint(@(x)vertical_A(m,x,Y),X,step);
axy=mixed_fivepoint(@(x,y)vertical_A(m,x,y),X,Y,step);
ref=struct('H',vertical_A(m,X,Y),'x',ax,'y',ay,'xx',axx,'xy',axy);
errors=struct();passed=true;
for key={'H','x','y','xx','xy'}
    name=key{1};tol=2e-10;if ismember(name,{'xx','xy'}),tol=2e-7;end
    errors.(name)=metric(d.(name),ref.(name),tol);passed=passed&&errors.(name).passed;
end
light=ipm_cr.model.tail(m,s,X,Y,false);
sourceOnly=ipm_cr.model.tail(m,s,X,Y,false,true);
sourceOnlyIdentical=isequaln(sourceOnly.x,d.x);
sourceOnlyHasOnlyX=isequal(fieldnames(sourceOnly),{'x'});
same=isequaln(light.H,d.H)&&isequaln(light.x,d.x)&&isequaln(light.y,d.y)&& ...
    isequaln(light.timeH,d.timeH)&&isequaln(light.timeX,d.timeX);
originExact=d.originValue==0&&d.originX==0&&all(d.H(Y==0)==0)&&all(d.x(Y==0)==0);
sharedReuse=shared_test(m,s,X,Y);
result=struct('X',X,'Y',Y,'firstDerivativeFivePointStep',h,'secondDerivativeFivePointStep',step, ...
    'reference','Independent elementary vertical A; fourth-order five-point differences for first/second derivatives and their tensor product for xy. No production derivative formula is used.', ...
    'errors',errors,'originAndWallIdentitiesExact',originExact,'lightFirstFieldsIdentical',same, ...
    'sourceOnlyXBitwiseIdentical',sourceOnlyIdentical,'sourceOnlyContainsOnlyX',sourceOnlyHasOnlyX, ...
    'sharedReuse',sharedReuse, ...
    'lightSecondDerivativesEmpty',isempty(light.xx)&&isempty(light.xy), ...
    'passed',passed&&originExact&&same&&sourceOnlyIdentical&&sourceOnlyHasOnlyX&& ...
        sharedReuse.passed&&isempty(light.xx)&&isempty(light.xy));
end

function result=scaled_test(m)
cases=struct('name',{'late_scale','late_scale_separate_y_square_underflow','combined_exponent_extreme'}, ...
    'ell',{56,56,400},'omega',{21,21,600},'shift',{6e24,6e24,.35*exp(400)}, ...
    'X',{[0,1,-1,2],[0,1,-1],[0,1,-1]},'Y',{[0,.3,1,2],[1e-139,2e-139,3e-139],[.3,1,2]});
rows=struct([]);passed=true;
for j=1:numel(cases)
    c=cases(j);s=scale0();s.logC_l=c.ell;s.logC_omega=c.omega;s.X_shift=c.shift;
    d=ipm_cr.model.tail(m,s,c.X,c.Y);
    sourceOnly=ipm_cr.model.tail(m,s,c.X,c.Y,false,true);
    sourceOnlyIdentical=isequaln(sourceOnly.x,d.x);
    sourceOnlyHasOnlyX=isequal(fieldnames(sourceOnly),{'x'});
    sharedReuse=shared_test(m,s,c.X,c.Y);
    x0=-c.shift*exp(-c.ell);f=@(x)vertical_coefficient(m,x);
    F=f(x0);Fx=fivepoint(f,x0,1e-4);Fxx=second_fivepoint(f,x0,1e-3);
    % y_phys is tiny in these cases. The omitted O(y_phys^2) relative term
    % is below 1e-47 (or underflows); each combined exponent stays intact.
    ref=struct('H',(exp(c.omega-2*c.ell)*c.Y.^2)*F, ...
        'x',(exp(c.omega-3*c.ell)*c.Y.^2)*Fx, ...
        'y',(2*exp(c.omega-2*c.ell)*c.Y)*F, ...
        'xx',(exp(c.omega-4*c.ell)*c.Y.^2)*Fxx, ...
        'xy',(2*exp(c.omega-3*c.ell)*c.Y)*Fx, ...
        'timeH',(exp(2*c.omega-3*c.ell)*c.Y.^2)*F, ...
        'timeX',(exp(2*c.omega-4*c.ell)*c.Y.^2)*Fx);
    errors=struct();ok=true;
    for key=fieldnames(ref)'
        name=key{1};tol=2e-10;if strcmp(name,'xx'),tol=2e-7;end
        errors.(name)=metric(d.(name),ref.(name),tol);ok=ok&&errors.(name).passed;
    end
    separate=exp(c.omega)*vertical_A(m,(c.X-c.shift)*exp(-c.ell),c.Y*exp(-c.ell));
    rows=append(rows,struct('name',c.name,'scale',s,'X',c.X,'Y',c.Y,'errors',errors, ...
        'computedH',d.H,'combinedExponentReferenceH',ref.H,'separatelyScaledDiagnosticH',separate, ...
        'sourceOnlyXBitwiseIdentical',sourceOnlyIdentical,'sourceOnlyContainsOnlyX',sourceOnlyHasOnlyX, ...
        'sharedReuse',sharedReuse, ...
        'referenceIsNonzero',any(ref.H~=0),'separateSquareUnderflows',any(c.Y~=0&(c.Y*exp(-c.ell)).^2==0), ...
        'passed',ok&&sourceOnlyIdentical&&sourceOnlyHasOnlyX&&sharedReuse.passed));
    passed=passed&&ok&&sourceOnlyIdentical&&sourceOnlyHasOnlyX&&sharedReuse.passed;
end
result=struct('records',rows,'reference','A=y_phys^2*F(x0)+O(y_phys^4), with omega and coordinate powers combined before exponentiation.', ...
    'roundoffPolicy','Absolute bound is relativeTolerance*referenceInfinity plus 128 ulps of that infinity, including subnormal values.', ...
    'passed',passed);
end

function result=shared_test(m,s,X,Y)
% The raw cache depends on model parameters, ell, shift and exact queries.
% Omega and the clock do not enter raw Gaussian values and may change.
rows=struct([]);passed=true;
for wantSecond=[true,false]
    single=ipm_cr.model.gaussian(m,s,X,Y,wantSecond);
    [dual,shared]=ipm_cr.model.gaussian(m,s,X,Y,wantSecond);
    plain=ipm_cr.model.tail(m,s,X,Y,wantSecond);
    explicitEmpty=ipm_cr.model.tail(m,s,X,Y,wantSecond,false,[]);
    reused=ipm_cr.model.tail(m,s,X,Y,wantSecond,false,shared);
    singleDouble=isequaln(single,dual);emptySame=isequaln(plain,explicitEmpty);
    reuseSame=isequaln(plain,reused);fields=fieldnames(plain);fieldEquality=struct();
    for k=1:numel(fields),name=fields{k};fieldEquality.(name)=isequaln(plain.(name),reused.(name));end
    identity=struct([]);allRejected=true;allowedScaleSame=true;
    if ~isempty(shared)
        names={'amplitude','centerX','widthX','widthY','skewWeight','ell','shift','X','Y'};
        for k=1:numel(names)
            mm=m;ss=s;xx=X;yy=Y;
            switch names{k}
                case 'amplitude',mm.amplitude=m.amplitude*1.01;
                case 'centerX',mm.centerX=m.centerX+.125;
                case 'widthX',mm.widthX=m.widthX*1.01;
                case 'widthY',mm.widthY=m.widthY*1.01;
                case 'skewWeight'
                    if m.skewWeight<=.5,mm.skewWeight=m.skewWeight+.1;
                    else,mm.skewWeight=m.skewWeight-.1;end
                case 'ell',ss.logC_l=s.logC_l+.125;
                case 'shift',ss.X_shift=s.X_shift+max(.125,1024*eps(max(1,abs(s.X_shift))));
                case 'X',xx(1)=X(1)+max(.125,1024*eps(max(1,abs(X(1)))));
                case 'Y',yy(1)=Y(1)+max(.125,1024*eps(max(1,abs(Y(1)))));
            end
            rejected=false;identifier='no_error';
            try
                ipm_cr.model.tail(mm,ss,xx,yy,wantSecond,false,shared);
            catch problem
                identifier=problem.identifier;rejected=strcmp(identifier,'ipm_cr:TailSharedIdentity');
            end
            identity=append(identity,struct('changedInput',names{k},'rejected',rejected,'identifier',identifier));
            allRejected=allRejected&&rejected;
        end
        changed=s;changed.logC_omega=s.logC_omega+.125;changed.physicalTime=.7;changed.physicalTimeLow=.03;
        allowed=ipm_cr.model.tail(m,changed,X,Y,wantSecond,false,shared);
        allowedPlain=ipm_cr.model.tail(m,changed,X,Y,wantSecond);
        allowedScaleSame=isequaln(allowed,allowedPlain);
    end
    ok=singleDouble&&emptySame&&reuseSame&&allRejected&&allowedScaleSame;
    rows=append(rows,struct('wantSecondDerivatives',wantSecond,'sharedAvailable',~isempty(shared), ...
        'gaussianSingleAndDoubleOutputBitwiseIdentical',singleDouble, ...
        'helperExplicitEmptyAndDefaultBitwiseIdentical',emptySame, ...
        'helperSharedAndUnsharedBitwiseIdentical',reuseSame,'helperFieldBitwiseEquality',fieldEquality, ...
        'identityGateExercised',~isempty(shared),'mismatchedIdentityChecks',identity, ...
        'omegaAndClockChangesReuseBitwiseIdentical',allowedScaleSame,'passed',ok));
    passed=passed&&ok;
end
result=struct('records',rows,'definition','Ephemeral raw Gaussian reuse only; empty shared selects the unchanged fallback.', ...
    'requiredIdentity','kind/model five parameters/ell/X_shift/exact X and Y; omega and physical clock are independent of raw fields.', ...
    'passed',passed);
end

function result=assembly_test(m)
s=scale0();s.logC_l=.2;s.logC_omega=.3;s.X_shift=.25;s.physicalTime=.4;s.physicalTimeLow=.03;
X=[-2,-.2,0,.7,2,1e10];Y=[0,.4,0,1,2,.8];
base=m;base.tailCorrection=false;g=ipm_cr.model.background(base,s,X,Y);
d=ipm_cr.model.tail(m,s,X,Y);b=ipm_cr.model.background(m,s,X,Y);
errors=struct();passed=true;
for key={'H','x','y','xx','xy'}
    name=key{1};expected=g.(name)+s.physicalTime*d.(name)+s.physicalTimeLow*d.(name);
    errors.(name)=metric(b.(name),expected,2e-13);passed=passed&&errors.(name).passed;
end
source=ipm_cr.model.source(m,s,X,Y);sourceError=metric(source,b.x,2e-13);
timeFields=isequaln(b.timeH,d.timeH)&&isequaln(b.timeX,d.timeX);
s.physicalTime=0;s.physicalTimeLow=0;zero=ipm_cr.model.background(m,s,X,Y);
original=ipm_cr.model.background(base,s,X,Y);unchanged=true;
for key={'H','x','y','xx','xy','originValue','originX'}
    name=key{1};unchanged=unchanged&&isequaln(zero.(name),original.(name));
end
sourceUnchanged=isequaln(ipm_cr.model.source(m,s,X,Y),ipm_cr.model.source(base,s,X,Y));
result=struct('errors',errors,'sourceError',sourceError,'timeFieldsAlreadyIncludeClockRate',timeFields, ...
    'zeroPhysicalTimeFieldsBitwiseIdentical',unchanged,'zeroPhysicalTimeSourceBitwiseIdentical',sourceUnchanged, ...
    'originHExactlyZero',b.H(3)==0&&b.timeH(3)==0, ...
    'passed',passed&&sourceError.passed&&timeFields&&unchanged&&sourceUnchanged&&b.H(3)==0&&b.timeH(3)==0);
end

function result=kinematic_test(m)
s=scale0();s.logC_l=.2;s.logC_omega=.3;s.X_shift=.2;s.physicalTime=.4;
cl=.17;cw=-.13;cr=.23;X=[-2,-.2,.7,2,30];Y=[0,.4,1,2,.8];
b=ipm_cr.model.background(m,s,X,Y);
expectedH=cw*b.H-(cl*X+cr).*b.x-cl*Y.*b.y+cr*b.originX+b.timeH;
expectedX=(cw-cl)*b.x-(cl*X+cr).*b.xx-cl*Y.*b.xy+b.timeX;
steps=[.04,.02,.01];rows=struct([]);
for h=steps
    offsets=[-2,-1,1,2];H=cell(1,4);Q=H;
    for k=1:4
        tau=offsets(k)*h;t=s;t.logC_l=s.logC_l+cl*tau;t.logC_omega=s.logC_omega+cw*tau;
        t.X_shift=exp(cl*tau)*s.X_shift+cr*expm1(cl*tau)/cl;
        t.physicalTime=s.physicalTime+exp(s.logC_omega-s.logC_l)*expm1((cw-cl)*tau)/(cw-cl);
        z=ipm_cr.model.background(m,t,X,Y);H{k}=z.H;Q{k}=z.x;
    end
    hdiff=(H{1}-8*H{2}+8*H{3}-H{4})/(12*h);
    qdiff=(Q{1}-8*Q{2}+8*Q{3}-Q{4})/(12*h);
    rows=append(rows,struct('differenceStep',h,'HError',metric(hdiff,expectedH,2e-8), ...
        'XError',metric(qdiff,expectedX,2e-8)));
end
result=struct('records',rows,'identity','H_tau=cw*H-(cl*X+cr)*Hx-cl*Y*Hy+cr*Hx(0)+timeH; (Hx)_tau=(cw-cl)*Hx-(cl*X+cr)*Hxx-cl*Y*Hxy+timeX', ...
    'physicalPoissonOverwritten',false,'passed',rows(end).HError.passed&&rows(end).XError.passed);
end

function result=initial_test(overrides)
[state,g,wall,method]=ipm_cr.evolve.init(overrides);
assert(isfield(state,'residual')&&method.background.tailCorrection,'ipm_cr:TailTestInitial');
assert(state.scale.physicalTime==0&&state.scale.physicalTimeLow==0,'ipm_cr:TailTestInitialTime');
base=method;base.background.tailCorrection=false;ctx=ipm_cr.context(g,wall,base);
[old,rho]=ipm_cr.evolve.residual(state.residual,state.scale,g,wall,base, ...
    [state.flow.c_l,state.flow.c_r],true,ctx);
rows=struct([]);passed=true;
for patch={'core','polar'}
    key=patch{1};grid=g.(key);active=grid.active;
    b=ipm_cr.model.background(method.background,state.scale,grid.X,grid.Y);
    errors=struct('rho',metric(state.rho.(key),rho.(key),2e-12), ...
        'q',metric(state.flow.q.(key),old.q.(key),2e-12), ...
        'activeCompleteF',metric(state.flow.F.(key)(active),old.F.(key)(active),2e-11), ...
        'activeResidualTimeTerm',metric(state.flow.residualRhs.(key)(active)+b.timeH(active),old.residualRhs.(key)(active),2e-11), ...
        'allArrayCompleteFDiagnostic',metric(state.flow.F.(key),old.F.(key),2e-11));
    rhoIdentical=isequaln(state.rho.(key),rho.(key));qIdentical=isequaln(state.flow.q.(key),old.q.(key));
    ok=errors.rho.passed&&errors.q.passed&&rhoIdentical&&qIdentical&& ...
        errors.activeCompleteF.passed&&errors.activeResidualTimeTerm.passed;
    rows=append(rows,struct('patch',key,'errors',errors,'activeSampleCount',nnz(active), ...
        'rhoBitwiseIdentical',rhoIdentical,'qBitwiseIdentical',qIdentical, ...
        'allActivePointsRetained',true,'passed',ok));passed=passed&&ok;
end
rates=[state.flow.c_l,state.flow.c_omega,state.flow.c_r,state.flow.Uref,state.flow.transportOrigin];
oldRates=[old.c_l,old.c_omega,old.c_r,old.Uref,old.transportOrigin];rateError=metric(rates,oldRates,2e-11);
timeTermsNonzero=any(state.flow.residualRhs.core~=old.residualRhs.core,'all')|| ...
    any(state.flow.residualRhs.polar~=old.residualRhs.polar,'all');
result=struct('tau',state.scale.canonicalTime,'physicalTime',state.scale.physicalTime, ...
    'records',rows,'ratesError',rateError,'etaExactlyZero',all(state.residual.core==0,'all')&&all(state.residual.polar==0,'all'), ...
    'timeTermsNonzero',timeTermsNonzero, ...
    'ghostScope','All-array differences are retained; ghosts are algebraic interpolants of eta_rhs, so only owning active DOFs have direct timeH cancellation.', ...
    'passed',passed&&rateError.passed&&timeTermsNonzero&& ...
        all(state.residual.core==0,'all')&&all(state.residual.polar==0,'all'));
end

function value=vertical_A(m,x,y)
rho=primitive(m,x,y);C=m.amplitude*m.widthX*m.widthY^2/(2*sqrt(pi));
value=(-4*C/m.widthY^2)*x.*y.^2.*rho./(1+x.^2+y.^2).^2;
end

function value=vertical_coefficient(m,x)
rho=primitive(m,x,zeros(size(x)));C=m.amplitude*m.widthX*m.widthY^2/(2*sqrt(pi));
value=(-4*C/m.widthY^2)*x.*rho./(1+x.^2).^2;
end

function [rho,qx,qy]=primitive(m,x,y)
z=(x-m.centerX)/m.widthX;v=exp(-(y/m.widthY).^2);w=m.skewWeight;
rho=m.amplitude*m.widthX*sqrt(pi)/2*((1-w)*erfc(-z)+w*erfc(-(z-.8)/1.25)).*v;
if nargout>1,qx=m.amplitude*((1-w)*exp(-z.^2)+(w/1.25)*exp(-((z-.8)/1.25).^2)).*v;end
if nargout>2,qy=-2*y/m.widthY^2.*rho;end
end

function d=fivepoint(f,x,h)
d=(f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))/(12*h);
end

function d=second_fivepoint(f,x,h)
d=(-f(x-2*h)+16*f(x-h)-30*f(x)+16*f(x+h)-f(x+2*h))/(12*h^2);
end

function d=mixed_fivepoint(f,x,y,h)
offsets=[-2,-1,1,2];weights=[1,-8,8,-1];d=zeros(size(x));
for j=1:4
    for k=1:4,d=d+weights(j)*weights(k)*f(x+offsets(j)*h,y+offsets(k)*h);end
end
d=d/(144*h^2);
end

function s=scale0()
s=struct('logC_l',0,'logC_omega',0,'X_shift',0,'physicalTime',0,'physicalTimeLow',0);
end

function [x,w]=gauss(n,a,b)
k=(1:n-1)';beta=k./sqrt(4*k.^2-1);[V,D]=eig(diag(beta,1)+diag(beta,-1));
[nodes,order]=sort(diag(D));x=(a+b)/2+(b-a)/2*nodes;
w=(b-a)*(V(1,order).^2)';
end

function info=metric(actual,expected,relativeTolerance)
assert(isequal(size(actual),size(expected)),'ipm_cr:TailTestMetricShape');
error=max(abs(actual(:)-expected(:)));scale=max(abs(expected(:)));
tolerance=relativeTolerance*scale+128*eps(scale);
relative=0;if scale~=0,relative=error/scale;elseif error~=0,relative=Inf;end
info=struct('absoluteInfinity',error,'referenceInfinity',scale,'relativeInfinity',relative, ...
    'absoluteTolerance',tolerance,'allFinite',all(isfinite(actual),'all')&&all(isfinite(expected),'all'), ...
    'passed',all(isfinite(actual),'all')&&all(isfinite(expected),'all')&&error<=tolerance);
end

function inventory=source_inventory(root)
files=dir(fullfile(root,'src','+ipm_cr','**','*.m'));paths={mfilename('fullpath')};
if ~endsWith(paths{1},'.m'),paths{1}=[paths{1},'.m'];end
for j=1:numel(files),paths{end+1}=fullfile(files(j).folder,files(j).name);end %#ok<AGROW>
for name={'ipm_cr.op.wenoMex'}
    native=which(name{1});if ~isempty(native),paths{end+1}=native;end %#ok<AGROW>
end
paths{end+1}=fullfile(root,'src','+ipm_cr','+op','wenoMex.cpp');
paths=sort(unique(paths));inventory=struct([]);
for j=1:numel(paths)
    f=dir(paths{j});inventory=append(inventory,struct('path',paths{j},'bytes',f.bytes,'sha256',ipm_cr.io.sha(paths{j})));
end
end

function persist(folder,report)
temporary=fullfile(folder,'tail_background.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'tail_background.mat'),'f');assert(ok,'ipm_cr:TailTestSave','%s',message);
temporary=fullfile(folder,'tail_background.json.writing');fid=fopen(temporary,'w');assert(fid>=0,'ipm_cr:TailTestSave');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'tail_background.json'),'f');assert(ok,'ipm_cr:TailTestSave','%s',message);
end

function a=append(a,b)
if isempty(a),a=b;else,a(end+1)=b;end
end
