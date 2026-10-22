function report=test_ipm_cr_background(outputDir,options)
%TEST_IPM_CR_BACKGROUND Independent centered-background arithmetic checks.
% H is compared with line integrals of independently coded Gx and Gy, not
% with subtraction of two primitive values or with the implementation's
% erf-increment recurrence. No physical Poisson field is forced to zero.
if nargin<2,options=struct();end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'time_evolution'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'test','results','background');end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),{'includeInitial','initialOptions'})), ...
    'ipm_cr:BackgroundTestOptions');
includeInitial=true;if isfield(options,'includeInitial'),includeInitial=options.includeInitial;end
initialOptions=struct();if isfield(options,'initialOptions'),initialOptions=options.initialOptions;end
validateattributes(includeInitial,{'logical'},{'scalar'});
outputDir=char(outputDir);assert(startsWith(outputDir,filesep),'ipm_cr:BackgroundTestOutput');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'background.json')),'ipm_cr:BackgroundTestPreserve','Choose a fresh output directory.');
timer=tic;cfg=ipm_c.config();p=cfg.physics;
model=struct('kind','gaussian_primitive','amplitude',p.gaussianAmplitude, ...
    'centerX',p.gaussianCenterX,'widthX',p.gaussianWidthX, ...
    'widthY',p.gaussianWidthY,'skewWeight',p.gaussianSkewWeight);
report=struct('kind','analytic_centered_background_arithmetic_and_kinematic_checks_v1', ...
    'status','running','timeAdvanced',0,'model',model,'arithmetic',struct([]), ...
    'kinematic',struct(),'initial',struct(),'sourceInventory',source_inventory(root), ...
    'sourceHashesUnchanged',false,'arithmeticPassed',false,'kinematicPassed',false, ...
    'initialPassed',false,'allRequestedChecksPassed',false,'elapsedSeconds',0, ...
    'HReference','Integral_0^X Gx(s,0) ds + Integral_0^Y Gy(X,s) ds using independent elementary analytic derivatives.', ...
    'kinematicScope','Prescribed pure gauge coordinate/amplitude motion with physical U=0; no nonzero-source IPM Darcy field is overwritten.', ...
    'limitations',{{'A kinematic eta=0 identity is not a solution of a nonzero-source IPM problem with its Darcy velocity suppressed.', ...
        'Exact initial residual and analytic-source checks do not certify future residual spatial resolution or the whole trajectory.'}});
persist(outputDir,report);
try
    cases=arithmetic_cases();rows=struct([]);
    for j=1:numel(cases)
        c=cases(j);b=ipm_cr.explicit.background(model,c.scale,c.X,c.Y);
        light=ipm_cr.explicit.background(model,c.scale,c.X,c.Y,false);
        sourceOnly=ipm_cr.explicit.backgroundSource(model,c.scale,c.X,c.Y);
        ref=reference_fields(model,c.scale,c.X,c.Y);href=zeros(size(c.X));
        for k=1:numel(c.X)
            x=c.X(k);y=c.Y(k);
            href(k)=integral(@(u)reference_component(model,c.scale,x*u,zeros(size(u)),'x')*x, ...
                0,1,'AbsTol',0,'RelTol',2e-13)+ ...
                integral(@(u)reference_component(model,c.scale,x+zeros(size(u)),y*u,'y')*y, ...
                0,1,'AbsTol',0,'RelTol',2e-13);
        end
        metrics=struct('H',component_metric(b.H,href), ...
            'sourceOnlyX',component_metric(sourceOnly,ref.x));
        for name={'x','y','xx','xy'}
            key=name{1};metrics.(key)=component_metric(b.(key),ref.(key));
        end
        origin=ipm_cr.explicit.background(model,c.scale,0,0);
        direct=ref.G-ref.originValue;
        lightIdentical=isequaln(light.H,b.H)&&isequaln(light.x,b.x)&&isequaln(light.y,b.y)&& ...
            isequaln(light.originValue,b.originValue)&&isequaln(light.originX,b.originX);
        passed=true;keys=fieldnames(metrics);
        for k=1:numel(keys)
            metric=metrics.(keys{k});passed=passed&&metric.allFinite&&metric.maximumNonzeroReferenceRelativeError<=1e-10&& ...
                metric.maximumZeroReferenceAbsoluteError<=1e-12;
        end
        assert(origin.H==0,'ipm_cr:BackgroundTestOrigin','The centered background origin must be exactly zero.');
        passed=passed&&lightIdentical&&isempty(light.xx)&&isempty(light.xy);
        rows=append(rows,struct('name',c.name,'scale',c.scale,'X',c.X,'Y',c.Y, ...
            'errors',metrics,'computedH',b.H,'integralReferenceH',href, ...
            'directPrimitiveSubtraction',direct, ...
            'directSubtractionError',component_metric(direct,href), ...
            'originHExactlyZero',true,'originValueError',component_metric(origin.originValue,ref.originValue), ...
            'originXError',component_metric(origin.originX,reference_component(model,c.scale,0,0,'x')), ...
            'lightHAndFirstDerivativesBitwiseIdentical',lightIdentical, ...
            'lightSecondDerivativesEmpty',isempty(light.xx)&&isempty(light.xy), ...
            'passed',passed));
    end
    report.arithmetic=rows;report.arithmeticPassed=all([rows.passed]);
    noneModel=struct('kind','none');noneScale=cases(1).scale;
    none=ipm_cr.explicit.background(noneModel,noneScale,[0,1],[0,.5],false);
    noneSource=ipm_cr.explicit.backgroundSource(noneModel,noneScale,[0,1],[0,.5]);
    report.noneModel=struct('sourceExactlyZero',all(noneSource==0,'all'), ...
        'HExactlyZero',all(none.H==0,'all'),'secondDerivativesEmpty',isempty(none.xx)&&isempty(none.xy));
    report.arithmeticPassed=report.arithmeticPassed&&report.noneModel.sourceExactlyZero&& ...
        report.noneModel.HExactlyZero&&report.noneModel.secondDerivativesEmpty;
    report.kinematic=kinematic_test(model);report.kinematicPassed=report.kinematic.passed;
    report.elapsedSeconds=toc(timer);persist(outputDir,report);
    if includeInitial
        [state,g,~,method]=ipm_cr.explicit.initialize(initialOptions);
        assert(isfield(state,'residual'),'ipm_cr:BackgroundTestInitial','The default new path must store eta separately.');
        etaMax=max([abs(state.residual.core(:));abs(state.residual.polar(:))]);
        qError=struct();bareError=struct();WError=struct();
        raw=ipm_cr.explicit.differentiate(state.rho,g);
        for patch={'core','polar'}
            key=patch{1};grid=g.(key);
            exact=reference_fields(method.background,state.scale,grid.X,grid.Y);
            bg=ipm_cr.explicit.background(method.background,state.scale,grid.X,grid.Y);
            qError.(key)=native_metric(state.flow.q.(key)-exact.x,grid);
            bareError.(key)=native_metric(raw.(key).x-exact.x,grid);
            WError.(key)=native_metric(state.rho.(key)-bg.H,grid);
        end
        qMax=max(qError.core.globalInfinity,qError.polar.globalInfinity);
        report.initial=struct('tau',state.scale.canonicalTime,'etaInfinity',etaMax, ...
            'etaExactlyZero',etaMax==0,'originWExactlyZero',state.rho.core(1,g.core.originIndex)==0, ...
            'analyticSourceEnabled',state.flow.analyticBackgroundSource, ...
            'qError',qError,'bareFullWDifferentiationDiagnostic',bareError, ...
            'WReconstructionError',WError,'poisson',state.flow.poisson, ...
            'passed',etaMax==0&&qMax<=1e-12&&state.flow.analyticBackgroundSource);
        report.initialPassed=report.initial.passed;
    end
    report.sourceHashesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
    assert(report.sourceHashesUnchanged,'ipm_cr:BackgroundTestSources');
    report.allRequestedChecksPassed=report.arithmeticPassed&&report.kinematicPassed&&(~includeInitial||report.initialPassed);
    if report.allRequestedChecksPassed,report.status='passed_requested_background_checks';
    else,report.status='completed_measurements_with_failed_checks';end
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function cases=arithmetic_cases()
cases=struct('name',{},'scale',{},'X',{},'Y',{});
cases(1)=struct('name','ordinary_and_tiny_increments', ...
    'scale',struct('logC_l',0,'logC_omega',0,'X_shift',.2), ...
    'X',[-1.1,-.1,-1e-12,0,1e-12,.3,1.2],'Y',[0,.5,.1,0,.3,.4,1]);
cases(2)=struct('name','X_lost_when_subtracting_large_shift', ...
    'scale',struct('logC_l',40,'logC_omega',80,'X_shift',.35*exp(40)), ...
    'X',[0,1,-1,0,.3],'Y',[0,0,0,1,.7]);
cases(3)=struct('name','joint_exponents_survive_separate_vertical_square_underflow', ...
    'scale',struct('logC_l',400,'logC_omega',600,'X_shift',.35*exp(400)), ...
    'X',[0,1,-1,0,.3],'Y',[0,0,0,1,.7]);
end

function result=kinematic_test(model)
s=struct('logC_l',.3,'logC_omega',.7,'X_shift',.2);
cl=.23;cw=-.17;cr=.31;
X=[-2,-1,.2,1,2,1e10];Y=[0,.4,.8,1,.1,1e10];
b=ipm_cr.explicit.background(model,s,X,Y);
expected=cw*b.H-(cl*X+cr).*b.x-cl*Y.*b.y+cr*b.originX;
steps=[.04,.02,.01];errors=zeros(size(steps));records=struct([]);
for j=1:numel(steps)
    h=steps(j);values=cell(1,4);offsets=[-2,-1,1,2];
    for k=1:4
        t=offsets(k)*h;scale=s;scale.logC_l=s.logC_l+cl*t;
        scale.logC_omega=s.logC_omega+cw*t;
        scale.X_shift=exp(cl*t)*s.X_shift+cr*expm1(cl*t)/cl;
        v=ipm_cr.explicit.background(model,scale,X,Y);values{k}=v.H;
    end
    derivative=(values{1}-8*values{2}+8*values{3}-values{4})/(12*h);
    m=component_metric(derivative,expected);errors(j)=max(abs(derivative-expected));
    records=append(records,struct('differenceStep',h,'errors',m, ...
        'etaRateResidualFromIndependentTimeDifference',derivative-expected));
end
scale=max(1,max(abs(expected)));
result=struct('pureGaugeRates',[cl,cw,cr],'scale',s,'X',X,'Y',Y, ...
    'identity','H_tau=cw*H-(cl*X+cr)*Gx-cl*Y*Gy+cr*Gx(0)', ...
    'etaZeroIsInvariantForPureGaugeKinematics',true, ...
    'physicalPoissonModified',false,'records',records,'absoluteInfinityErrors',errors, ...
    'observedDifferenceOrders',log2(errors(1:end-1)./errors(2:end)), ...
    'passed',errors(end)/scale<=1e-8);
end

function value=reference_component(model,s,X,Y,name)
b=reference_fields(model,s,X,Y);value=b.(name);
end

function b=reference_fields(model,s,X,Y)
% Independent elementary formula; H itself is deliberately absent.
ell=s.logC_l;omega=s.logC_omega;A=model.amplitude;sx=model.widthX;sy=model.widthY;
z=((X-s.X_shift)*exp(-ell)-model.centerX)/sx;
z0=(-s.X_shift*exp(-ell)-model.centerX)/sx;
eta2=(Y*exp(-ell)/sy).^2;weights=[1-model.skewWeight,model.skewWeight];
shifts=[0,.8];stretches=[1,1.25];P=0;P0=0;q=0;qq=0;
for k=1:2
    t=(z-shifts(k))/stretches(k);t0=(z0-shifts(k))/stretches(k);
    e=exp(-t.^2-eta2);P=P+weights(k)*erfc(-t);P0=P0+weights(k)*erfc(-t0);
    q=q+weights(k)/stretches(k)*e;
    qq=qq-2*weights(k)*t/(sx*stretches(k)^2).*e;
end
C=A*sx*sqrt(pi)/2;
b=struct('x',A*exp(omega-ell)*q, ...
    'y',-2*C/sy^2*exp(omega-2*ell)*Y.*P.*exp(-eta2), ...
    'xx',A*exp(omega-2*ell)*qq, ...
    'xy',-2*A/sy^2*exp(omega-3*ell)*Y.*q, ...
    'G',C*exp(omega)*P.*exp(-eta2),'originValue',C*exp(omega)*P0);
end

function info=component_metric(actual,expected)
d=abs(actual-expected);nonzero=expected~=0;
relative=zeros(size(expected));relative(nonzero)=d(nonzero)./abs(expected(nonzero));
zeroError=0;if any(~nonzero,'all'),zeroError=max(d(~nonzero));end
info=struct('absoluteInfinity',max(d(:)), ...
    'maximumNonzeroReferenceRelativeError',max(relative(:)), ...
    'maximumZeroReferenceAbsoluteError',zeroError, ...
    'allFinite',all(isfinite(actual),'all')&&all(isfinite(expected),'all'));
end

function info=native_metric(value,p)
mask=p.active;indices=find(mask);[v,k]=max(abs(value(indices)));index=indices(k);
near=mask&hypot(p.X,p.Y)<=3;remote=mask&hypot(p.X,p.Y)>3;
info=struct('globalInfinity',v,'nearInfinity',max(abs(value(near))), ...
    'remoteInfinity',0,'maximumLocation',struct('X',p.X(index),'Y',p.Y(index),'r',hypot(p.X(index),p.Y(index))));
if any(remote,'all'),info.remoteInfinity=max(abs(value(remote)));end
end

function inventory=source_inventory(root)
files=dir(fullfile(root,'test','+ipm_cr','**','*.m'));paths={mfilename('fullpath')};
if ~endsWith(paths{1},'.m'),paths{1}=[paths{1},'.m'];end
for j=1:numel(files),paths{end+1}=fullfile(files(j).folder,files(j).name);end %#ok<AGROW>
native=which('ipm_c.field.weno5UnitPairMex');if ~isempty(native),paths{end+1}=native;end
paths{end+1}=fullfile(root,'test','+ipm_cr','+explicit','wenoPairSmoothMex.cpp');
native=which('ipm_cr.explicit.wenoPairSmoothMex');if ~isempty(native),paths{end+1}=native;end
paths=sort(unique(paths));inventory=struct([]);
for j=1:numel(paths)
    item=dir(paths{j});inventory=append(inventory,struct('path',paths{j},'bytes',item.bytes, ...
        'sha256',ipm_c.output.sha256(paths{j})));
end
end

function persist(folder,report)
temporary=fullfile(folder,'background.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'background.mat'),'f');assert(ok,'ipm_cr:BackgroundTestSave','%s',message);
temporary=fullfile(folder,'background.json.writing');fid=fopen(temporary,'w');assert(fid>=0,'ipm_cr:BackgroundTestSave');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'background.json'),'f');assert(ok,'ipm_cr:BackgroundTestSave','%s',message);
end

function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
