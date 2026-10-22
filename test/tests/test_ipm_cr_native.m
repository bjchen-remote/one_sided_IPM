function report=test_ipm_cr_native(outputDir,options)
%TEST_IPM_CR_NATIVE Load and check the optional native reconstruction kernel.
% Polynomial exactness, constants and frozen-density velocity derivatives
% are local mathematical checks, not a global IPM accuracy certificate.
if nargin<2,options=struct();end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'time_evolution'));
if nargin<1||isempty(outputDir)
    [~,token]=fileparts(tempname(fullfile(root,'test','results')));
    outputDir=fullfile(root,'test','results',['native_test_',token]);
end
outputDir=char(outputDir);if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'native_test.json')),'ipm_cr:NativeTestPreserve');
timer=tic;report=struct('kind','native_reconstruction_local_mathematical_checks_v1', ...
    'status','running','passed',false,'globalAccuracyCertified',false, ...
    'poissonCalls',0,'timeAdvanced',0,'matlabVersion',version);
try
    binary=which('ipm_cr.explicit.wenoPairSmoothMex');
    assert(~isempty(binary)&&exist(binary,'file')==3,'ipm_cr:NativeTestBinary', ...
        'Build the native kernel before running this test.');
    source=fullfile(root,'test','+ipm_cr','+explicit','wenoPairSmoothMex.cpp');
    report.native=struct('binary',binary,'binarySha256',ipm_c.output.sha256(binary), ...
        'source',source,'sourceSha256',ipm_c.output.sha256(source));
    report.build=build_evidence(root,report.native,options);
    names={'ipm_cr.explicit.transport','ipm_cr.explicit.geometry', ...
        'ipm_cr.explicit.synchronize','build_ipm_cr_native','test_ipm_cr_native'};
    report.sources=struct([]);
    for k=1:numel(names)
        path=which(names{k});entry=struct('path',path,'sha256',ipm_c.output.sha256(path));
        if isempty(report.sources),report.sources=entry;else,report.sources(end+1)=entry;end %#ok<AGROW>
    end

    % A degree-five continuation gives an independent exact derivative at
    % every node, including both boundary nodes and their exterior ghosts.
    x=-1:1/32:1;dx=x(2)-x(1);q=zeros(6,numel(x));exact=q;
    for degree=0:5
        q(degree+1,:)=x.^degree;
        if degree>0,exact(degree+1,:)=degree*x.^(degree-1);end
    end
    mask=true(6,1);
    [plus,minus]=ipm_cr.explicit.wenoPairSmoothMex(q,dx,1e-12,5,mask,mask,0);
    errors=max(max(abs(plus-exact),[],2),max(abs(minus-exact),[],2));
    report.polynomial=struct('degrees',0:5,'maximumErrors',errors', ...
        'boundaryError',max(abs([plus(:,[1,end])-exact(:,[1,end]), ...
            minus(:,[1,end])-exact(:,[1,end])]),[],'all'), ...
        'tolerance',5e-10,'allNodesIncluded',true);
    assert(max(errors)<=report.polynomial.tolerance,'ipm_cr:NativeTestPolynomial', ...
        'The linear kernel failed polynomial differentiation, including boundary ghosts.');

    constant=2^40*ones(2,numel(x));mask=false(2,1);constantErrors=zeros(1,2);
    for mode=0:1
        [a,b]=ipm_cr.explicit.wenoPairSmoothMex(constant,dx,1e-12,5,mask,mask,mode);
        constantErrors(mode+1)=max(abs([a,b]),[],'all');
    end
    assert(all(constantErrors==0),'ipm_cr:NativeTestConstant','A constant must differentiate to exactly zero.');
    % Dyadic quadratic samples survive this large constant addition exactly.
    % The candidate reconstructions agree for this polynomial; this avoids
    % falsely claiming arbitrary-field offset invariance of RMS weights.
    quadratic=.5*x+.125*x.^2;derivative=.5+.25*x;offset=2^30;
    [a,b]=ipm_cr.explicit.wenoPairSmoothMex(quadratic,dx,1e-12,5,false,false,1);
    [shiftA,shiftB]=ipm_cr.explicit.wenoPairSmoothMex(quadratic+offset,dx,1e-12,5,false,false,1);
    [defaultA,defaultB]=ipm_cr.explicit.wenoPairSmoothMex(quadratic,dx,1e-12,5,false,false);
    translationError=max(abs([shiftA-a,shiftB-b]),[],'all');
    derivativeError=max(abs([shiftA-derivative,shiftB-derivative]),[],'all');
    report.constants=struct('value',2^40,'exactZeroErrorsByMode',constantErrors, ...
        'smoothQuadraticOffset',offset,'translationDerivativeError',translationError, ...
        'shiftedExactDerivativeError',derivativeError,'sixArgumentDefaultIdentical', ...
        isequal(a,defaultA)&&isequal(b,defaultB),'arbitraryFieldOffsetInvarianceClaimed',false);
    assert(all(isfinite([a,b,shiftA,shiftB]),'all')&&translationError<1e-10&& ...
        derivativeError<1e-10&&report.constants.sixArgumentDefaultIdentical, ...
        'ipm_cr:NativeTestSmooth','Smooth reconstruction failed the finite/constant-shift/default-mode check.');

    g=ipm_cr.explicit.geometry();c=g.core;p=g.polar;
    % Localized Gaussian residual: exactly zero at the origin and decays
    % to zero in the far field, instead of a centered-density plateau.
    gaussian=@(X,Y)(X.^2+Y.^2).*exp(-((X-.45).^2+Y.^2)/.4^2);
    rho=struct('core',gaussian(c.X,c.Y),'polar',gaussian(p.X,p.Y));
    rho=ipm_cr.explicit.synchronize(rho,g);
    u=struct('core',.24+.04*exp(-(c.X.^2+c.Y.^2)), ...
        'polar',.24+.04*exp(-(p.X.^2+p.Y.^2)));
    v=struct('core',.03*c.Y.*exp(-(c.X.^2+c.Y.^2)), ...
        'polar',.03*p.Y.*exp(-(p.X.^2+p.Y.^2)));
    sigma=.1;nu=sigma/min(c.dx,c.dy);cl=-.07;delta=.11;h=2e-5;
    report.transport=struct([]);
    for reconstruction={'linear5','smooth_weno5'}
        tr=struct('reconstruction',reconstruction{1},'smoothSpeedCoefficient',nu);
        prepared=ipm_cr.explicit.transport(rho,g,tr);
        [vx,vy]=velocity(u,v,g,cl,delta);
        [B,rate,position,partials]=prepared(vx,vy);
        checks=struct();
        for name={'c_l','delta'}
            dc=0;dd=0;if strcmp(name{1},'c_l'),dc=h;else,dd=h;end
            [ax,ay]=velocity(u,v,g,cl+dc,delta+dd);upper=prepared(ax,ay);
            [ax,ay]=velocity(u,v,g,cl-dc,delta-dd);lower=prepared(ax,ay);
            finiteDifference=struct('core',(upper.core-lower.core)/(2*h), ...
                'polar',(upper.polar-lower.polar)/(2*h));
            checks.(name{1})=difference(finiteDifference,partials.(name{1}),g);
            assert(checks.(name{1}).maximumMixedError<2e-7, ...
                'ipm_cr:NativeTestPartials','%s %s partial failed central differences.',reconstruction{1},name{1});
        end
        poisoned=vy;poisoned.core(1,:)=1e6;poisoned.polar(:,[1,end])=-1e6;
        [wallB,wallRate]=prepared(vx,poisoned);
        wallExact=isequal(B,wallB)&&isequal(rate,wallRate);
        assert(wallExact,'ipm_cr:NativeTestWall','A physical-wall normal speed affected transport/CFL.');
        z=struct('core',zeros(size(c.X)),'polar',zeros(size(p.X)));
        [~,floorRate,floorPosition]=prepared(z,z);
        % At an active Cartesian interior node both physical directions
        % retain sigma; an active annular interior node also retains both.
        coreBound=sigma*(1/c.dx+1/c.dy);
        polarBound=sigma/p.drMetric(p.activeRows(1))+ ...
            sigma/(p.r(p.activeRows(1))*p.dtheta);
        assert(floorRate>=max(coreBound,polarBound)*(1-1e-12)&& ...
            floorPosition.maximumCoreRate>0&&floorPosition.maximumPolarRate>0&& ...
            position.allActiveOverlapIncluded&&~position.outerDegreeSwitching, ...
            'ipm_cr:NativeTestFloor','The full active mapped CFL omitted a speed floor or an overlap patch.');
        % Production transport removes a common origin before weights. The
        % floating-point error from adding a large constant remains visible.
        shifted=struct('core',rho.core+2^20,'polar',rho.polar+2^20);
        shiftedB=ipm_cr.explicit.transport(shifted,vx,vy,g,tr);
        shift=difference(shiftedB,B,g);
        assert(shift.maximumMixedError<2e-6,'ipm_cr:NativeTestTransportShift', ...
            'Production centered transport is unstable under a representable large constant shift.');
        row=struct('reconstruction',reconstruction{1},'physicalSpeedFloor',sigma, ...
            'smoothSpeedCoefficient',nu,'rate',rate,'position',position, ...
            'partialStep',h,'partials',checks,'wallNormalExactlyIgnored',wallExact, ...
            'zeroVelocityFloorRate',floorRate,'floorPosition',floorPosition, ...
            'independentFloorLowerBounds',[coreBound,polarBound], ...
            'largeConstantShiftError',shift,'finiteBothFullPatches', ...
            all(isfinite(B.core),'all')&&all(isfinite(B.polar),'all'));
        if isempty(report.transport),report.transport=row;else,report.transport(end+1)=row;end %#ok<AGROW>
    end
    report.sourceAndBinaryUnchanged=strcmp(report.native.binarySha256,ipm_c.output.sha256(binary))&& ...
        strcmp(report.native.sourceSha256,ipm_c.output.sha256(source));
    assert(report.sourceAndBinaryUnchanged,'ipm_cr:NativeTestSourceChanged');
    report.status='completed';report.passed=true;report.elapsedSeconds=toc(timer);persist(outputDir,report);
catch problem
    report.status='failed';report.elapsedSeconds=toc(timer);
    report.error=struct('identifier',problem.identifier,'message',problem.message,'stack',problem.stack);
    persist(outputDir,report);rethrow(problem);
end
fprintf('NATIVE local checks passed=%d polynomial %.3e smooth shift %.3e %.3fs\n', ...
    report.passed,max(errors),translationError,report.elapsedSeconds);
end

function [vx,vy]=velocity(u,v,g,cl,delta)
vx=struct('core',u.core+cl*g.core.X+delta,'polar',u.polar+cl*g.polar.X+delta);
vy=struct('core',v.core+cl*g.core.Y,'polar',v.polar+cl*g.polar.Y);
end

function result=difference(a,b,g)
result=struct();result.maximumMixedError=0;
for patch={'core','polar'}
    name=patch{1};mask=g.(name).active;d=abs(a.(name)-b.(name));
    scale=max(1,max(abs(b.(name)),[],'all'));
    result.(name)=struct('fullAbsoluteInfinity',max(d,[],'all'), ...
        'activeAbsoluteInfinity',max(d(mask)),'referenceInfinity',max(abs(b.(name)),[],'all'), ...
        'mixedError',max(d,[],'all')/scale);
    result.maximumMixedError=max(result.maximumMixedError,result.(name).mixedError);
end
end

function evidence=build_evidence(root,native,options)
if isfield(options,'buildReport'),paths={char(options.buildReport)};
else
    listing=dir(fullfile(root,'test','results','native_build*','native_build.json'));
    paths=arrayfun(@(a)fullfile(a.folder,a.name),listing,'UniformOutput',false);
end
evidence=struct();
for k=numel(paths):-1:1
    build=jsondecode(fileread(paths{k}));
    if isfield(build,'binarySha256')&&strcmp(build.binarySha256,native.binarySha256)&& ...
            strcmp(build.sourceSha256,native.sourceSha256)&&strcmp(build.status,'completed')
        evidence=struct('path',paths{k},'sha256',ipm_c.output.sha256(paths{k}), ...
            'binarySha256',build.binarySha256,'sourceSha256',build.sourceSha256);return
    end
end
error('ipm_cr:NativeTestBuildEvidence','No completed source/binary SHA-matching build report found.');
end

function persist(folder,report)
save(fullfile(folder,'native_test.mat'),'report');
fid=fopen(fullfile(folder,'native_test.json'),'w');assert(fid>=0,'ipm_cr:NativeTestReport');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');
end
