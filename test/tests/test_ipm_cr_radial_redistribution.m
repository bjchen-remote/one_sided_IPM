function report=test_ipm_cr_radial_redistribution(outputDir,options)
%TEST_IPM_CR_RADIAL_REDISTRIBUTION Static invariants and analytic accuracy.
% No IVP, Poisson or gauge is advanced. Original level_2 checkpoint geometry
% is optional; when available its SHA and disabled bitwise identity are checked.
if nargin<2,options=struct();end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'test'),fullfile(root,'time_evolution'));
o=struct('baselineStudyDir',fullfile(root,'test','results','tau10_order_20261022_1323','study'), ...
    'radialRedistribution',struct('enabled',true,'factor',2,'preserveRadius',4,'fullRefinementRadius',20));
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(o))), ...
    'ipm_cr:RedistributionTestOptions','Options are baselineStudyDir (empty disables optional checkpoint comparison) and radialRedistribution.');
names=fieldnames(options);for k=1:numel(names),o.(names{k})=options.(names{k});end
assert(isstruct(o.radialRedistribution)&&isscalar(o.radialRedistribution)&&o.radialRedistribution.enabled&& ...
    o.radialRedistribution.preserveRadius==4&&o.radialRedistribution.fullRefinementRadius==20, ...
    'ipm_cr:RedistributionTestContract','This acceptance test retains r<=4 and tests the transition to full refinement at r=20.');
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'test','results','radial_redistribution');end
outputDir=char(outputDir);assert(startsWith(outputDir,filesep),'ipm_cr:RedistributionTestOutput','Use an absolute output folder.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'radial_redistribution.json')),'ipm_cr:RedistributionTestFresh','Use a fresh result folder.');
timer=tic;report=struct('kind','static_prefix_preserving_radial_redistribution_v1','status','initializing', ...
    'options',o,'records',struct([]),'orders',struct([]),'baselineCheckpoint',struct(), ...
    'sourceInventory',source_inventory(root),'sourceHashesUnchanged',false, ...
    'timeAdvanced',0,'poissonCalls',0,'gaugeCalls',0,'passed',false, ...
    'globalIVPOrderCertified',false,'trajectoryAccuracyCertified',false, ...
    'scope','Static geometry, exchange, same-velocity CFL, gradient and production mapped-transport tests. All native active points retained.', ...
    'orderRule','At least fourth order for resolved native-field errors; roundoff-limited pairs are explicitly unresolved, not a false failure or an order certificate.');
persist(outputDir,report);
try
    [baseline,report.baselineCheckpoint]=load_optional_baseline(o.baselineStudyDir);
    for level=1:3
        h=.06/2^(level-1);nr=192*2^(level-1)+1;nt=64*2^(level-1)+1;
        geom=struct('coreSpacing',h,'grid',struct('nr',nr,'ntheta',nt,'firstSpacing',h));
        if level==2&&~isempty(baseline),geom=baseline.options;end
        old=ipm_cr.explicit.geometry(geom);
        off=geom;off.radialRedistribution=o.radialRedistribution;off.radialRedistribution.enabled=false;
        disabled=ipm_cr.explicit.geometry(off);
        disabledIdentity=geometry_identity(old,disabled);
        if level==2&&~isempty(baseline)
            report.baselineCheckpoint.disabledMatchesCheckpoint=geometry_identity(baseline.geometry,disabled);
        end
        on=geom;on.radialRedistribution=o.radialRedistribution;
        g=ipm_cr.explicit.geometry(on);
        invariant=prefix_invariants(old,g,o.radialRedistribution.preserveRadius);
        polynomial=exchange_polynomials(g);constant=constant_check(g);
        oldCfl=cfl_test(old);newCfl=cfl_test(g);
        kept=nnz(old.polar.r<=4);reference=cfl_array(old);changed=cfl_array(g);
        completePrefix=kept-3;
        retainedRateError=max(abs(reference(1:completePrefix,:)-changed(1:completePrefix,:)),[],'all');
        retainedRateScale=max(1,max(abs(reference(1:completePrefix,:)),[],'all'));
        edgeRateError=max(abs(reference(completePrefix+1:kept,:)-changed(completePrefix+1:kept,:)),[],'all');
        assert(retainedRateError<=128*eps*retainedRateScale, ...
            'ipm_cr:RedistributionRetainedCfl','CFL at the same preserved physical nodes changed beyond roundoff.');
        [~,i24]=min(abs(old.polar.r-2.4));
        point24=struct('r',old.polar.r(i24),'maximumOldAngularRate',max(reference(i24,:)), ...
            'maximumNewAngularRate',max(changed(i24,:)),'maximumDifference',max(abs(reference(i24,:)-changed(i24,:))));
        fields=struct();
        for name={'shifted_gaussian','transition_gaussian'}
            key=name{1};fields.(key)=analytic_test(g,key);
        end
        row=struct('level',level,'h',h,'requestedNr',nr,'requestedTheta',nt, ...
            'disabledIdentity',disabledIdentity,'invariants',invariant,'polynomial',polynomial,'constant',constant, ...
            'geometry',geometry_summary(g),'oldGeometry',geometry_summary(old), ...
            'oldCfl',oldCfl,'newCfl',newCfl,'sameVelocityMaximumRateRatio',newCfl.rate/oldCfl.rate, ...
            'preservedPointCflRoundoffError',retainedRateError,'strictCflPrefixRows',completePrefix, ...
            'lastThreePreservedRowsCflDifference',edgeRateError, ...
            'lastThreeRowsInterpretation','Their seven-point metric stencil crosses r=4; equality is not part of the prefix contract.', ...
            'preservedR24',point24,'fields',fields);
        report.records=append(report.records,row);report.elapsedSeconds=toc(timer);persist(outputDir,report);
    end
    report.orders=spatial_orders(report.records);
    report.allResolvedOrdersAtLeastFour=all([report.orders.resolvedPairsMeetFourthOrder]);
    report.fourthOrderEvidenceResolvedForAllPairs=all([report.orders.allPairsResolved]);
    report.sourceHashesUnchanged=isequaln(report.sourceInventory,source_inventory(root));
    assert(report.sourceHashesUnchanged,'ipm_cr:RedistributionTestSources','Numerical sources changed during this static test.');
    report.structuralChecksPassed=true;
    report.spatialQualificationPassed=report.allResolvedOrdersAtLeastFour&&report.fourthOrderEvidenceResolvedForAllPairs;
    report.status='completed_static_redistribution_measurements';report.passed=report.spatialQualificationPassed;
    report.orderFailureInterpretation='A measured order below four remains an unresolved qualification failure; this static three-level sequence alone does not identify an implementation error or prove an asymptotic regime.';
catch problem
    report.status='failed';report.errorIdentifier=problem.identifier;
    report.errorMessage=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(timer);persist(outputDir,report);rethrow(problem);
end
report.elapsedSeconds=toc(timer);persist(outputDir,report);
end

function info=geometry_identity(a,b)
keys={'core','polar','exchange','geometryKey','gridOptions','ghostWidth','interpolationWidth','activeDofCount'};
info=struct('allExact',true,'comparisonExcludesBuildSeconds',true);
for k=1:numel(keys)
    key=keys{k};info.(key)=isequaln(a.(key),b.(key));info.allExact=info.allExact&&info.(key);
end
assert(info.allExact,'ipm_cr:RedistributionDisabled','Disabled/default geometry core, polar, exchange and key must be byte-identical apart from buildSeconds.');
end

function info=prefix_invariants(a,b,R)
n=nnz(a.polar.r<=R);assert(n>=7&&numel(b.polar.r)>=n,'ipm_cr:RedistributionPrefix','Retained prefix must fit a seven-point stencil.');
info=struct('retainedRadius',R,'retainedAnnularNodeCount',n, ...
    'radialPrefixExact',isequaln(a.polar.r(1:n),b.polar.r(1:n)), ...
    'coreExact',isequaln(a.core,b.core),'thetaExact',isequaln(a.polar.theta,b.polar.theta), ...
    'innerHaloRadiusExact',isequaln(a.polar.innerHaloRadius,b.polar.innerHaloRadius), ...
    'innerActiveRadiusExact',isequaln(a.polar.innerActiveRadius,b.polar.innerActiveRadius), ...
    'polarGhostCoordinatesExact',isequaln(a.polar.X(a.polar.ghostIndices),b.polar.X(b.polar.ghostIndices))&& ...
        isequaln(a.polar.Y(a.polar.ghostIndices),b.polar.Y(b.polar.ghostIndices)), ...
    'polarFromCoreExact',isequaln(a.exchange.polarFromCore,b.exchange.polarFromCore), ...
    'coreFromPolarPhysicalDonorsExact',isequaln(canonical_map(a),canonical_map(b)), ...
    'completeSevenPointPrefixRows',n-3, ...
    'drMetricPrefixExact',isequaln(a.polar.drMetric(1:n-3),b.polar.drMetric(1:n-3)));
flags={'radialPrefixExact','coreExact','thetaExact','innerHaloRadiusExact','innerActiveRadiusExact', ...
    'polarGhostCoordinatesExact','polarFromCoreExact','coreFromPolarPhysicalDonorsExact','drMetricPrefixExact'};
info.passed=true;for k=1:numel(flags),info.passed=info.passed&&info.(flags{k});end
assert(info.passed,'ipm_cr:RedistributionPrefix','Redistribution changed a retained inner node, core, collar, complete-prefix metric or physical interpolation donor/weight.');
info.productionSpacingPrefixMaximumDifference=max(abs(a.polar.Jr(1:n-3)/(a.polar.nr-1)-b.polar.Jr(1:n-3)/(b.polar.nr-1)));
end

function values=canonical_map(g)
[row,col,w]=find(g.exchange.coreFromPolar);donor=g.exchange.polarDonorIndices(col);
[ir,it]=ind2sub(size(g.polar.X),donor);
values=sortrows([row,g.polar.r(ir),g.polar.theta(it)',w],[1,2,3]);
end

function info=constant_check(g)
f=struct('core',1.25*ones(size(g.core.X)),'polar',1.25*ones(size(g.polar.X)));
s=ipm_cr.explicit.synchronize(f,g);d=ipm_cr.explicit.differentiate(s,g);
value=max([abs(s.core(:)-1.25);abs(s.polar(:)-1.25)]);
derivative=max([abs(d.core.x(:));abs(d.core.y(:));abs(d.polar.x(:));abs(d.polar.y(:))]);
assert(value<=1e-10&&derivative<=1e-10,'ipm_cr:RedistributionConstant','Synchronization and differentiation must preserve a moderate constant.');
info=struct('synchronizationInfinity',value,'gradientInfinity',derivative,'passed',true);
end

function [vx,vy]=velocity(g)
vx=struct('core',.4*g.core.X+1.2,'polar',.4*g.polar.X+1.2);
vy=struct('core',.4*g.core.Y,'polar',.4*g.polar.Y);
end

function info=cfl_test(g)
f=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));[vx,vy]=velocity(g);
[~,rate,position]=ipm_cr.explicit.transport(f,vx,vy,g,struct('reconstruction','linear5','smoothSpeedCoefficient',.1));
info=struct('rate',rate,'position',position,'velocity','vx=.4X+1.2, vy=.4Y; physical wall vy=0', ...
    'allActiveNodesIncluded',true,'reactionRate',0);
end

function rate=cfl_array(g)
[vx,vy]=velocity(g);v=vx.polar;w=vy.polar;w(:,[1,end])=0;
co=g.polar.cosTheta;si=g.polar.sinTheta;sigma=.1*min(g.core.dx,g.core.dy);
vr=v.*co+w.*si;vt=(-v.*si+w.*co)./g.polar.r;vt(:,[1,end])=0;
a=hypot(vt,sigma./g.polar.r);a(:,[1,end])=0;
rate=hypot(vr,sigma)./(g.polar.Jr/(g.polar.nr-1))+a/g.polar.dtheta;
end

function info=analytic_test(g,name)
[a,ax,ay]=gaussian(g.core.X,g.core.Y,name);[b,bx,by]=gaussian(g.polar.X,g.polar.Y,name);
f=ipm_cr.explicit.synchronize(struct('core',a,'polar',b),g);
d=ipm_cr.explicit.differentiate(f,g);[vx,vy]=velocity(g);
[B,~,~]=ipm_cr.explicit.transport(f,vx,vy,g,struct('reconstruction','linear5','smoothSpeedCoefficient',.1));
exactB=struct('core',-vx.core.*ax-vy.core.*ay,'polar',-vx.polar.*bx-vy.polar.*by);
[manual,radial,exactRadial]=manual_polar_transport(f,g,vx,vy,bx,by);
consistency=max(abs(manual-B.polar),[],'all');scale=max(1,max(abs(B.polar),[],'all'));
assert(consistency<=256*eps*scale,'ipm_cr:RedistributionLogicalMetric','Production mapped transport disagrees with direct radial/theta MEX D/J reconstruction.');
info=struct('nativeGradient',error_metrics(struct('core',hypot(d.core.x-ax,d.core.y-ay), ...
    'polar',hypot(d.polar.x-bx,d.polar.y-by)),g), ...
    'transport',error_metrics(struct('core',B.core-exactB.core,'polar',B.polar-exactB.polar),g), ...
    'mappedCentralRadial',polar_metric(radial-exactRadial,g), ...
    'productionLogicalMetricConsistencyInfinity',consistency, ...
    'referenceGradientInfinity',max([hypot(ax(:),ay(:));hypot(bx(:),by(:))]), ...
    'referenceTransportInfinity',max([abs(exactB.core(:));abs(exactB.polar(:))]));
end

function [value,radial,exact]=manual_polar_transport(f,g,vx,vy,bx,by)
p=g.polar;origin=f.core(1,g.core.originIndex);q=f.polar-origin;
[rp,rm]=ipm_cr.explicit.wenoPairSmoothMex(q',1/(p.nr-1),1e-12,1,false(p.ntheta,1),false(p.ntheta,1),0);
[tp,tm]=ipm_cr.explicit.wenoPairSmoothMex(q,p.dtheta,1e-12,1,false(p.nr,1),false(p.nr,1),0);
v=vx.polar;w=vy.polar;w(:,[1,end])=0;co=p.cosTheta;si=p.sinTheta;
vr=v.*co+w.*si;vt=(-v.*si+w.*co)./p.r;vt(:,[1,end])=0;sigma=.1*min(g.core.dx,g.core.dy);
at=hypot(vt,sigma./p.r);at(:,[1,end])=0;
value=-.5*(vr.*(rp'+rm')+hypot(vr,sigma).*(rp'-rm'))./p.Jr-.5*(vt.*(tp+tm)+at.*(tp-tm));
radial=.5*(rp'+rm')./p.Jr;exact=bx.*co+by.*si;
end

function [v,x,y]=gaussian(X,Y,name)
if strcmp(name,'shifted_gaussian'),cx=73;sx=34.5;sy=46;else,cx=12;sx=6;sy=8;end
v=exp(-((X-cx)/sx).^2-(Y/sy).^2);x=-2*(X-cx)/sx^2.*v;y=-2*Y/sy^2.*v;
end

function info=error_metrics(e,g)
masks=struct('global',struct('core',g.core.active,'polar',g.polar.active), ...
    'transition',struct('core',g.core.active&hypot(g.core.X,g.core.Y)>=4&hypot(g.core.X,g.core.Y)<=20, ...
        'polar',g.polar.active&g.polar.r>=4&g.polar.r<=20), ...
    'far',struct('core',false(size(g.core.X)),'polar',g.polar.active&g.polar.r>=20));
info=struct();names=fieldnames(masks);
for j=1:numel(names)
    key=names{j};mask=masks.(key);v=[abs(e.core(mask.core));abs(e.polar(mask.polar))];
    X=[g.core.X(mask.core);g.polar.X(mask.polar)];Y=[g.core.Y(mask.core);g.polar.Y(mask.polar)];
    if isempty(v),entry=struct('absoluteInfinity',NaN,'rms',NaN,'sampleCount',0,'maximumLocation',[NaN,NaN]);
    else,[z,i]=max(v);entry=struct('absoluteInfinity',z,'rms',sqrt(mean(v.^2)),'sampleCount',numel(v),'maximumLocation',[X(i),Y(i)]);end
    info.(key)=entry;
end
end

function info=polar_metric(error,g)
e=struct('core',zeros(size(g.core.X)),'polar',error);info=error_metrics(e,g);
end

function info=geometry_summary(g)
r=g.polar.r;[~,i]=min(abs(r-97));i=max(4,min(numel(r)-3,i));
info=struct('nr',numel(r),'actualFullAxisNr',g.polar.originalNr,'ntheta',numel(g.polar.theta),'requestedGrid',g.gridOptions, ...
    'innerHaloRadius',g.polar.innerHaloRadius,'innerActiveRadius',g.polar.innerActiveRadius, ...
    'activeDofCount',g.activeDofCount,'spacingAt97',struct('r',r(i),'rightGap',r(i+1)-r(i), ...
        'metricSpacing',g.polar.Jr(i)/(numel(r)-1),'sevenPointSpan',r(i+3)-r(i-3)), ...
    'maximumAdjacentGapRatio',max(max(diff(r(2:end))./diff(r(1:end-1))),max(diff(r(1:end-1))./diff(r(2:end)))));
if isfield(g,'radialRedistribution'),info.radialRedistribution=g.radialRedistribution;end
if isfield(g.polar,'redistribution'),info.redistribution=g.polar.redistribution;end
end

function rows=spatial_orders(records)
rows=struct([]);
for field={'shifted_gaussian','transition_gaussian'}
    for quantity={'nativeGradient','transport','mappedCentralRadial'}
        for region={'global','transition','far'}
            e=zeros(1,3);scale=1;
            for j=1:3
                f=records(j).fields.(field{1});e(j)=f.(quantity{1}).(region{1}).absoluteInfinity;
                if strcmp(quantity{1},'transport'),scale=max(scale,f.referenceTransportInfinity);else,scale=max(scale,f.referenceGradientInfinity);end
            end
            floor=1e-11*scale;resolved=e(1:2)>floor&e(2:3)>floor;
            p=log2(e(1:2)./e(2:3));finite=all(isfinite(e))&&all(e>=0);
            % Finite small field differences below the diagnostic floor are
            % reported unresolved; they do not fabricate an asymptotic p.
            passed=finite&&all(p(resolved)>=4);
            rows=append(rows,struct('field',field{1},'quantity',quantity{1},'region',region{1}, ...
                'errors',e,'observedOrders',p,'resolvedPairs',resolved,'resolutionFloor',floor, ...
                'allPairsResolved',all(resolved),'resolvedPairsMeetFourthOrder',passed));
        end
    end
end
end

function [baseline,info]=load_optional_baseline(studyDir)
baseline=[];info=struct('requestedStudyDir',char(studyDir),'available',false,'comparisonPerformed',false);
if isempty(studyDir),info.reason='Explicitly disabled optional historical checkpoint comparison.';return;end
folder=fullfile(char(studyDir),'level_2');cursor=fullfile(folder,'latest.json');
if ~isfile(cursor),info.reason='Optional original level_2 checkpoint is unavailable; same-config disabled/default identity is still required.';return;end
p=jsondecode(fileread(cursor));assert(isfile(p.path)&&strcmp(ipm_c.output.sha256(p.path),p.sha256), ...
    'ipm_cr:RedistributionBaselineSha','Historical accepted checkpoint SHA mismatch.');
data=load(p.path,'checkpoint');c=data.checkpoint;
assert(c.acceptedOnly&&c.state.step==p.sequence&&c.state.scale.canonicalTime==p.tau, ...
    'ipm_cr:RedistributionBaselineState','Historical pointer and complete accepted state disagree.');
baseline=struct('geometry',c.geometry,'options',c.explicitOptions.geometry);
info=struct('requestedStudyDir',char(studyDir),'available',true,'comparisonPerformed',true, ...
    'checkpointPath',p.path,'checkpointSHA256',p.sha256,'tau',p.tau,'step',p.sequence,'authority',c.authority);
end

function inventory=source_inventory(root)
files=[dir(fullfile(root,'test','+ipm_cr','**','*')); ...
    dir(fullfile(root,'time_evolution','+ipm','+mesh','**','*')); ...
    dir(fullfile(root,'time_evolution','+ipm_c','+mesh','**','*'))];paths={mfilename('fullpath')};
if ~endsWith(paths{1},'.m'),paths{1}=[paths{1},'.m'];end
for j=1:numel(files)
    if ~files(j).isdir&&any(endsWith(files(j).name,{'.m','.cpp',['.',mexext]})),paths{end+1}=fullfile(files(j).folder,files(j).name);end %#ok<AGROW>
end
for symbol={'ipm_c.mesh.build','ipm_c.mesh.uniformCoreAxis','ipm.mesh.fdMatrix','ipm_c.output.sha256','ipm_cr.explicit.wenoPairSmoothMex'}
    path=which(symbol{1});assert(~isempty(path)&&isfile(path),'ipm_cr:RedistributionSource','Missing source/native %s.',symbol{1});paths{end+1}=path; %#ok<AGROW>
end
paths=sort(unique(paths));inventory=struct([]);
for j=1:numel(paths),item=dir(paths{j});inventory=append(inventory,struct('path',paths{j},'bytes',item.bytes,'sha256',ipm_c.output.sha256(paths{j})));end
end

function persist(folder,report)
temporary=fullfile(folder,'radial_redistribution.writing.mat');save(temporary,'report','-v7.3');
[ok,message]=movefile(temporary,fullfile(folder,'radial_redistribution.mat'),'f');assert(ok,'ipm_cr:RedistributionTestSave','%s',message);
temporary=fullfile(folder,'radial_redistribution.json.writing');fid=fopen(temporary,'w');assert(fid>=0,'ipm_cr:RedistributionTestSave','Cannot open JSON output.');
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');clear cleanup
[ok,message]=movefile(temporary,fullfile(folder,'radial_redistribution.json'),'f');assert(ok,'ipm_cr:RedistributionTestSave','%s',message);
end

function result=exchange_polynomials(g)
c=g.core;p=g.polar;e=g.exchange;
cx=c.X(e.coreGhostIndices);cy=c.Y(e.coreGhostIndices);
result.coreFromPolar=polynomial_map(e.coreFromPolar, ...
    p.r(p.activeRows),p.theta(p.activeColumns),hypot(cx,cy),atan2(cy,cx), ...
    'r','theta');
result.polarFromCore=polynomial_map(e.polarFromCore, ...
    c.y(c.activeRows),c.x(c.activeColumns),p.Y(e.polarGhostIndices), ...
    p.X(e.polarGhostIndices),'y','x');
end


function info=polynomial_map(M,axis0,axis1,query0,query1,name0,name1)
axis0=axis0(:);axis1=axis1(:);query0=query0(:);query1=query1(:);
nq=size(M,1);assert(size(M,2)==numel(axis0)*numel(axis1), ...
    'ipm_cr:ConnectionTestMapShape','A tensor map must match its active donor packing.');
degrees=(0:6)';tensorPowers=[1,1;2,1;1,2;2,2;3,2;2,3;3,3];
axisErrors=zeros(nq,7,2);tensorErrors=zeros(nq,size(tensorPowers,1));
rowMax=zeros(nq,1);rowAbs=zeros(nq,1);rowNeg=zeros(nq,1);
rowT1=zeros(nq,1);negativeCount=0;zeroSpan=[0,0];
for j=1:nq
    [~,columns,w]=find(M(j,:));w=w(:);columns=columns(:);
    i0=mod(columns-1,numel(axis0))+1;i1=floor((columns-1)/numel(axis0))+1;
    [z0,t0,collapsed0]=normalized(axis0(i0),query0(j));
    [z1,t1,collapsed1]=normalized(axis1(i1),query1(j));
    zeroSpan=zeroSpan+[collapsed0,collapsed1];
    rowMax(j)=max(abs(w));rowAbs(j)=sum(abs(w));
    rowNeg(j)=sum(abs(w(w<0)));negativeCount=negativeCount+nnz(w<0);
    rowT1(j)=abs(sum(w)-1);
    for k=1:7
        degree=degrees(k);
        axisErrors(j,k,1)=abs(sum(w.*z0.^degree)-t0^degree);
        axisErrors(j,k,2)=abs(sum(w.*z1.^degree)-t1^degree);
    end
    for k=1:size(tensorPowers,1)
        powers=tensorPowers(k,:);
        tensorErrors(j,k)=abs(sum(w.*z0.^powers(1).*z1.^powers(2))- ...
            t0^powers(1)*t1^powers(2));
    end
end
perRow=max([max(axisErrors(:,:,1),[],2),max(axisErrors(:,:,2),[],2), ...
    max(tensorErrors,[],2),rowT1],[],2);
assert(all(perRow<=1e-10*max(1,rowAbs)), ...
    'ipm_cr:ConnectionTestPolynomial', ...
    'A local polynomial or constant exchange failed its 1e-10 row-norm-scaled tolerance.');
info=struct('rowCount',nq,'donorCount',size(M,2),'nonzeroCount',nnz(M), ...
    'donorAxes',{{name0,name1}},'degrees',degrees', ...
    'normalization','Each query uses the min/max of its actual nonzero donor support, separately on each coordinate, mapped to [0,1]; an exact-hit collapsed axis uses z=0.', ...
    'collapsedSupportAxisCounts',zeroSpan,'T1InfinityResidual',max(rowT1), ...
    'axis0PolynomialInfinityErrors',max(axisErrors(:,:,1),[],1), ...
    'axis1PolynomialInfinityErrors',max(axisErrors(:,:,2),[],1), ...
    'tensorPowers',tensorPowers,'tensorPolynomialInfinityErrors',max(tensorErrors,[],1), ...
    'maximumRowNormScaledPolynomialError',max(perRow./max(1,rowAbs)), ...
    'operatorInfinityNorm',norm(M,inf),'negativeWeightCount',negativeCount, ...
    'rowsWithNegativeWeights',nnz(rowNeg>0),'maximumRowNegativeAbsoluteMass',max(rowNeg), ...
    'maximumAbsoluteWeight',max(rowMax),'rowMaximumAbsoluteWeight',rowMax, ...
    'rowAbsoluteWeightSum',rowAbs,'rowNegativeAbsoluteWeightSum',rowNeg, ...
    'rowT1AbsoluteResidual',rowT1,'negativeWeightsImplyInstabilityClaimed',false);
end


function [z,t,collapsed]=normalized(nodes,query)
lower=min(nodes);span=max(nodes)-lower;collapsed=span==0;
if collapsed
    assert(abs(query-lower)<=64*eps(max(1,abs(lower))), ...
        'ipm_cr:ConnectionTestExactHit','A collapsed donor axis must be an exact query hit.');
    z=zeros(size(nodes));t=0;
else,z=(nodes-lower)/span;t=(query-lower)/span;end
end


function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
