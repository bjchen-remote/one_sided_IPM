function report=check_shared_spatial_sdirk(beforeStepFile)
%CHECK_SHARED_SPATIAL_SDIRK Shared SDIRK arithmetic and unique-origin manufacture.
% The oracle is the immutable old stepRk with only its native-flow entry
% routed to this test's manufactured equation. No physical PDE is advanced.
if nargin<1
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    beforeStepFile=fullfile(root,'runs','ipm_c_ipm_architecture_20261020_01', ...
        'before_runtime','+ipm','+evolve','stepRk.m');
end
assert(isfile(beforeStepFile),'ipm:SharedSdirkOracle','The immutable old stepRk source is required.');
temporary=tempname;package=fullfile(temporary,'+ipm_sdirk_oracle');mkdir(package);
currentPackage=fullfile(temporary,'+ipm_sdirk_current');mkdir(currentPackage);
source=fileread(beforeStepFile);
source=strrep(source,'ipm.evolve.flow(','ipm_sdirk_oracle.flow(');
write_text(fullfile(package,'stepRk.m'),source);
write_text(fullfile(package,'flow.m'),sprintf([ ...
    'function [F,flow]=flow(rho,ops,scale)\n', ...
    '[F,flow]=ops.manufacturedEvaluate(rho,ops,scale);\nend\n']));
current=fileread(which('ipm.evolve.stepRk'));
current=strrep(current,'ipm.evolve.flow(','ipm_sdirk_current.flow(');
write_text(fullfile(currentPackage,'stepRk.m'),current);
write_text(fullfile(currentPackage,'flow.m'),fileread(fullfile(package,'flow.m')));
addpath(temporary);cleanup=onCleanup(@()remove_oracle(temporary)); %#ok<NASGU>
method='sdirk2_alexander_practical_v2';contract=ipm.evolve.sdirkContract(method);
ops=fixture_ops(contract);rho=reshape(linspace(.1,.9,20),4,5);
scale=ipm.evolve.initialScale(ops);h=.05;
[reference,referenceFlow,referenceScale,~,referenceInfo]= ...
    ipm_sdirk_oracle.stepRk(rho,h,ops,method,scale);
[fallback,fallbackFlow,fallbackScale,~,fallbackInfo]= ...
    ipm_sdirk_current.stepRk(rho,h,ops,method,scale);
assert(isequaln(reference,fallback)&&isequaln(referenceFlow,fallbackFlow)&& ...
    isequaln(referenceScale,fallbackScale)&&isequaln(referenceInfo,fallbackInfo), ...
    'ipm:SharedSdirkCartesianFallback','The Cartesian fallback must preserve old arithmetic exactly.');
hooked=ops;hooked.spatialGeometryKey=struct('kind','manufactured_cartesian_identity', ...
    'x',ops.x,'y',ops.y,'operatorVersion',1,'remeshCount',0);
hooked.spatialOperators=struct('evaluate',@manufactured, ...
    'derivativeX',@(R)R*ops.Dx','pack',@(R)R(:),'unpack',@(x)reshape(x,size(rho)), ...
    'preconditioner',@(flow,scales,gh,c)diagonal_preconditioner(flow,scales,gh,c,@(R)R(:)));
[identity,identityFlow,identityScale,~,identityInfo]=ipm.evolve.stepRk(rho,h,hooked,method,scale);
assert(isequaln(reference,identity)&&isequaln(referenceFlow,identityFlow)&& ...
    isequaln(referenceScale,identityScale)&&isequaln(referenceInfo,identityInfo), ...
    'ipm:SharedSdirkIdentityHooks','Identity spatial hooks must preserve the old SDIRK stages.');

nr=5;nt=7;origin=0.3;positive=reshape(linspace(.1,.9,(nr-1)*nt),nr-1,nt);
initial=[repmat(origin,1,nt);positive];polar=rmfield(ops,{'x','y','Dx'});
pack=@(R)pack_origin(R,nr,nt);unpack=@(x)unpack_origin(x,nr,nt);
polar.spatialGeometryKey=struct('kind','manufactured_unique_origin', ...
    'r',linspace(0,1,nr)','theta',linspace(0,pi,nt), ...
    'operatorVersion',1,'remeshCount',0);
polar.spatialOperators=struct('evaluate',@manufactured,'derivativeX',@(R)R, ...
    'pack',pack,'unpack',unpack, ...
    'preconditioner',@(flow,scales,gh,c)diagonal_preconditioner(flow,scales,gh,c,pack));
assert(numel(pack(initial))==1+(nr-1)*nt&&isequal(unpack(pack(initial)),initial));
forced=polar;forced.timeIntegrationContract.preconditionerDefectTrials=0;
[jvValue,~,~,~,jvInfo]=ipm.evolve.stepRk(initial,.05,forced,method,scale);
assert(all(jvValue(1,:)==jvValue(1,1))&& ...
    all([jvInfo.implicit.stages.completeJvGmresUpdates]>=1), ...
    'ipm:SharedSdirkPackedJv','The unique-origin manufacture must also exercise full central Jv GMRES.');
counts=[4,8,16];errors=zeros(1,3);clockErrors=errors;records=cell(1,3);T=.4;
for level=1:numel(counts)
    value=initial;s=ipm.evolve.initialScale(polar);cache=[];cacheHits=0;
    for step=1:counts(level)
        [value,~,s,cache,info]=ipm.evolve.stepRk(value,T/counts(level),polar,method,s,cache);
        assert(all(value(1,:)==value(1,1)),'ipm:SharedSdirkUniqueOrigin');
        cacheHits=cacheHits+double(info.initialStageCacheHit);
    end
    errors(level)=max(abs(value-exp(-.7*T)*initial),[],'all');
    clockErrors(level)=abs(s.physicalTime-expm1(-.7*T)/(-.7));
    assert(cacheHits==counts(level)-1,'ipm:SharedSdirkAcceptedCache');
    records{level}=struct('substeps',counts(level),'densityError',errors(level), ...
        'physicalClockError',clockErrors(level),'initialStageCacheHits',cacheHits);
end
orders=log2(errors(1:end-1)./errors(2:end));
clockOrders=log2(clockErrors(1:end-1)./clockErrors(2:end));
assert(all(orders>1.8)&all(orders<2.2)&all(clockOrders>1.8)&all(clockOrders<2.2), ...
    'ipm:SharedSdirkOrder','The common two-stage method must retain second-order density and clock convergence.');
changed=polar;changed.spatialGeometryKey.r(2)=changed.spatialGeometryKey.r(2)+.01;
[~,~,~,~,keyMiss]=ipm.evolve.stepRk(value,.025,changed,method,s,cache);
changed=polar;changed.spatialOperators.derivativeX=@(R)2*R;
[~,~,~,~,operatorMiss]=ipm.evolve.stepRk(value,.025,changed,method,s,cache);
changed=polar;changed.remeshCount=1;
[~,~,~,~,remeshMiss]=ipm.evolve.stepRk(value,.025,changed,method,s,cache);
assert(~keyMiss.initialStageCacheHit&&~operatorMiss.initialStageCacheHit&& ...
    ~remeshMiss.initialStageCacheHit,'ipm:SharedSdirkCacheBinding', ...
    'Geometry, spatial operator identity and remesh changes must invalidate the initial cache.');
report=struct('kind','shared_sdirk_cartesian_oracle_and_unique_origin_manufacture_v1', ...
    'passed',true,'oracleStepFile',beforeStepFile,'cartesianFallbackBitwise',true, ...
    'identityHooksBitwise',true,'uniqueOriginDimension',numel(pack(initial)), ...
    'uniqueOriginCentralJvStages',jvInfo.implicit.stages, ...
    'contract',contract,'records',{records},'densityOrders',orders,'physicalClockOrders',clockOrders, ...
    'geometryOperatorAndRemeshCacheControls',true,'physicalPdeAdvanced',false);
end

function ops=fixture_ops(contract)
ops=struct('dynamicScaleGeometry','isotropic','remeshCount',0, ...
    'x',linspace(-1,1,5),'y',linspace(0,1,4)','Dx',speye(5), ...
    'timeIntegrationContract',contract,'manufacturedEvaluate',@manufactured);
ops.rescaling=struct('enabled',true,'lengthGauge','wall_gradient_integral', ...
    'cOmegaGauge','wall_omega_global_spline_maximum','gaugeRateGuess',[0,0]);
end
function [F,flow]=manufactured(rho,~,~)
F=-.7*rho;
flow=struct('c_l',0,'c_omega',-.7,'c_r',0, ...
    'transportU1',zeros(size(rho)),'transportU2',zeros(size(rho)),'generatorCflRate',0);
end
function [L,U,rowBalance]=diagonal_preconditioner(flow,scales,gh,c,pack)
v=pack(scales);n=numel(v);B=flow.c_omega*speye(n);M=speye(n)-gh*B;
M=M*spdiags(v,0,n,n);rowBalance=1./max(sum(abs(M),2),realmin);
[L,U]=ilu(spdiags(rowBalance,0,n,n)*M,struct('type','ilutp','droptol',c.iluDropTolerance));
end
function x=pack_origin(R,nr,nt)
assert(isequal(size(R),[nr,nt])&&all(R(1,:)==R(1,1)),'ipm:SharedSdirkOriginPack');
positive=R(2:end,:).';x=[R(1,1);positive(:)];
end
function R=unpack_origin(x,nr,nt)
assert(isequal(size(x),[1+(nr-1)*nt,1]),'ipm:SharedSdirkOriginUnpack');
R=[repmat(x(1),1,nt);reshape(x(2:end),nt,nr-1).'];
end
function write_text(path,value)
fid=fopen(path,'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',value);
end
function remove_oracle(path)
clear('ipm_sdirk_oracle.stepRk','ipm_sdirk_oracle.flow', ...
    'ipm_sdirk_current.stepRk','ipm_sdirk_current.flow');
rmpath(path);rmdir(path,'s');
end
