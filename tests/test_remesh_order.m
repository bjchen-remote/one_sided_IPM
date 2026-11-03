function report=test_remesh_order(outputDir,controls)
%TEST_REMESH_ORDER Nonzero-eta same-time remesh error against direct new data.
% Each level uses one new-grid FEM context. Both complete RHS evaluations
% use that same context; no old-grid RHS, time step or old zero residual is
% used as a reference. Polynomial degree alone is not an order certificate.
if nargin<2,controls=struct();end
assert(isstruct(controls)&&isscalar(controls)&& ...
    all(ismember(fieldnames(controls),{'wallSeconds'})),'ipm_cr:RemeshOrderOptions');
budget=120;if isfield(controls,'wallSeconds'),budget=controls.wallSeconds;end
validateattributes(budget,{'numeric'},{'scalar','real','finite','positive'});
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1||isempty(outputDir),outputDir=fullfile(root,'results','remesh_order');end
outputDir=char(java.io.File(char(outputDir)).getCanonicalPath());
if ~isfolder(outputDir),mkdir(outputDir);end
path=fullfile(outputDir,'remesh_order.json');
assert(~isfile(path),'ipm_cr:RemeshOrderFresh','Use a fresh output directory.');
timer=tic;levels=[12,24,48];fields={'W','q','F','Q'};records=cell(1,numel(levels));
report=struct('kind','complete_same_time_nonpolynomial_remesh_order_v1', ...
    'status','running','passed',false,'targetOrder',4,'minimumObservedOrder',3.9,'levels',levels, ...
    'timeAdvanced',0,'oldGridRhsEvaluated',false,'analyticEtaReference',true, ...
    'etaDefinition','.02*sin(.8*X+.5*Y).*exp(-.3*(X.^2+Y.^2))', ...
    'reference','direct analytic eta on the identical new grid with the same complete FEM/transport/IPC RHS', ...
    'roundoffRule','Errors at or below 32 ulp of max(1,reference field scale) do not supply an observed order.', ...
    'continuumAccuracyCertified',false,'adaptiveTrajectoryQualified',false);
try
    [xc,yc]=meshgrid(linspace(-.75,.75,49),linspace(0,.75,25));
    [tp,rp]=meshgrid(linspace(0,pi,33),unique([linspace(1.2,4,19),logspace(log10(4),10,41)]));
    report.commonRegions=struct('coreX',[-.75,.75],'coreY',[0,.75], ...
        'polarR',[1.2,1e10],'polarTheta',[0,pi], ...
        'coreQueryCount',numel(xc),'polarQueryCount',numel(rp), ...
        'definition','fixed physical queries in each independent patch, active donors only; overlap retained');
    for j=1:numel(levels)
        started=tic;L=levels(j);oldg=ipm_cr.mesh.build(geometry(L));
        a=axes_of(oldg);h=1.2/L;
        a.x=a.x+.2*h*(a.x/1.2).*(1-(a.x/1.2).^2);
        a.y=a.y+.2*h*(a.y/1.2).*(1-a.y/1.2);
        a.dxMetric=metric(a.x,2)';a.dyMetric=metric(a.y,1);
        options=struct('geometry',struct('axes',a),'poissonBackend','fem', ...
            'solverOptions',struct('rOuter',1e10), ...
            'femOptions',struct('degree',5,'quadOrder',8,'stride',[], ...
                'nodalLayout','native','outerBoundary','forced_dtn', ...
                'exteriorModel','finite_source_disk_zero_exterior', ...
                'cornerCorrection',true,'linearSolver','sparse_direct'), ...
            'remesh',struct('enabled',false));
        [direct,c]=ipm_cr.init(options);newg=c.mesh;old=direct;
        for patch={'core','polar'}
            k=patch{1};old.eta.(k)=eta(oldg.(k).X,oldg.(k).Y);
            direct.eta.(k)=eta(newg.(k).X,newg.(k).Y);
            b=ipm_cr.model.background(c.opts.background,old.scale,oldg.(k).X,oldg.(k).Y);
            old.W.(k)=old.eta.(k)+b.H;
            old.eta.(k)(oldg.(k).ghost)=NaN;
        end
        old.ctxId='';old.ctxKey=[];
        [transported,transfer]=ipm_cr.remesh('transfer',old,oldg,newg);
        transported.ctxId=c.id;transported.ctxKey=c.key;
        [transported.rhs,transported.W]=ipm_cr.rhs(transported,c);
        [direct.rhs,direct.W]=ipm_cr.rhs(direct,c);
        assert(isequaln(transported.scale,direct.scale)&& ...
            isequaln(transported.originalAuthority,direct.originalAuthority)&& ...
            transported.step==direct.step,'ipm_cr:RemeshOrderFixedState');
        assert(transfer.routes.core.crossPatch==0&&transfer.routes.polar.crossPatch==0&& ...
            ~transfer.ghostValuesRead,'ipm_cr:RemeshOrderIndependentPatches');
        cm=ipm_cr.mesh.tensorMap(newg.core.y(newg.core.activeRows), ...
            newg.core.x(newg.core.activeColumns),yc(:),xc(:),7);
        pm=ipm_cr.mesh.tensorMap(newg.polar.r(newg.polar.activeRows), ...
            newg.polar.theta,rp(:),tp(:),7);
        maps=struct('core',cm,'polar',pm);
        query=struct('core',struct('X',xc(:),'Y',yc(:)), ...
            'polar',struct('X',rp(:).*cos(tp(:)),'Y',rp(:).*sin(tp(:))));
        stats=c.stats();
        r=struct('L',L,'h',1/L,'seconds',0,'coreSize',size(newg.core.X), ...
            'polarSize',size(newg.polar.X),'activeDofCount',newg.activeDofCount, ...
            'minimumMetricSpacing',[min(newg.core.dxMetric),min(newg.core.dyMetric)], ...
            'maximumNodeMotion',[max(abs(newg.core.x-oldg.core.x)),max(abs(newg.core.y-oldg.core.y))], ...
            'transfer',transfer,'poissonCalls',stats.rhsPoissonCalls);
        for f=1:numel(fields)
            name=fields{f};actual=field(transported,name);reference=field(direct,name);
            for patch={'core','polar'}
                k=patch{1};g=newg.(k);ids=g.activeIndices;
                error=actual.(k)-reference.(k);
                r.(name).native.(k)=measure(error(ids),reference.(k)(ids),g.X(ids),g.Y(ids));
                q=maps.(k);values=reference.(k)(ids);
                r.(name).common.(k)=measure(apply(q,error(ids)),apply(q,values),query.(k).X,query.(k).Y);
            end
            r.(name).native.global=combine(r.(name).native);
            r.(name).common.global=combine(r.(name).common);
        end
        r.seconds=toc(started);records{j}=r;report.records=[records{1:j}];persist();
        fprintf('L=%d dof=%d W %.4e q %.4e F %.4e Q %.4e seconds %.3f\n',L, ...
            newg.activeDofCount,r.W.native.global.absoluteInfinity,r.q.native.global.absoluteInfinity, ...
            r.F.native.global.absoluteInfinity,r.Q.native.global.absoluteInfinity,r.seconds);
        clear c oldg newg old transported direct cm pm maps
        assert(toc(timer)<=budget,'ipm_cr:RemeshOrderBudget','The independent remesh order test exceeded its budget.');
    end
    qualified=true;
    for f=1:numel(fields)
        name=fields{f};
        for scope={'native','common'}
            key=scope{1};errors=zeros(1,numel(levels));floors=errors;
            for j=1:numel(levels)
                m=records{j}.(name).(key).global;
                errors(j)=m.absoluteInfinity;floors(j)=m.roundoffScale;
            end
            resolved=errors(1:end-1)>floors(1:end-1)&errors(2:end)>floors(2:end);
            observed=NaN(1,numel(levels)-1);
            observed(resolved)=log(errors(find(resolved))./errors(find(resolved)+1))/log(2);
            passed=any(resolved)&&all(observed(resolved)>=report.minimumObservedOrder);
            report.order.(name).(key)=struct('absoluteErrors',errors,'roundoffScales',floors, ...
                'resolvedPairs',resolved,'observedOrders',observed,'orderQualified',passed);
            qualified=qualified&&passed;
        end
    end
    report.orderQualified=qualified;persist();
    assert(qualified,'ipm_cr:RemeshOrder', ...
        'Complete remesh W/q/F/Q needs near-fourth order (>=3.9) on resolved native and common-query pairs; see remesh_order.json.');
    report.status='passed';report.passed=true;
catch problem
    report.status='failed';report.failure=struct('identifier',problem.identifier, ...
        'message',problem.message,'stack',problem.stack);persist();rethrow(problem)
end
persist();
    function persist()
        report.elapsedSeconds=toc(timer);ipm_cr.io.json(path,report);
    end
end

function geo=geometry(L)
a=1.2;geo=struct('coreHalfWidth',a,'coreHeight',a,'coreSpacing',a/L, ...
    'polarStart',.65,'ghostWidth',3,'interpolationWidth',7, ...
    'grid',struct('nr',4*L+1,'ntheta',2*L+1,'rmax',1e10,'firstSpacing',a/(2*L), ...
    'radialMap','uniform_core','coreRadius',2,'coreTailFraction',.01,'stencilWidth',7), ...
    'coreMap',struct('enabled',true,'intervals',[2*L,L],'minSpacing',[a/L,a/L]/sinh(1)));
end

function a=axes_of(g)
a=struct('x',g.core.x,'y',g.core.y,'r',g.polar.r,'theta',g.polar.theta, ...
    'dxMetric',g.core.dxMetric,'dyMetric',g.core.dyMetric,'drMetric',g.polar.drMetric);
end

function value=metric(axis,logicalLength)
value=ipm_cr.mesh.fd(linspace(0,logicalLength,numel(axis)),1,7)*axis(:)*logicalLength/(numel(axis)-1);
end

function v=eta(x,y)
v=.02*sin(.8*x+.5*y).*exp(-.3*(x.^2+y.^2));
end

function v=field(s,name)
if strcmp(name,'W'),v=s.W;else,v=s.rhs.(name);end
end

function v=apply(M,donor)
[row,column,weight]=find(M);[~,anchor]=max(abs(M),[],2);anchor=full(anchor);
base=donor(anchor);v=base+accumarray(row,weight.*(donor(column)-base(row)),[size(M,1),1]);
end

function m=measure(error,reference,x,y)
assert(all(isfinite(error))&&all(isfinite(reference)),'ipm_cr:RemeshOrderFinite');
[value,index]=max(abs(error));scale=max(abs(reference));
m=struct('absoluteInfinity',value,'absoluteRms',sqrt(mean(error.^2)), ...
    'referenceInfinity',scale,'mixedRelativeInfinity',value/max(1,scale), ...
    'maximumLocation',struct('X',x(index),'Y',y(index)), ...
    'roundoffScale',32*eps(max(1,scale)),'pointCount',numel(error));
end

function m=combine(parts)
if parts.core.absoluteInfinity>=parts.polar.absoluteInfinity,m=parts.core;m.patch='core';
else,m=parts.polar;m.patch='polar';end
m.referenceInfinity=max(parts.core.referenceInfinity,parts.polar.referenceInfinity);
m.mixedRelativeInfinity=m.absoluteInfinity/max(1,m.referenceInfinity);
m.roundoffScale=max(parts.core.roundoffScale,parts.polar.roundoffScale);
m.pointCount=parts.core.pointCount+parts.polar.pointCount;
m.absoluteRms=sqrt((parts.core.absoluteRms^2*parts.core.pointCount+ ...
    parts.polar.absoluteRms^2*parts.polar.pointCount)/m.pointCount);
end
