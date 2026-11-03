function report=test_mainline_mesh(outputDir)
%TEST_MAINLINE_MESH Exact nonuniform half-density and active-only exchange.
% Geometry, exchange and analytic scalar-transport checks; no Poisson,
% gauge, PDE step or remesh call.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'result','mainline_mesh');end
if ~isfolder(outputDir),mkdir(outputDir);end
m=ipm_cr.mesh.build(struct('preset','medium'));
c=ipm_cr.mesh.build(struct('preset','coarse','grid',struct('rmax',1e10)));
assert(isequal(size(m.core.X),[129,241])&&isequal(size(m.polar.X),[549,193]));
assert(isequal(size(c.core.X),[65,121])&&isequal(size(c.polar.X),[275,97]));
assert(m.activeDofCount==134988&&c.activeDofCount==33514);
assert(isequal(c.core.x,m.core.x(1:2:end))&&isequal(c.core.y,m.core.y(1:2:end))&& ...
    isequal(c.polar.r,m.polar.r(1:2:end))&&isequal(c.polar.theta,m.polar.theta(1:2:end)), ...
    'ipm_cr:MeshExactHalf','Every coarse axis must be an exact medium-node subsequence.');
assert(all(c.core.dxMetric==2*m.core.dxMetric(1:2:end))&& ...
    all(c.core.dyMetric==2*m.core.dyMetric(1:2:end))&& ...
    all(c.polar.drMetric==2*m.polar.drMetric(1:2:end)));
assert(c.core.dxMetric(c.core.originIndex)==m.core.dxMetric(m.core.originIndex)*2&& ...
    c.core.mapParameters.logicalSpacing(1)==2*m.core.mapParameters.logicalSpacing(1));
assert(max(diff(m.core.x))/min(diff(m.core.x))>10&& ...
    max(diff(m.core.y))/min(diff(m.core.y))>10);
records=struct([]);
for j=1:2
    if j==1,g=m;else,g=c;end
    [a,item]=check_geometry(g);
    if isempty(records),records=repmat(item,1,2);end
    records(j)=item;
    replay=ipm_cr.mesh.build(struct('axes',a));
    assert(isequaln(replay.geometryKey,g.geometryKey), ...
        'ipm_cr:MeshExplicitReplay','Explicit axes must reconstruct the same mathematical geometry.');
end
a.dxMetric(1)=1.01*a.dxMetric(1);
changed=ipm_cr.mesh.build(struct('axes',a));
assert(~isequaln(changed.geometryKey,c.geometryKey),'ipm_cr:MeshMetricIdentity', ...
    'A changed metric must invalidate the mathematical geometry identity.');
% Retain the public uniform and analytic-sinh generator paths.
legacy=ipm_cr.mesh.build(struct());
assert(~isfield(legacy.core,'coreMap')&&legacy.core.x(legacy.core.originIndex)==0);
mapped=ipm_cr.mesh.build(struct('coreMap',struct('enabled',true,'kind','sinh', ...
    'intervals',[80,40],'minSpacing',[.02,.02])));
assert(mapped.core.coreMap.enabled&&all(mapped.core.dxMetric>0));
expect(@()ipm_cr.mesh.build(struct('preset','coarse','coreSpacing',.06)));
expect(@()ipm_cr.mesh.build(struct('preset','coarse','grid',struct('rmax',1e9))));
bad=a;bad.drMetric(1)=-1;expect(@()ipm_cr.mesh.build(struct('axes',bad)));
bad=a;bad.theta(2)=bad.theta(2)*1.1;expect(@()ipm_cr.mesh.build(struct('axes',bad)));
bad=rmfield(a,'dyMetric');expect(@()ipm_cr.mesh.build(struct('axes',bad)));
report=struct('kind','ipm_cr_nonuniform_exact_half_geometry_v1','passed',true, ...
    'sourceGeometrySHA256',m.axisSource.sourceGeometrySHA256, ...
    'exactCoordinateHalving',true,'metricsScaleWithLogicalStep',true, ...
    'explicitAxesReplay',true,'metricIdentityEnforced',true, ...
    'legacyConstructorsRetained',true,'records',records, ...
    'interfaceDerivatives',interface_derivatives(c,m), ...
    'timeAdvanced',0,'poissonCalls',0,'trajectoryQualified',false);
ipm_cr.io.json(fullfile(outputDir,'mainline_mesh.json'),report);
end

function [a,info]=check_geometry(g)
c=g.core;p=g.polar;e=g.exchange;
assert(c.x(c.originIndex)==0&&c.y(1)==0&&p.r(1)>0&&p.r(end)==1e10);
assert(~any(c.ghost(e.coreDonorIndices))&&~any(p.ghost(e.polarDonorIndices)));
assert(isequal(c.dxMetric,c.Jx*c.dxi)&&isequal(c.dyMetric,c.Jy*c.deta));
assert(c.mapParameters.originMetricSpacing(1)==c.dxMetric(c.originIndex)&& ...
    c.mapParameters.originMetricSpacing(2)==c.dyMetric(1));
constant=struct('core',3*ones(size(c.X)),'polar',3*ones(size(p.X)));
assert(isequal(ipm_cr.mesh.fill(constant,g),constant));
localZero=struct('core',zeros(size(c.X)),'polar',zeros(size(p.X)));
localZero.polar(:,1:5)=1;
angles=(.2:.031:.8)*pi;radius=.99*p.r(end);
[zero,map]=ipm_cr.mesh.sample(localZero,g,radius*cos(angles),radius*sin(angles));
assert(~any(zero,'all')&&isequal(zero,ipm_cr.mesh.sample(localZero,map)), ...
    'ipm_cr:MeshLocalAnchor','A remote nonzero source contaminated an exact local zero stencil.');
map.referenceOnly=true;assert(isequal(zero,ipm_cr.mesh.sample(localZero,map)));
[one,scalarMap]=ipm_cr.mesh.sample(constant,g,.13,.17);
scalarMap.referenceOnly=true;
assert(one==3&&ipm_cr.mesh.sample(constant,scalarMap)==3, ...
    'ipm_cr:MeshScalarReference','One-query native/reference interpolation must both preserve constants.');
native=struct('core',sin(c.X)+.3*c.Y,'polar',sin(p.X)+.3*p.Y./(1+p.Y));
filled=ipm_cr.mesh.fill(native,g);poison=native;
poison.core(c.ghost)=1e12;poison.polar(p.ghost)=-1e12;
assert(isequal(filled,ipm_cr.mesh.fill(poison,g))&& ...
    isequal(filled,ipm_cr.mesh.fill(filled,g))&& ...
    isequal(filled.core(c.active),native.core(c.active))&& ...
    isequal(filled.polar(p.active),native.polar(p.active)));
% D7 physical derivatives reproduce low polynomials on the actual axes.
axisError=zeros(1,4);axes={c.x,c.y,p.r,p.theta};operators={c.Dx,c.Dy,p.Dr,p.Dt};
for j=1:4
    x=axes{j}(:);D=operators{j};span=x(end)-x(1);z=(x-x(1))/span;
    for k=1:6
        exact=k*z.^(k-1)/span;
        error=max(abs(D*z.^k-exact))/max(1,norm(D,inf));
        axisError(j)=max(axisError(j),error);
    end
end
assert(all(axisError<1e-10),'ipm_cr:MeshPolynomialDerivative','D7 polynomial reproduction failed.');
% Physical q and Q derivatives must not depend on algebraic ghost values.
d=ipm_cr.mesh.diff(native,g);dp=ipm_cr.mesh.diff(poison,g);
for patch={'core','polar'}
    key=patch{1};mask=g.(key).active;
    for direction={'x','y'}
        name=direction{1};
        assert(isequal(d.(key).(name)(mask),dp.(key).(name)(mask)), ...
            'ipm_cr:MeshActiveDerivative','Artificial ghost values changed an active derivative.');
    end
end
a=struct('x',c.x,'y',c.y,'r',p.r,'theta',p.theta, ...
    'dxMetric',c.dxMetric,'dyMetric',c.dyMetric,'drMetric',p.drMetric);
info=struct('coreShape',size(c.X),'polarShape',size(p.X), ...
    'fullNodes',numel(c.X)+numel(p.X),'activeNodes',g.activeDofCount, ...
    'activeDonorsOnly',true,'constantExchangeExact',true, ...
    'ghostPoisonIndependent',true,'activeValuesPreserved',true, ...
    'derivativeRowNormScaledError',axisError);
end

function expect(fn)
failed=false;
try
    fn();
catch
    failed=true;
end
assert(failed,'ipm_cr:MeshExpectedFailure','Invalid geometry was accepted.');
end

function result=interface_derivatives(coarse,medium)
% Continuous analytic values, gradients and Hessians supply the reference.
% The same true scalar field is sampled independently on each actual grid.
parameters=[-.25,.35,.65,.75,.2;3.9,.85,2.8,1.25,-.2];
records=struct([]);
for fixture=1:2
    errors=zeros(2,2);Berrors=zeros(2,1);
    for level=1:2
        if level==1,g=coarse;else,g=medium;end
        f=struct();vx=f;vy=f;truth=f;reference=f;
        p=parameters(fixture,:);
        for patch={'core','polar'}
            key=patch{1};x=g.(key).X;y=g.(key).Y;
            u=(x-p(1))/p(3);v=(y-p(2))/p(4);
            z=exp(-(u.^2+2*p(5)*u.*v+v.^2));
            a=-2*(u+p(5)*v)/p(3);b=-2*(v+p(5)*u)/p(4);
            fx=a.*z;fy=b.*z;fxx=(a.^2-2/p(3)^2).*z;
            fxy=(a.*b-2*p(5)/(p(3)*p(4))).*z;
            f.(key)=z;vx.(key)=.8+.1*x;vy.(key)=.15*y;
            truth.(key)=-vx.(key).*fx-vy.(key).*fy;
            reference.(key)=-.1*fx-vx.(key).*fxx-vy.(key).*fxy;
        end
        f=ipm_cr.mesh.fill(f,g);
        op=ipm_cr.op.prepare(f,g,struct('reconstruction','linear5','smoothSpeedCoefficient',.1));
        out=op(vx,vy,struct('cfl',false,'partials',false));
        derivative=ipm_cr.mesh.diff(ipm_cr.mesh.fill(out.B,g),g);
        for patch={'core','polar'}
            key=patch{1};part=g.(key);band=false(size(part.X));
            if strcmp(key,'core')
                ix=part.activeColumns;iy=part.activeRows;
                band(:,ix([1:3,end-2:end]))=true;band(iy(end-2:end),:)=true;
            else,band(part.activeRows(1:3),:)=true;end
            band=band&part.active;
            e=abs(derivative.(key).x-reference.(key));
            errors(level,1)=max(errors(level,1),max(e(part.active)));
            errors(level,2)=max(errors(level,2),max(e(band)));
            e=abs(out.B.(key)-truth.(key));Berrors(level)=max(Berrors(level),max(e(part.active)));
        end
    end
    orders=log2(errors(1,:)./errors(2,:));
    assert(all(orders>=4),'ipm_cr:MeshInterfaceOrder', ...
        'The actual-grid scalar B_x full/interface finite-pair orders must remain >=4.');
    row=struct('parameters',parameters(fixture,:),'levels',{{'coarse','medium'}}, ...
        'BxErrors',errors,'BxOrders',orders,'BErrors',Berrors, ...
        'scope','All-active and moving three-layer interface bands; analytic prescribed velocity, not full IPM Q.');
    if isempty(records),records=row;else,records(end+1)=row;end %#ok<AGROW>
end
result=struct('passed',true,'records',records,'minimumOrder',4, ...
    'physicalDerivativeClosure','same_patch_active_D7_v1','completeSolverQualified',false);
end
