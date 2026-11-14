classdef diagnostics
    %DIAGNOSTICS Passive field readout and overlap diagnostics.
    methods (Static)
        function values=fields(state,g,method,x,y,spline)
            %FIELDS Read complete accepted fields; no comparison-region policy.
            assert(isnumeric(x)&&isreal(x)&&~isempty(x)&&isnumeric(y)&&isreal(y)&& ...
                isequal(size(x),size(y))&&all(isfinite(x),'all')&&all(isfinite(y),'all'), ...
                'ipm_cr:SampleQueries','Use matching nonempty finite real coordinate arrays.');
            assert(isfield(state.rhs,'diagnosticsComplete')&&state.rhs.diagnosticsComplete, ...
                'ipm_cr:SampleDiagnostics','Sample a complete accepted endpoint.');
            if ~isfield(state,'eta')
                f=state.rhs;u1=struct('core',f.physical.core.u1,'polar',f.physical.polar.u1);
                u2=struct('core',f.physical.core.u2,'polar',f.physical.polar.u2);
                values=struct('W',ipm_cr.mesh.interpolate.sample(state.W,g,x,y), ...
                    'q',ipm_cr.mesh.interpolate.sample(f.q,g,x,y),'Q',ipm_cr.mesh.interpolate.sample(f.Q,g,x,y), ...
                    'F',ipm_cr.mesh.interpolate.sample(f.F,g,x,y), ...
                    'u1Centered',ipm_cr.mesh.interpolate.sample(u1,g,x,y)-f.Uref, ...
                    'u2',ipm_cr.mesh.interpolate.sample(u2,g,x,y),'rates',rate_values(state), ...
                    'scalars',scalar_values(state));return;
            end
            if nargin<6
                spline=struct('x',ipm_cr.space.basis.build(g.core.x), ...
                    'y',ipm_cr.space.basis.build(g.core.y), ...
                    't',ipm_cr.space.basis.build(log(g.polar.r)), ...
                    'theta',ipm_cr.space.basis.build(g.polar.theta));
            end
            b=ipm_cr.model.background(method.background,state.scale,x,y);
            eta=ipm_cr.mesh.interpolate.sample(state.eta,g,x,y);
            d=ipm_cr.space.basis.gradient(state.eta,g,spline);
            qr=ipm_cr.mesh.exchange(struct('core',d.core.x,'polar',d.polar.x),g);
            d=ipm_cr.space.basis.gradient(state.rhs.residualRhs,g,spline);
            Qr=ipm_cr.mesh.exchange(struct('core',d.core.x,'polar',d.polar.x),g);
            f=state.rhs;Ht=f.c_omega*b.H-(f.c_l*x+f.c_r).*b.x-f.c_l*y.*b.y+f.c_r*b.originX;
            Htx=(f.c_omega-f.c_l)*b.x-(f.c_l*x+f.c_r).*b.xx-f.c_l*y.*b.xy;
            % These explicit-background time terms already include dt_physical/dtau.
            Ht=Ht+b.timeH;Htx=Htx+b.timeX;
            u1=struct('core',f.physical.core.u1,'polar',f.physical.polar.u1);
            u2=struct('core',f.physical.core.u2,'polar',f.physical.polar.u2);
            values=struct('W',b.H+eta,'q',b.x+ipm_cr.mesh.interpolate.sample(qr,g,x,y), ...
                'Q',Htx+ipm_cr.mesh.interpolate.sample(Qr,g,x,y), ...
                'F',Ht+ipm_cr.mesh.interpolate.sample(f.residualRhs,g,x,y), ...
                'u1Centered',ipm_cr.mesh.interpolate.sample(u1,g,x,y)-f.Uref, ...
                'u2',ipm_cr.mesh.interpolate.sample(u2,g,x,y), ...
                'rates',rate_values(state),'scalars',scalar_values(state));
        end

        function info=overlap(state,g,options)
            %AUDIT Independent duplicate-patch audit for an accepted explicit state.
            % This offline/checkpoint diagnostic never participates in the spatial RHS,
            % changes a field, advances time, or declares a prescribed accuracy passed.
            % Each patch uses ONLY its own ACTIVE donor subarray and a local 7-by-7
            % tensor polynomial. No unique-owner interpolation or ghost values are read.
            % Queries are inside BOTH active domains. Near their boundaries the seven
            % available active nodes provide a one-sided stencil; this is reported.
            %
            % OPTIONS: radialSamples=11, angularSamples=65, includeWall=true,
            % includeQueryCoordinates=false. All interior ray/radius samples use fixed
            % fractional offsets, independently of either patch's collocation nodes.

            if nargin<3,options=struct();end
            assert(isstruct(state)&&isscalar(state)&&all(isfield(state,{'W','rhs'}))&& ...
                isstruct(g)&&isscalar(g)&&all(isfield(g,{'core','polar'})), ...
                'ipm_cr:ExplicitAuditState','Supply one complete accepted explicit state and its geometry.');
            opts=resolve_options(options);
            c=g.core;p=g.polar;cr=c.activeRows;cc=c.activeColumns;
            staircase=isfield(c,'stair');
            pr=p.activeRows;pc=p.activeColumns;
            cx=c.x(cc);cy=c.y(cr);rr=p.r(pr);tt=p.theta(pc);
            assert((staircase||all(c.active(cr,cc),'all'))&&all(p.active(pr,pc),'all')&& ...
                min([numel(cx),numel(cy),numel(rr),numel(tt)])>=7, ...
                'ipm_cr:ExplicitAuditDonors','Both independent active donor tensors must fit seven-point stencils.');
            key=struct('coreX',cx(:),'coreY',cy(:),'polarR',rr(:), ...
                'polarTheta',tt(:),'radialSamples',opts.radialSamples, ...
                'angularSamples',opts.angularSamples,'includeWall',opts.includeWall);
            if staircase,key.stair=struct('edges',c.stair.edges,'top',c.stair.top,'active',c.active);end
            persistent retained
            hit=[];
            for j=1:numel(retained)
                if isequaln(retained{j}.key,key),hit=j;break;end
            end
            if isempty(hit)
                cache=prepare_maps(key,c);
                if isempty(retained),retained={cache};else,retained{end+1}=cache;end
                if numel(retained)>2,retained=retained(end-1:end);end
                cacheReused=false;
            else
                cache=retained{hit};cacheReused=true;
            end
            flow=state.rhs;
            assert(all(isfield(flow,{'q','F','Q','physical'}))&& ...
                all(isfield(flow.physical,{'core','polar'}))&& ...
                all(isfield(flow.physical.core,{'u1','u2'}))&& ...
                all(isfield(flow.physical.polar,{'u1','u2'})), ...
                'ipm_cr:ExplicitAuditFields','The accepted flow must contain q, F, Q and both physical velocities.');
            fields=struct();
            [fields.rho,~]=compare(state.W);
            [fields.q,~]=compare(flow.q);
            [fields.F,~]=compare(flow.F);
            [fields.Q,~]=compare(flow.Q);
            u1=struct('core',flow.physical.core.u1,'polar',flow.physical.polar.u1);
            u2=struct('core',flow.physical.core.u2,'polar',flow.physical.polar.u2);
            [fields.u1,velocity1]=compare(u1);
            [fields.u2,velocity2]=compare(u2);
            coreSpeed=hypot(velocity1.core,velocity2.core);
            polarSpeed=hypot(velocity1.polar,velocity2.polar);
            velocityError=hypot(velocity1.core-velocity1.polar,velocity2.core-velocity2.polar);
            [velocityMaximum,velocityIndex]=max(velocityError);
            velocity=norm_record(velocityMaximum,max(coreSpeed),max(polarSpeed));
            velocity.absoluteRms=sqrt(mean(velocityError.^2));
            velocity.maximumLocation=query_location(cache,velocityIndex);
            velocity.definition='Euclidean difference of the two physical Cartesian velocity samples';
            info=struct('kind','independent_active_overlap_local7_checkpoint_audit_v1', ...
                'timeAdvanced',0,'fieldsModified',false,'fields',fields,'physicalVelocity',velocity, ...
                'queryCount',numel(cache.x),'interiorQueryCount',nnz(~cache.wall), ...
                'wallQueryCount',nnz(cache.wall),'queryBounds',cache.bounds, ...
                'activeDonorBounds',cache.donorBounds,'interpolationWidth',7, ...
                'queryPattern','fractionally offset log-radius rays plus independent wall radii', ...
                'stencilPolicy','seven available ACTIVE nodes; one-sided near active boundaries', ...
                'oneSidedQueryCounts',cache.oneSidedCounts, ...
                'uniqueOwnerInterpolationUsed',false,'ghostValuesRead',false, ...
                'geometryQueryMapCacheReused',cacheReused,'fieldsCached',false, ...
                'acceptanceThresholdApplied',false,'continuumAccuracyCertified',false, ...
                'normalization','absolute infinity error divided independently by each sampled patch infinity scale; mixed scale uses max(1,core,polar)', ...
                'cfl',struct('definition','max over BOTH active interiors of alpha_x/dx+alpha_y/dy or alpha_r/drMetric+alpha_theta/(r*dtheta), then add abs(c_omega) once; alpha=hypot(velocity,configured physical floor), physical wall normal term zero', ...
                    'activeOverlapIncludedOnBothPatches',true,'onlyAlgebraicGhostCollarsExcluded',true, ...
                    'smallVelocityCutoffApplied',false, ...
                    'requiredCheckScope','cached initial state, every complete explicit RK stage and fresh accepted endpoint', ...
                    'allStageCheckPerformedByThisAudit',false, ...
                    'activeCoreCount',nnz(c.active),'activePolarCount',nnz(p.active), ...
                    'excludedCoreGhostCount',nnz(~c.active),'excludedPolarGhostCount',nnz(~p.active), ...
                    'spatialAccuracyCertifiedByCfl',false));
            if isfield(state,'scale')&&isfield(state.scale,'canonicalTime')
                info.canonicalTime=state.scale.canonicalTime;
            end
            if isfield(state,'step'),info.step=state.step;end
            if isfield(flow,'rate'),info.cfl.reportedEndpointTotalRate=flow.rate;end
            if isfield(flow,'cflPosition'),info.cfl.reportedEndpointPosition=flow.cflPosition;end
            if opts.includeQueryCoordinates
                info.queries=struct('x',cache.x,'y',cache.y,'r',cache.r, ...
                    'theta',cache.theta,'physicalWall',cache.wall);
            end

                function [record,values]=compare(field)
                    assert(isstruct(field)&&isscalar(field)&&all(isfield(field,{'core','polar'}))&& ...
                        isnumeric(field.core)&&isreal(field.core)&&isequal(size(field.core),size(c.X))&& ...
                        isnumeric(field.polar)&&isreal(field.polar)&&isequal(size(field.polar),size(p.X)), ...
                        'ipm_cr:ExplicitAuditFields','Every compared field must use both native patch array shapes.');
                    a=field.core(cr,cc);b=field.polar(pr,pc);
                    if staircase,a=field.core(c.activeIndices);end
                    assert(all(isfinite(a),'all')&&all(isfinite(b),'all'), ...
                        'ipm_cr:ExplicitAuditFields','Active audit donor values must be finite.');
                    av=cache.coreMap*a(:);bv=cache.polarMap*b(:);
                    assert(all(isfinite(av))&&all(isfinite(bv)), ...
                        'ipm_cr:ExplicitAuditInterpolation','Independent audit interpolation must remain finite.');
                    error=abs(av-bv);[maximum,index]=max(error);
                    record=norm_record(maximum,max(abs(av)),max(abs(bv)));
                    record.absoluteRms=sqrt(mean(error.^2));
                    record.maximumLocation=query_location(cache,index);
                    record.coreValueAtMaximum=av(index);record.polarValueAtMaximum=bv(index);
                    record.interiorAbsoluteInfinity=max_or_nan(error(~cache.wall));
                    record.wallAbsoluteInfinity=max_or_nan(error(cache.wall));
                    values=struct('core',av,'polar',bv);
                end
        end
    end
end

function v=rate_values(s)
v=struct('names',{{'c_l','c_omega','c_r','delta','Uref','transportOrigin'}}, ...
    'values',[s.rhs.c_l,s.rhs.c_omega,s.rhs.c_r,s.rhs.delta,s.rhs.Uref,s.rhs.transportOrigin]);
end

function v=scalar_values(s)
tail=zeros(1,2);
if isfield(s.scale,'physicalTimeTail'),tail=reshape(s.scale.physicalTimeTail,1,[]);end
v=struct('names',{{'logC_l','logC_omega','physicalTimeCompensated','X_shift','densityOrigin'}}, ...
    'values',[s.scale.logC_l,s.scale.logC_omega,s.scale.physicalTime+s.scale.physicalTimeLow,s.scale.X_shift,s.scale.densityOrigin], ...
    'clockHigh',s.scale.physicalTime,'clockLow',s.scale.physicalTimeLow, ...
    'clockTail',tail,'clockWords',[s.scale.physicalTime,s.scale.physicalTimeLow,tail]);
end

function opts=resolve_options(options)
assert(isstruct(options)&&isscalar(options),'ipm_cr:ExplicitAuditOptions', ...
    'Use a scalar audit options struct.');
opts=struct('radialSamples',11,'angularSamples',65,'includeWall',true, ...
    'includeQueryCoordinates',false);
names=fieldnames(options);
for j=1:numel(names)
    assert(isfield(opts,names{j}),'ipm_cr:ExplicitAuditOptions', ...
        'Unknown independent audit option %s.',names{j});
    opts.(names{j})=options.(names{j});
end
for name={'radialSamples','angularSamples'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','integer','finite','>=',2});
end
validateattributes(opts.includeWall,{'logical'},{'scalar'});
validateattributes(opts.includeQueryCoordinates,{'logical'},{'scalar'});
end

function cache=prepare_maps(key,core)
cx=key.coreX;cy=key.coreY;rr=key.polarR;tt=key.polarTheta;
assert(all(isfinite([cx;cy;rr;tt]))&&all(diff(cx)>0)&&all(diff(cy)>0)&& ...
    all(diff(rr)>0)&&all(diff(tt)>0)&&cx(1)<0&&cx(end)>0&&cy(1)==0&& ...
    rr(1)>0&&tt(1)==0&&abs(tt(end)-pi)<=16*eps(pi), ...
    'ipm_cr:ExplicitAuditGeometry','Use active Cartesian and origin-free half-annular donor axes.');
theta=pi*((0:key.angularSamples-1)'+.371)/key.angularSamples;
staircase=isfield(key,'stair');
if staircase,stairHi=ipm_cr.mesh.staircase('ray',key.stair,[theta;0;pi]);end
x=[];y=[];r=[];angle=[];wall=false(0,1);
fraction=((0:key.radialSamples-1)'+.413)/key.radialSamples;
for j=1:numel(theta)
    a=theta(j);co=cos(a);si=sin(a);
    if co>0,upperX=cx(end)/co;elseif co<0,upperX=cx(1)/co;else,upperX=Inf;end
    upper=min([rr(end),upperX,cy(end)/si]);lower=rr(1);
    if staircase,upper=min(upper,stairHi(j));end
    if upper<=lower*(1+128*eps),continue;end
    rq=exp(log(lower)+fraction*log(upper/lower));
    x=[x;rq*co];y=[y;rq*si];r=[r;rq];angle=[angle;a*ones(size(rq))]; %#ok<AGROW>
    wall=[wall;false(size(rq))]; %#ok<AGROW>
end
if key.includeWall
    for j=1:2
        a=(j-1)*pi;
        if a==0,upper=min(rr(end),cx(end));else,upper=min(rr(end),-cx(1));end
        lower=rr(1);
        if staircase,upper=min(upper,stairHi(key.angularSamples+j));end
        if upper<=lower*(1+128*eps),continue;end
        rq=exp(log(lower)+fraction*log(upper/lower));
        x=[x;rq*cos(a)];y=[y;zeros(size(rq))];r=[r;rq]; %#ok<AGROW>
        angle=[angle;a*ones(size(rq))];wall=[wall;true(size(rq))]; %#ok<AGROW>
    end
end
assert(~isempty(x)&&all(x>=cx(1)&x<=cx(end)&y>=cy(1)&y<=cy(end))&& ...
    all(r>=rr(1)&r<=rr(end)&angle>=tt(1)&angle<=tt(end)), ...
    'ipm_cr:ExplicitAuditOverlap','The two independent active tensors need a nonempty common query domain.');
if staircase
    coreMap=ipm_cr.mesh.staircase('donors',core,x,y,7);
else
    coreMap=ipm_cr.mesh.interpolate.tensor(cy,cx,y,x,7);
end
polarMap=ipm_cr.mesh.interpolate.tensor(rr,tt,r,angle,7);
oneSidedCounts=struct('coreX',one_sided_count(cx,x), ...
    'coreY',one_sided_count(cy,y),'polarR',one_sided_count(rr,r), ...
    'polarTheta',one_sided_count(tt,angle));
if staircase
    [iq,donor]=find(coreMap);nodes=core.activeIndices(donor);
    [j,i]=ind2sub(size(core.active),nodes);
    [~,cxNearest]=min(abs(x-core.x),[],2);
    [~,cyNearest]=min(abs(y-core.y.'),[],2);
    % An active staircase shifts some stencils away from their centered
    % line even when the query is away from the bounding rectangle edges.
    oneSidedCounts.coreX=numel(unique(iq(abs(i-cxNearest(iq))>3)));
    oneSidedCounts.coreY=numel(unique(iq(abs(j-cyNearest(iq))>3)));
end
cache=struct('key',key,'coreMap',coreMap,'polarMap',polarMap, ...
    'x',x,'y',y,'r',r,'theta',angle,'wall',wall, ...
    'bounds',struct('x',[min(x),max(x)],'y',[min(y),max(y)], ...
        'r',[min(r),max(r)],'theta',[min(angle),max(angle)]), ...
    'donorBounds',struct('coreX',cx([1,end])','coreY',cy([1,end])', ...
        'polarR',rr([1,end])','polarTheta',tt([1,end])'), ...
    'oneSidedCounts',oneSidedCounts);
end

function count=one_sided_count(axis,queries)
nearest=interp1(axis,(1:numel(axis))',queries,'nearest');
count=nnz(nearest<4|nearest>numel(axis)-3);
end

function record=norm_record(error,coreScale,polarScale)
record=struct('absoluteInfinity',error,'coreInfinityScale',coreScale, ...
    'polarInfinityScale',polarScale, ...
    'normalizedByCoreScale',error/max(coreScale,realmin), ...
    'normalizedByPolarScale',error/max(polarScale,realmin), ...
    'mixedRelativeInfinity',error/max([1,coreScale,polarScale]), ...
    'coreScaleIsZero',coreScale==0,'polarScaleIsZero',polarScale==0);
end

function location=query_location(cache,index)
location=struct('queryIndex',index,'x',cache.x(index),'y',cache.y(index), ...
    'r',cache.r(index),'theta',cache.theta(index),'physicalWall',cache.wall(index));
end

function value=max_or_nan(values)
if isempty(values),value=NaN;else,value=max(values);end
end
