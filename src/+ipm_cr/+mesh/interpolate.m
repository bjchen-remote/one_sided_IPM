classdef interpolate
    %INTERPOLATE Local interpolation weights and sampled field evaluation.
    methods (Static)
        function M=local(nodes,queries,width,extrapolate)
            %LOCALMAP Sparse local polynomial value map on arbitrary increasing nodes.
            % Barycentric weights use the actual local nodes after translation/scaling;
            % Chebyshev endpoint weights are inapplicable to the stretched carrier axis.
            % EXTRAPOLATE explicitly continues the nearest endpoint polynomial. The
            % default forbids extrapolation; remesh uses this continuation to retain
            % each independently evolved patch when its active boundary moves.
            if nargin<4,extrapolate=false;end
            validateattributes(extrapolate,{'logical'},{'scalar'});
            nodes=nodes(:);queries=queries(:);n=numel(nodes);nq=numel(queries);
            assert(all(isfinite(nodes))&&all(diff(nodes)>0)&&all(isfinite(queries))&& ...
                width<=n,'ipm_cr:CarrierInterpolation','Finite increasing interpolation axes are required.');
            lowerTolerance=64*eps(max(1,abs(nodes(1))));
            upperTolerance=64*eps(max(1,abs(nodes(end))));
            if ~extrapolate
                assert(all(queries>=nodes(1)-lowerTolerance&queries<=nodes(end)+upperTolerance), ...
                    'ipm_cr:CarrierInterpolation','Local source interpolation cannot extrapolate.');
                queries=max(nodes(1),min(nodes(end),queries));
            end
            % Clip only the stencil selector when extending; evaluate at the actual
            % query, rather than flattening the continued field at the old boundary.
            selector=max(nodes(1),min(nodes(end),queries));
            if mod(width,2)==0
                % Even stencils center on the bracketing interval. The stencil then
                % changes at a donor node where both polynomials agree.
                left=discretize(selector,nodes);
                first=left-(width/2-1);
            else
                nearest=interp1(nodes,(1:n)',selector,'nearest');
                first=nearest-floor(width/2);
            end
            first=max(1,min(n-width+1,first));
            columns=first+(0:width-1);weights=zeros(nq,width);
            for k=1:nq
                x=nodes(columns(k,:));q=queries(k);
                % Use each node's ulp, rather than the outer-domain scale, to identify
                % exact local hits near the origin of a disk with radius 1e10.
                hit=find(abs(q-x)<=8*eps(max(1,abs(x))),1);
                if ~isempty(hit),weights(k,hit)=1;continue;end
                span=x(end)-x(1);z=(x-x(1))/span;zq=(q-x(1))/span;
                bary=ones(width,1);
                for j=1:width
                    other=[1:j-1,j+1:width];bary(j)=1/prod(z(j)-z(other));
                end
                bary=bary/max(abs(bary));v=bary./(zq-z);v=v/sum(v);
                assert(all(isfinite(v)),'ipm_cr:CarrierInterpolation','Nonfinite local barycentric weights.');
                weights(k,:)=v';
            end
            rows=repmat((1:nq)',1,width);
            M=sparse(rows(:),columns(:),weights(:),nq,n);
        end

        function M=tensor(r,theta,rq,tq,width,extrapolate)
            %TENSORMAP Sparse rowwise local interpolation of an r-by-theta field.
            % Optional endpoint polynomial continuation is explicit; sampling and source
            % maps keep the default no-extrapolation contract.
            if nargin<6,extrapolate=false;end
            assert(numel(rq)==numel(tq),'ipm_cr:CarrierInterpolation','Query axes must have matching lengths.');
            R=ipm_cr.mesh.interpolate.local(r,rq,width,extrapolate);
            T=ipm_cr.mesh.interpolate.local(theta,tq,width,extrapolate);
            nq=numel(rq);nr=numel(r);nt=numel(theta);
            [ri,rc,rv]=find(R');[ti,tc,tv]=find(T'); %#ok<ASGLU>
            % Each query has at most width nonzeros. Build compact rowwise lists after
            % dropping exact-node zero weights, avoiding a global tensor Kronecker map.
            rcols=zeros(nq,width);rvals=rcols;tcols=rcols;tvals=rcols;
            counts=zeros(nq,1);
            for k=1:numel(rc),j=rc(k);counts(j)=counts(j)+1;rcols(j,counts(j))=ri(k);rvals(j,counts(j))=rv(k);end
            counts=zeros(nq,1);
            for k=1:numel(tc),j=tc(k);counts(j)=counts(j)+1;tcols(j,counts(j))=ti(k);tvals(j,counts(j))=tv(k);end
            rows=cell(width,1);cols=rows;vals=rows;
            for j=1:width
                valid=tcols(:,j)>0;
                row=repmat(find(valid),1,width);
                col=rcols(valid,:)+nr*(tcols(valid,j)-1);
                value=rvals(valid,:).*tvals(valid,j);
                use=rcols(valid,:)>0&value~=0;
                % A single surviving query makes MATLAB preserve row orientation.
                % Normalize every block before concatenating unequal stencil counts.
                rows{j}=reshape(row(use),[],1);cols{j}=reshape(col(use),[],1);
                vals{j}=reshape(value(use),[],1);
            end
            M=sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),nq,nr*nt);
        end

        function [values,map]=sample(fields,g,x,y)
            %INTERPOLATE Local degree-six values with one active-core ownership rule.
            % FIELDS.core is y-by-x; FIELDS.polar is r-by-theta, without a polar origin.
            % Inside the active Cartesian cells use Cartesian data; elsewhere use
            % the annulus. Full arrays, including synchronized halos, supply stencils.
            % Optional MAP contains geometry-only sparse core/polar maps for reuse:
            % values=ipm_cr.mesh.interpolate.sample(fields,map) uses the identical anchored apply.
            % Invalid queries return NaN. Only tiny wall/outer-radius roundoff is repaired.
            if nargin==2
                map=g;values=apply_map(fields,map);return
            end
            assert(isstruct(g)&&all(isfield(g,{'core','polar'})), ...
                'ipm_cr:ExplicitInterpolationGrid','Supply the hybrid transport geometry.');
            c=g.core;p=g.polar;xAxis=c.x(:)';yAxis=c.y(:);r=p.r(:);theta=p.theta(:)';
            assert(isstruct(fields)&&isscalar(fields)&&all(isfield(fields,{'core','polar'}))&& ...
                isnumeric(fields.core)&&isreal(fields.core)&&isequal(size(fields.core),[numel(yAxis),numel(xAxis)])&& ...
                isnumeric(fields.polar)&&isreal(fields.polar)&&isequal(size(fields.polar),[numel(r),numel(theta)]), ...
                'ipm_cr:ExplicitInterpolationField','Use complete native core/polar scalar fields.');
            assert(isnumeric(x)&&isreal(x)&&isnumeric(y)&&isreal(y)&&isequal(size(x),size(y)), ...
                'ipm_cr:ExplicitInterpolationQueries','Use matching real query arrays.');
            assert(all(isfield(c,{'activeRows','activeColumns'})), ...
                'ipm_cr:ExplicitInterpolationGrid','The geometry must identify the active Cartesian cells.');
            columns=c.activeColumns;rows=c.activeRows;
            bounds=[xAxis(columns(1)),xAxis(columns(end)),yAxis(rows(1)),yAxis(rows(end))];
            key=struct('x',xAxis,'y',yAxis,'r',r,'theta',theta,'bounds',bounds);
            if isfield(c,'stair')
                key.stair=struct('edges',c.stair.edges,'top',c.stair.top);
            end
            persistent retained
            hit=[];
            for k=1:numel(retained)
                if isequaln(retained{k}.key,key)&&isequaln(retained{k}.x,x)&&isequaln(retained{k}.y,y)
                    hit=k;break;
                end
            end
            if isempty(hit)
                map=prepare_map(key,x,y);
                entry=struct('key',key,'x',x,'y',y,'map',map);
                if isempty(retained),retained={entry};else,retained{end+1}=entry;end
                % Retain only geometry/query maps, never a sampled field or value.
                if numel(retained)>4,retained=retained(end-3:end);end
            else
                map=retained{hit}.map;
            end
            values=apply_map(fields,map);
        end

        function D = metricDerivative(axis,derivativeOrder,stencilWidth)
            %IPM.MESH.FDMATRIX Local-polynomial finite-difference matrix on any monotone axis.
            %   D = IPM.MESH.FDMATRIX(X,M,WIDTH) differentiates nodal data M times.  Every
            %   row is exact for polynomials through degree WIDTH-1, including the
            %   one-sided boundary closures.  WIDTH=7 therefore gives sixth-order first
            %   derivatives and at least fifth-order second derivatives on a smooth grid.

            if nargin < 3
                stencilWidth = 7;
            end
            validateattributes(axis,{'numeric'}, ...
                {'vector','real','finite','increasing'},mfilename,'axis');
            validateattributes(derivativeOrder,{'numeric'}, ...
                {'scalar','integer','nonnegative'},mfilename,'derivativeOrder');
            validateattributes(stencilWidth,{'numeric'}, ...
                {'scalar','integer','>=',derivativeOrder+1},mfilename,'stencilWidth');

            axis = axis(:);
            numberOfNodes = numel(axis);
            if stencilWidth > numberOfNodes
                error('ipm:FiniteDifferenceStencil', ...
                    'The stencil width cannot exceed the number of axis nodes.');
            end

            rows = zeros(numberOfNodes*stencilWidth,1);
            columns = zeros(numberOfNodes*stencilWidth,1);
            values = zeros(numberOfNodes*stencilWidth,1);
            halfWidth = floor(stencilWidth/2);
            cursor = 0;
            for node = 1:numberOfNodes
                first = min(max(node-halfWidth,1),numberOfNodes-stencilWidth+1);
                stencil = first:first+stencilWidth-1;
                offsets = axis(stencil)-axis(node);
                scale = max(abs(offsets));
                if scale == 0
                    error('ipm:FiniteDifferenceScale','A finite-difference stencil collapsed.');
                end
                normalized = offsets/scale;
                moments = normalized'.^((0:stencilWidth-1)');
                target = zeros(stencilWidth,1);
                target(derivativeOrder+1) = ...
                    factorial(derivativeOrder)/scale^derivativeOrder;
                coefficients = moments\target;
                slots = cursor+(1:stencilWidth);
                rows(slots) = node;
                columns(slots) = stencil;
                values(slots) = coefficients;
                cursor = cursor+stencilWidth;
            end
            D = sparse(rows,columns,values,numberOfNodes,numberOfNodes);
        end
    end
end

function value=apply_map(fields,map)
% Each query is anchored at its strongest LOCAL donor. A remote nonzero
% platform must not inject a row-sum defect into a stencil of exact zeros.
native=~isfield(map,'referenceOnly')||~map.referenceOnly;
value=local_apply(map.core,fields.core(:),map.coreAnchor,native)+ ...
    local_apply(map.polar,fields.polar(:),map.polarAnchor,native);
value(~map.valid)=NaN;value=reshape(value,map.shape);
end

function value=local_apply(M,donor,anchor,native)
binary=which('ipm_cr.mesh.sampleMex');
if native&&~isempty(binary)&&exist(binary,'file')==3
    value=ipm_cr.mesh.sampleMex(M,donor,anchor);return
end
% Bounded-memory reference path, identical to the native per-row formula.
value=zeros(size(M,1),1);
for first=1:2048:numel(value)
    rows=first:min(first+2047,numel(value));index=anchor(rows);
    reference=zeros(numel(rows),1);use=index>0;reference(use)=donor(index(use));
    [i,j,w]=find(M(rows,:));
    i=i(:);j=j(:);w=w(:);
    delta=w.*(donor(j)-reshape(reference(i),[],1));
    value(rows)=reference+accumarray(i,delta,[numel(rows),1],@sum,0);
end
end

function map=prepare_map(key,x,y)
shape=size(x);x=x(:);y=y(:);n=numel(x);width=7;
assert(min([numel(key.x),numel(key.y),numel(key.r),numel(key.theta)])>=width&& ...
    all(isfinite([key.x(:);key.y(:);key.r(:);key.theta(:)]))&& ...
    all(diff(key.x)>0)&&all(diff(key.y)>0)&&all(diff(key.r)>0)&& ...
    all(diff(key.theta)>0)&&key.r(1)>0, ...
    'ipm_cr:ExplicitInterpolationAxes','Finite increasing axes with seven nodes and positive annular radii are required.');
wallTolerance=8*eps(max(1,abs(x)));
snap=isfinite(x)&isfinite(y)&y<0&y>=-wallTolerance;y(snap)=0;
radius=hypot(x,y);finite=isfinite(x)&isfinite(y)&y>=0;
withinDisk=radius<=key.r(end)*(1+16*eps);
coreOwner=finite&withinDisk&x>=key.bounds(1)&x<=key.bounds(2)& ...
    y>=key.bounds(3)&y<=key.bounds(4);
if isfield(key,'stair')
    coreOwner=finite&withinDisk&ipm_cr.mesh.staircase('contains',key.stair,x,y);
end
polarOwner=finite&withinDisk&~coreOwner&radius>=key.r(1)*(1-16*eps);
core=spalloc(n,numel(key.x)*numel(key.y),nnz(coreOwner)*width^2);
polar=spalloc(n,numel(key.r)*numel(key.theta),nnz(polarOwner)*width^2);
indices=find(coreOwner);
if ~isempty(indices)
    local=ipm_cr.mesh.interpolate.tensor(key.y,key.x,y(indices),x(indices),width);
    [row,column,value]=find(local);core=sparse(indices(row),column,value,n,size(core,2));
end
indices=find(polarOwner);
if ~isempty(indices)
    rq=max(key.r(1),min(key.r(end),radius(indices)));
    tq=max(key.theta(1),min(key.theta(end),atan2(y(indices),x(indices))));
    local=ipm_cr.mesh.interpolate.tensor(key.r,key.theta,rq,tq,width);
    [row,column,value]=find(local);polar=sparse(indices(row),column,value,n,size(polar,2));
end
valid=coreOwner|polarOwner;
rowNorm=full(sum(abs(core),2)+sum(abs(polar),2));
if isempty(rowNorm),operatorNorm=0;else,operatorNorm=max(rowNorm);end
[~,coreAnchor]=max(abs(core),[],2);[~,polarAnchor]=max(abs(polar),[],2);
coreAnchor=full(coreAnchor);polarAnchor=full(polarAnchor);
coreAnchor(~coreOwner)=0;polarAnchor(~polarOwner)=0;
map=struct('core',core,'polar',polar,'shape',shape,'valid',valid, ...
    'coreAnchor',coreAnchor,'polarAnchor',polarAnchor, ...
    'coreOwner',coreOwner,'polarOwner',polarOwner,'allValid',all(valid), ...
    'activeCoreBounds',key.bounds,'interpolationWidth',width, ...
    'operatorInfinityNorm',operatorNorm,'fieldCacheRetained',false);
end
