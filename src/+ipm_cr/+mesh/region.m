classdef region
    %REGION Wall-connected superlevels, sections and remesh settings.
    methods (Static)
        function shapes=measure(s,g,p)
            % Piecewise-linear superlevel geometry of the original active-owner field.
            % Both levels share one triangulation covering the lowest superlevel set.
            assert(isstruct(s)&&isfield(s,'rhs')&&isfield(s.rhs,'q')&& ...
                isequal(size(s.rhs.q.core),size(g.core.X))&& ...
                isequal(size(s.rhs.q.polar),size(g.polar.X)), ...
                'ipm_cr:RemeshRegion','Use native q arrays matching both patches.');
            q=s.rhs.q;ca=g.core.active;pa=g.polar.active;
            cx=g.core.X(ca);cy=g.core.Y(ca);px=g.polar.X(pa);py=g.polar.Y(pa);
            qc=q.core(ca);qp=q.polar(pa);[pc,ic]=max(qc);[pp,ip]=max(qp);
            shape=struct('valid',false,'reason','nonpositive_peak','level',NaN, ...
                'definition','positive_peak_connected_component_projection_v1', ...
                'method','active_owner_local_linear_triangles_subcell_crossings', ...
                'peak',max(pc,pp),'peakXY',[NaN,NaN],'ownedPeak',NaN,'ownedPeakXY',[NaN,NaN], ...
                'peakDifference',NaN,'peakRelativeDifference',NaN,'threshold',NaN, ...
                'xBounds',[NaN,NaN],'yBounds',[NaN,NaN],'widths',[NaN,NaN],'extents',[NaN,NaN], ...
                'triggerExtents',[NaN,NaN], ...
                'monitorRadius',NaN,'monitorNodes',0,'triangleCount',0,'componentNodes',0, ...
                'otherPositiveNodes',0,'crossings',zeros(0,2),'wallConnected',false, ...
                'globalActiveCount',numel(qc)+numel(qp),'ghostDonorsUsed',false, ...
                'overlapOwner','core_inside_closed_active_rectangle');
            if pc>=pp,shape.peakXY=[cx(ic),cy(ic)];else,shape.peakXY=[px(ip),py(ip)];end
            shapes=fill_levels(shape,p.level);
            if ~(shape.peak>0),return;end
            threshold=p.level(1)*shape.peak;
            radius=max(hypot(cx,cy));above=qp>=threshold;
            if any(above),radius=max(radius,max(hypot(px(above),py(above))));end
            rr=g.polar.r(g.polar.activeRows);
            outer=find(rr>radius+64*eps(max(1,radius)),1);
            if isempty(outer)
                shape.reason='level_set_not_enclosed_by_active_domain';
                shapes=fill_levels(shape,p.level);return
            end
            shape.monitorRadius=rr(outer);
            inside=px>=min(cx)&px<=max(cx)&py>=min(cy)&py<=max(cy);
            if isfield(g.core,'stair')
                inside=ipm_cr.mesh.staircase('contains',g.core.stair,px,py);
                shape.overlapOwner='core_inside_closed_active_staircase';
            end
            take=~inside&hypot(px,py)<=shape.monitorRadius+64*eps(max(1,shape.monitorRadius));
            points=[cx,cy;px(take),py(take)];values=[qc;qp(take)];
            [points,keep]=unique(points,'rows','stable');values=values(keep);
            shape.monitorNodes=size(points,1);[shape.ownedPeak,seed]=max(values);
            shape.ownedPeakXY=points(seed,:);shape.peakDifference=shape.peak-shape.ownedPeak;
            shape.peakRelativeDifference=shape.peakDifference/shape.peak;
            triangles=delaunay(points(:,1),points(:,2));shape.triangleCount=size(triangles,1);
            [edges,~,edgeId]=unique(sort([triangles(:,[1,2]);triangles(:,[2,3]);triangles(:,[3,1])],2),'rows');
            triangleEdges=reshape(edgeId,[],3);
            shapes=fill_levels(shape,p.level);
            for k=1:numel(shapes)
                shapes(k)=measure_level(shapes(k),points,values,seed,edges,triangleEdges);
            end
        end

        function p=policy(overrides)
            % Resolve and validate the shape-based remesh settings.
            assert(isstruct(overrides)&&isscalar(overrides),'ipm_cr:RemeshOptions','Use scalar remesh settings.');
            p=struct('enabled',true,'level',[.25,.80],'threshold',.25,'checkEveryTau',.5, ...
                'minCheckEveryTau',.005,'candidateAxes',[],'stopTau',[]);
            names=fieldnames(overrides);
            assert(all(ismember(names,fieldnames(p))),'ipm_cr:RemeshOptions','Unknown remesh option.');
            for k=1:numel(names),p.(names{k})=overrides.(names{k});end
            validateattributes(p.enabled,{'logical'},{'scalar'});
            skip={'enabled','candidateAxes','level','stopTau'};
            allNames=fieldnames(p);
            for k=1:numel(allNames)
                key=allNames{k};if ismember(key,skip),continue;end
                validateattributes(p.(key),{'numeric'},{'scalar','real','finite','positive'});
            end
            validateattributes(p.level,{'numeric'},{'vector','real','finite','>',0,'<',1,'increasing'});
            assert(numel(p.level)<=2,'ipm_cr:RemeshOptions','Use one or two increasing superlevel fractions.');
            p.level=reshape(p.level,1,[]);
            if ~isempty(p.stopTau)
                validateattributes(p.stopTau,{'numeric'},{'scalar','real','finite','nonnegative'});
            end
            assert(p.minCheckEveryTau<=p.checkEveryTau,'ipm_cr:RemeshOptions', ...
                'minCheckEveryTau must not exceed checkEveryTau.');
            assert(isempty(p.candidateAxes)||(isstruct(p.candidateAxes)&&isscalar(p.candidateAxes)), ...
                'ipm_cr:RemeshOptions','candidateAxes must be empty or an explicit-axis structure.');
        end
    end
end

function shapes=fill_levels(shape,levels)
shapes=repmat(shape,size(levels));
for k=1:numel(shapes)
    shapes(k).level=levels(k);
    shapes(k).threshold=levels(k)*shape.peak;
end
end

function shape=measure_level(shape,points,values,seed,edges,triangleEdges)
threshold=shape.threshold;
if shape.ownedPeak<threshold,shape.reason='peak_absent_from_owner_field';return;end
positive=values>=threshold;joined=edges(positive(edges(:,1))&positive(edges(:,2)),:);
n=size(points,1);adjacency=sparse([joined(:,1);joined(:,2)],[joined(:,2);joined(:,1)],1,n,n);
component=false(n,1);component(seed)=true;queue=zeros(n,1);queue(1)=seed;head=1;tail=1;
while head<=tail
    neighbors=find(adjacency(:,queue(head)));head=head+1;
    fresh=neighbors(~component(neighbors));component(fresh)=true;
    queue(tail+(1:numel(fresh)))=fresh;tail=tail+numel(fresh);
end
shape.componentNodes=nnz(component);shape.otherPositiveNodes=nnz(positive&~component);
shape.wallConnected=any(component&points(:,2)==0);
if ~shape.wallConnected,shape.reason='main_component_detached_from_wall';return;end
isCrossed=xor(component(edges(:,1)),component(edges(:,2)));
crossed=edges(isCrossed,:);
a=crossed(:,1);b=crossed(:,2);
fraction=(threshold-values(a))./(values(b)-values(a));
crossings=nan(size(edges,1),2);
crossings(isCrossed,:)=points(a,:)+fraction.*(points(b,:)-points(a,:));
shape.crossings=unique(crossings(isCrossed,:),'rows','stable');
support=[shape.crossings;points(component,:)];
shape.xBounds=[min(support(:,1)),max(support(:,1))];
shape.yBounds=[min(support(:,2)),max(support(:,2))];
shape.widths=[diff(shape.xBounds),diff(shape.yBounds)];
shape.extents=[max(abs(shape.xBounds)),shape.yBounds(2)];
height=peak_section(crossings,triangleEdges,isCrossed,shape.ownedPeakXY);
if isempty(height),shape.reason='peak_normal_section_not_enclosed';return;end
shape.triggerExtents=[shape.extents(1),height];
shape.valid=true;shape.reason='ok';
end

function height=peak_section(crossings,triangleEdges,isCrossed,peak)
% Every mixed triangle cuts the contour in a segment. Intersect that
% segment with the wall-normal ray and take the first exit above the peak.
cut=triangleEdges(any(isCrossed(triangleEdges),2),:);
cut(~isCrossed(cut))=0;
cut=sort(cut,2,'descend');
a=crossings(cut(:,1),:);b=crossings(cut(:,2),:);
x=peak(1);tol=16*eps(max(1,abs(x)));
between=x>=min(a(:,1),b(:,1))-tol&x<=max(a(:,1),b(:,1))+tol;
dx=b(:,1)-a(:,1);slanted=between&dx~=0;vertical=between&dx==0;
t=max(0,min(1,(x-a(slanted,1))./dx(slanted)));
y=[a(slanted,2)+t.*(b(slanted,2)-a(slanted,2));a(vertical,2);b(vertical,2)];
height=min(y(y>peak(2)));
end
