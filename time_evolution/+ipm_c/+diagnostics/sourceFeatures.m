function [features,info]=sourceFeatures(q,g,levels)
%SOURCEFEATURES Complete polar-source physical P1 areas and radial bounds.
% Self-contained geometry adapter; no experiments namespace is used at runtime.
% Call with COMPLETE q=D_X(rho), not abs(q), masked q, or one component.
% Each annular quadrilateral uses the same lower-left/upper-right P1
% diagonal as the original ipm area monitor. The origin is one fan vertex;
% no duplicated zero-area origin quadrilaterals are ever constructed.
% All disconnected positive-level components and threshold plateaus count.
%
% FEATURES.featureBounds=[0,R_l] is a padded radial envelope for radialCdf.
% FEATURES.radialProjectionBounds and radialIntervals retain the actual
% clipped radial projection, including gaps. They are not 2-D component IDs.
% Areas are exact for the P1 field on the physical straight-edge polygonal
% mesh. They are not exact curved semicircle integrals or node counts.
% This module measures only: it neither triggers nor performs a remesh.
if nargin<3,levels=[.1,.2,.5];end
assert(isnumeric(q)&&isreal(q)&&ismatrix(q)&&all(isfinite(q),'all'), ...
    'ipm_c:FeatureSource','Supply a finite complete real source field.');
assert(isstruct(g)&&isscalar(g)&&all(isfield(g,{'r','theta','X','Y'})), ...
    'ipm_c:FeatureGrid','Supply one polar geometry with r, theta, X and Y.');
r=g.r(:);theta=g.theta(:)';nr=numel(r);nt=numel(theta);
assert(nr>=2&&nt>=3&&isequal(size(q),[nr,nt])&& ...
    isequal(size(g.X),[nr,nt])&&isequal(size(g.Y),[nr,nt])&& ...
    all(isfinite(r))&&r(1)==0&&all(diff(r)>0)&& ...
    all(isfinite(theta))&&theta(1)==0&&theta(end)==pi&& ...
    all(isfinite(g.X),'all')&&all(isfinite(g.Y),'all'), ...
    'ipm_c:FeatureGrid','The complete source and coordinate arrays must match finite increasing polar axes.');
dt=pi/(nt-1);
assert(max(abs(diff(theta)-dt))<=128*eps(pi), ...
    'ipm_c:FeatureTheta','The candidate retains uniform theta on [0,pi].');
assert(all(g.X(1,:)==0)&&all(g.Y(1,:)==0)&&all(q(1,:)==q(1,1)), ...
    'ipm_c:FeatureOrigin','The origin must be one physical point and one source value.');
coordinateTolerance=128*eps(max(1,r(end)));
assert(max(abs(g.X-r*cos(theta)),[],'all')<=coordinateTolerance&& ...
    max(abs(g.Y-r*sin(theta)),[],'all')<=coordinateTolerance, ...
    'ipm_c:FeatureCoordinates','Coordinates must belong to the declared polar half-disk.');
assert(isnumeric(levels)&&isreal(levels)&&isvector(levels)&&~isempty(levels)&& ...
    all(isfinite(levels))&&all(levels>0&levels<1)&&all(diff(levels)>0), ...
    'ipm_c:FeatureLevels','Use increasing relative positive-peak levels strictly between zero and one.');
levels=reshape(levels,1,[]);nlevels=numel(levels);
[peak,peakIndex]=max(q(:));q=double(q);
features=struct('levels',levels,'areas',zeros(1,nlevels), ...
    'featureBounds',NaN(nlevels,2),'radialProjectionBounds',NaN(nlevels,2), ...
    'radialIntervals',{cell(1,nlevels)},'peak',double(peak),'peakIndex',peakIndex, ...
    'usableForCdf',false,'areaDefinition','physical_triangle_P1_superlevel_clipping', ...
    'featureBoundsDefinition','[0,R_l] padded envelope of every component; not the exact radial union', ...
    'areaTriggerRelative',.2,'areaTriggerImplemented',false);
info=struct('kind','complete_source_physical_triangle_features_v1', ...
    'allDisconnectedComponentsIncluded',true,'nodeCountsUsedAsArea',false, ...
    'twoDimensionalComponentCountComputed',false,'originPhysicalPointCount',1, ...
    'originFanTriangleCount',nt-1,'triangleCount',0,'zeroAreaTrianglesConstructed',0, ...
    'maximumTrianglesPerProcessedBatch',0,'nominalBatchTriangleBudget',32768, ...
    'fullPolygonArea',.5*(nt-1)*sin(dt)*r(end)^2, ...
    'curvedSemicircleArea',.5*pi*r(end)^2, ...
    'polygonRelativeAreaDefect',1-sin(dt)/dt, ...
    'pdeFieldChanged',false,'sourceMaskUsed',false, ...
    'remeshImplemented',false,'spaceOrTrajectoryQualified',false);
if peak<=0,return;end
areas=zeros(1,nlevels);compensation=areas;unions=cell(1,nlevels);

% One triangle per first-ring sector, with a single physical origin.
fanX=[zeros(nt-1,1),g.X(2,1:end-1)',g.X(2,2:end)'];
fanY=[zeros(nt-1,1),g.Y(2,1:end-1)',g.Y(2,2:end)'];
fanQ=[q(1,1)*ones(nt-1,1),q(2,1:end-1)',q(2,2:end)'];
consume(fanX,fanY,fanQ);
% Bound workspace by radial slabs; there is no amplitude-based field crop.
rowsPerBatch=max(1,floor(32768/(2*(nt-1))));
for first=2:rowsPerBatch:nr-1
    last=min(nr-1,first+rowsPerBatch-1);rows=first:last;
    [x,y,values]=annular_triangles(rows);
    consume(x,y,values);
end
features.areas=areas+compensation;
for level=1:nlevels
    interval=unions{level};
    assert(~isempty(interval)&&all(isfinite(interval),'all'), ...
        'ipm_c:FeatureProjection','Every positive relative level must have a finite radial projection.');
    features.radialIntervals{level}=interval;
    features.radialProjectionBounds(level,:)=[interval(1,1),interval(end,2)];
    features.featureBounds(level,:)=[0,interval(end,2)];
end
features.usableForCdf=all(isfinite(features.areas))&&all(features.areas>0)&& ...
    all(features.featureBounds(:,2)>0);
assert(info.triangleCount==(2*nr-3)*(nt-1)&& ...
    all(isfinite(features.areas))&&all(features.areas>=0), ...
    'ipm_c:FeatureMeasurement','The complete physical triangulation and finite nonnegative areas must be retained.');

    function [x,y,values]=annular_triangles(rows)
        x00=g.X(rows,1:end-1);x10=g.X(rows+1,1:end-1);
        x11=g.X(rows+1,2:end);x01=g.X(rows,2:end);
        y00=g.Y(rows,1:end-1);y10=g.Y(rows+1,1:end-1);
        y11=g.Y(rows+1,2:end);y01=g.Y(rows,2:end);
        q00=q(rows,1:end-1);q10=q(rows+1,1:end-1);
        q11=q(rows+1,2:end);q01=q(rows,2:end);
        x=[x00(:),x10(:),x11(:);x00(:),x11(:),x01(:)];
        y=[y00(:),y10(:),y11(:);y00(:),y11(:),y01(:)];
        values=[q00(:),q10(:),q11(:);q00(:),q11(:),q01(:)];
    end
    function consume(x,y,values)
        triangleArea=.5*abs((x(:,2)-x(:,1)).*(y(:,3)-y(:,1))- ...
            (x(:,3)-x(:,1)).*(y(:,2)-y(:,1)));
        assert(all(isfinite(triangleArea))&&all(triangleArea>0), ...
            'ipm_c:FeatureTriangle','Every constructed physical triangle must have positive finite area.');
        info.triangleCount=info.triangleCount+numel(triangleArea);
        info.maximumTrianglesPerProcessedBatch=max(info.maximumTrianglesPerProcessedBatch,numel(triangleArea));
        sorted=sort(values,2);
        for level=1:nlevels
            threshold=levels(level)*double(peak);
            fraction=triangle_fraction(sorted,threshold);
            subtotal=sum(triangleArea.*fraction);
            sumHigh=areas(level)+subtotal;tail=sumHigh-areas(level);
            compensation(level)=compensation(level)+ ...
                (areas(level)-(sumHigh-tail))+(subtotal-tail);
            areas(level)=sumHigh;
            exists=any(values>=threshold,2);
            if any(exists)
                bounds=clipped_radial_ranges(x(exists,:),y(exists,:),values(exists,:),threshold);
                unions{level}=merge_intervals([unions{level};bounds]);
            end
        end
    end
end

function fraction=triangle_fraction(sorted,threshold)
% Original ipm closed-form area fraction, including exact-threshold plateaus.
a=sorted(:,1);b=sorted(:,2);c=sorted(:,3);fraction=zeros(size(a));
full=threshold<=a;fraction(full)=1;partial=~full&threshold<c;
lower=partial&threshold<=b;upper=partial&threshold>b;
fraction(lower)=1-((threshold-a(lower))./(b(lower)-a(lower))).* ...
    ((threshold-a(lower))./(c(lower)-a(lower)));
fraction(upper)=((c(upper)-threshold)./(c(upper)-a(upper))).* ...
    ((c(upper)-threshold)./(c(upper)-b(upper)));
fraction=min(max(fraction,0),1);
end

function bounds=clipped_radial_ranges(x,y,q,threshold)
% The maximum radius is at a clipped polygon vertex. Its minimum can lie
% inside the new clipping edge, so inspecting threshold vertices alone is
% insufficient. Include nearest points on every retained/new polygon edge.
n=size(x,1);minimum=Inf(n,1);maximum=zeros(n,1);
firstX=zeros(n,1);firstY=firstX;secondX=firstX;secondY=firstX;
hasFirst=false(n,1);hasSecond=hasFirst;
edges=[1,2;2,3;3,1];
for edge=1:3
    a=edges(edge,1);b=edges(edge,2);insideA=q(:,a)>=threshold;insideB=q(:,b)>=threshold;
    keep=insideA|insideB;cross=insideA~=insideB;
    xa=x(:,a);ya=y(:,a);xb=x(:,b);yb=y(:,b);
    if any(cross)
        t=(threshold-q(cross,a))./(q(cross,b)-q(cross,a));t=min(max(t,0),1);
        xc=xa(cross)+t.*(xb(cross)-xa(cross));
        yc=ya(cross)+t.*(yb(cross)-ya(cross));
        positions=find(cross);second=hasFirst(positions);first=~second;
        firstX(positions(first))=xc(first);firstY(positions(first))=yc(first);
        secondX(positions(second))=xc(second);secondY(positions(second))=yc(second);
        hasSecond(positions(second))=true;hasFirst(positions(first))=true;
        replaceA=~insideA(positions);replaceB=~insideB(positions);
        xa(positions(replaceA))=xc(replaceA);ya(positions(replaceA))=yc(replaceA);
        xb(positions(replaceB))=xc(replaceB);yb(positions(replaceB))=yc(replaceB);
    end
    minimum(keep)=min(minimum(keep),segment_radius(xa(keep),ya(keep),xb(keep),yb(keep)));
    maximum(keep)=max(maximum(keep),max(hypot(xa(keep),ya(keep)),hypot(xb(keep),yb(keep))));
end
assert(all(hasFirst==hasSecond),'ipm_c:FeatureClipEdge', ...
    'A partially clipped triangle must have both endpoints of its new clipping edge.');
if any(hasSecond)
    minimum(hasSecond)=min(minimum(hasSecond),segment_radius( ...
        firstX(hasSecond),firstY(hasSecond),secondX(hasSecond),secondY(hasSecond)));
end
assert(all(isfinite(minimum))&&all(isfinite(maximum))&&all(minimum<=maximum), ...
    'ipm_c:FeatureClipRange','Clipped radial ranges must be finite and ordered.');
bounds=[minimum,maximum];
end

function radius=segment_radius(xa,ya,xb,yb)
dx=xb-xa;dy=yb-ya;denominator=dx.^2+dy.^2;t=zeros(size(dx));
nonzero=denominator>0;
t(nonzero)=-(xa(nonzero).*dx(nonzero)+ya(nonzero).*dy(nonzero))./denominator(nonzero);
t=min(max(t,0),1);radius=hypot(xa+t.*dx,ya+t.*dy);
end

function merged=merge_intervals(intervals)
if isempty(intervals),merged=zeros(0,2);return;end
intervals=sortrows(intervals,[1,2]);right=cummax(intervals(:,2));
first=[1;find(intervals(2:end,1)>right(1:end-1))+1];
last=[first(2:end)-1;size(intervals,1)];
merged=[intervals(first,1),right(last)];
end
