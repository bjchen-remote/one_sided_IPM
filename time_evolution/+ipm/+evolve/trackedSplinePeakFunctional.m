function peak=trackedSplinePeakFunctional(values,ops,initialSeedX)
%TRACKEDSPLINEPEAKFUNCTIONAL A continued local peak of one fixed C2 spline.
% This is a different observable from the three-node quadratic/global-node
% peak. Natural cubic boundary conditions are imposed at fixed physical
% coordinates pinX +/- peakTrackingHalfWidth, independent of remeshing.
% Endpoint samples are linear degree-five interpolants of the wall samples.
% Pdot=w*qdot and adot=(wx*qdot)/K use exact linear spline functionals.
% Only a nondegenerate stationary branch overlapping its previous basin is
% continued. A missing/ambiguous branch is rejected, never replaced by the
% tallest competing node or stationary point.
x=ops.x(:)';q=values(:)';r=ops.rescaling;
if nargin<3,initialSeedX=r.pinX;end
definition='tracked_natural_cubic_spline';
assert(numel(q)==numel(x) && numel(x)>=9 && all(isfinite(q)) && ...
    all(isfinite(x)) && all(diff(x)>0), ...
    'ipm:DegenerateNormalization','Invalid tracked-spline wall samples.');
required={'pinX','peakTrackingHalfWidth','peakPositionRestoringRate', ...
    'peakEvaluationTime','peakTrackingTimeScale','cOmegaGaugeFloor'};
assert(all(isfield(r,required)), ...
    'ipm:DegenerateNormalization','Tracked spline controls are incomplete.');
window=r.pinX+[-1,1]*r.peakTrackingHalfWidth;
assert(all(isfinite(window)) && window(1)>x(1) && window(2)<x(end) && ...
    window(2)>window(1) && isfinite(r.peakEvaluationTime) && ...
    r.peakEvaluationTime>=0 && isfinite(r.peakTrackingTimeScale) && ...
    r.peakTrackingTimeScale>0, ...
    'ipm:DegenerateNormalization','Invalid fixed spline window or time.');
g=spline_geometry(x,window);
y=g.E*q';
second=zeros(numel(g.knots),1);
second(2:end-1)=g.L'\(g.L\(g.B*y));
critical=zeros(0,5); % x, cell, t, value, K=-Qxx
for cellIndex=1:numel(g.h)
    H=g.h(cellIndex);m0=second(cellIndex);m1=second(cellIndex+1);
    c0=y(cellIndex+1)-y(cellIndex)+H^2*(-2*m0-m1)/6;
    c1=H^2*m0;c2=H^2*(m1-m0)/2;
    roots=real_quadratic_roots(c2,c1,c0);
    for j=1:numel(roots)
        t=roots(j);tol=max(128*eps(max(1,abs(t))), ...
            128*eps(max([1,abs(g.knots(cellIndex)), ...
            abs(g.knots(cellIndex+1))]))/H);
        if t < -tol || t > 1+tol,continue;end
        t=min(max(t,0),1);
        a=g.knots(cellIndex)+H*t;
        if a<=window(1) || a>=window(2),continue;end
        P=(1-t)*y(cellIndex)+t*y(cellIndex+1)+H^2/6* ...
            (((1-t)^3-(1-t))*m0+(t^3-t)*m1);
        K=-((1-t)*m0+t*m1);
        critical(end+1,:)=[a,cellIndex,t,P,K]; %#ok<AGROW>
    end
end
if isempty(critical)
    error('ipm:TrackedPeakBranchLost','The fixed spline has no stationary peak.');
end
critical=sortrows(critical,1);
keep=true(size(critical,1),1);
for j=2:size(critical,1)
    if abs(critical(j,1)-critical(j-1,1)) <= ...
            256*eps(max([1,abs(critical(j,1)),abs(critical(j-1,1))]))
        keep(j)=false; % Same C2 stationary point at a shared cell endpoint.
    end
end
critical=critical(keep,:);
floors=512*eps(1)*max(max(abs(y)),realmin)./g.h(critical(:,2)).^2;
valid=critical(:,4)>r.cOmegaGaugeFloor & critical(:,5)>floors;
candidates=critical(valid,:);
if isempty(candidates)
    error('ipm:TrackedPeakBranchLost', ...
        'The continued spline has no resolved positive concave stationary peak.');
end
allMaxima=candidates;
anchor=[];
if isfield(r,'trackedPeakAnchor'),anchor=r.trackedPeakAnchor;end
if isempty(anchor)
    assert(isnumeric(initialSeedX) && isscalar(initialSeedX) && ...
        isfinite(initialSeedX) && initialSeedX>window(1) && initialSeedX<window(2), ...
        'ipm:TrackedPeakAnchor','Invalid initial peak-branch seed.');
    expected=initialSeedX;allowed=Inf;
else
    ipm.evolve.validateTrackedPeakAnchor(anchor,x,window,[],true);
    deltaTau=r.peakEvaluationTime-anchor.tau;
    expected=anchor.x*exp(-r.peakPositionRestoringRate*deltaTau);
    if ~isfinite(expected) || expected<=window(1) || expected>=window(2)
        error('ipm:TrackedPeakBranchLost','The predicted peak left its fixed window.');
    end
    currentH=g.h(max(1,min(numel(g.h), ...
        find(g.knots<=expected,1,'last'))));
    cellScale=max(anchor.cellWidth,currentH);
    rounding=1024*eps(max([1,abs(anchor.x),abs(expected)]));
    if isequal(anchor.gridX,x)
        % The phase ODE predicts a exactly. Permit a conservative fraction
        % of a cell on normal stages; the allowance shrinks with stage time,
        % so an O(cell) jump cannot survive repeated dt halving.
        allowed=max(rounding,min(.25*cellScale, ...
            4*cellScale*abs(deltaTau)/r.peakTrackingTimeScale));
    else
        % Explicit grid transfer: interpolation may move a at the same tau.
        % Keep branch identity/basin and a quarter-cell spatial guard.
        allowed=max(rounding,.25*cellScale);
    end
    inOldBasin=candidates(:,1)>anchor.basin(1)-allowed & ...
        candidates(:,1)<anchor.basin(2)+allowed;
    candidates=candidates(inOldBasin,:);
    if isempty(candidates)
        error('ipm:TrackedPeakBranchLost', ...
            'The previous stationary-peak basin contains no continued maximum.');
    end
end
[distances,order]=sort(abs(candidates(:,1)-expected));
selected=candidates(order(1),:);
coordinateTol=1024*eps(max([1,abs(expected),abs(selected(1))]));
other=allMaxima(abs(allMaxima(:,1)-selected(1))>coordinateTol,1);
if isempty(other),separation=Inf;else,separation=min(abs(other-selected(1)));end
if distances(1)>allowed+coordinateTol || ...
        (numel(distances)>1 && distances(2)-distances(1)<=coordinateTol) || ...
        (isfinite(separation) && distances(1)>=.5*separation-coordinateTol)
    error('ipm:TrackedPeakBranchLost', ...
        ['No unambiguous stationary branch continues the predicted peak ' ...
        '(distance %.3e, allowed %.3e).'],distances(1),allowed);
end
a=selected(1);cellIndex=selected(2);t=selected(3);H=g.h(cellIndex);
minima=critical(critical(:,5)<0,1);
left=minima(minima<a);right=minima(minima>a);
basin=[window(1),window(2)];
if ~isempty(left),basin(1)=left(end);end
if ~isempty(right),basin(2)=right(1);end
if ~isempty(anchor)
    if expected<basin(1)-allowed || expected>basin(2)+allowed || ...
            min(basin(2),anchor.basin(2))<=max(basin(1),anchor.basin(1))
        error('ipm:TrackedPeakBranchLost', ...
            'Old and new stationary-peak basins do not continue one another.');
    end
end
[w,wx,wxx]=functional_rows(g,cellIndex,t);
indices=find(w~=0 | wx~=0 | wxx~=0);
P=w(indices)*q(indices)';qx=wx(indices)*q(indices)';
qxx=wxx(indices)*q(indices)';K=-qxx;
curvatureFloor=512*eps(1)*max(max(abs(y)),realmin)/H^2;
if ~isfinite(P) || P<=r.cOmegaGaugeFloor || ...
        ~isfinite(K) || K<=curvatureFloor
    error('ipm:TrackedPeakBranchLost','The selected spline branch is degenerate.');
end
[~,node]=min(abs(x-a));
neighbors=max(1,node-1):min(numel(x),node+1);
neighbors(neighbors==node)=[];
valueScale=max(sum(abs(w(indices).*q(indices))),realmin);
derivativeScale=max(sum(abs(wx(indices).*q(indices))),realmin);
if abs(qx)>1e-8*derivativeScale || ...
        abs(P-selected(4))>1024*eps(1)*max(1,valueScale)
    error('ipm:TrackedPeakBranchLost', ...
        'Independent spline evaluation and its transpose functional disagree.');
end
newAnchor=struct('version',1,'definition',definition,'x',a, ...
    'tau',r.peakEvaluationTime,'gridX',x,'window',window, ...
    'basin',basin,'cellWidth',H);
peak=struct('definition',definition,'index',node,'indices',indices, ...
    'weights',w(indices),'firstWeights',wx(indices), ...
    'secondWeights',wxx(indices),'x',a,'value',P,'curvature',qxx/2, ...
    'stationarityResidual',qx,'normalizedStationarityResidual',qx/derivativeScale, ...
    'curvatureCondition',abs(qxx/2)*H^2/max(abs(P),realmin), ...
    'weightCondition',1/sum(abs(w(indices))), ...
    'activeNodeMargin',(q(node)-max(q(neighbors)))/max(P,realmin), ...
    'vertexInteriorFraction',min(t,1-t), ...
    'gridValue',q(node),'relativeSubgridLift',(P-q(node))/P, ...
    'anchor',newAnchor,'predictedX',expected,'continuationDistance',distances(1), ...
    'continuationAllowance',allowed,'window',window,'basin',basin, ...
    'nearestCompetingMaximumDistance',separation, ...
    'cellIndex',cellIndex,'splineKnots',g.knots,'splineSampleMap',g.E);
end

function g=spline_geometry(x,window)
persistent cached
if ~isempty(cached) && isequal(cached.x,x) && isequal(cached.window,window)
    g=cached;return
end
leftTolerance=128*eps(max(1,abs(window(1))));
rightTolerance=128*eps(max(1,abs(window(2))));
inside=find(x>window(1)+leftTolerance & x<window(2)-rightTolerance);
assert(numel(inside)>=7,'ipm:DegenerateNormalization', ...
    'The fixed cubic-spline window needs at least seven interior samples.');
knots=[window(1),x(inside),window(2)];n=numel(knots);h=diff(knots);
E=sparse(n,numel(x));E(1,:)=sample_row(x,window(1));
E(end,:)=sample_row(x,window(2));
for j=1:numel(inside),E(j+1,inside(j))=1;end
A=spalloc(n-2,n-2,3*(n-2));B=spalloc(n-2,n,3*(n-2));
for j=2:n-1
    row=j-1;hl=h(j-1);hr=h(j);
    A(row,row)=2*(hl+hr);
    if row>1,A(row,row-1)=hl;end
    if row<n-2,A(row,row+1)=hr;end
    B(row,j-1:j+1)=[6/hl,-6*(1/hl+1/hr),6/hr];
end
L=chol(A,'lower');
g=struct('x',x,'window',window,'knots',knots,'h',h(:), ...
    'E',E,'A',A,'B',B,'L',L);
cached=g;
end

function row=sample_row(x,query)
row=sparse(1,numel(x));
[distance,index]=min(abs(x-query));
if distance<=128*eps(max(1,abs(query))),row(index)=1;return;end
cellIndex=find(x<=query,1,'last');
first=min(max(cellIndex-2,1),numel(x)-5);indices=first:first+5;
offset=x(indices)-query;scale=max(abs(offset));
V=(offset/scale).^((0:5)');
weights=V\[1;0;0;0;0;0];
row(indices)=weights';
end

function roots=real_quadratic_roots(a,b,c)
scale=max(abs([a,b,c]));roots=[];
if scale==0,return;end
if abs(a)<=64*eps(1)*scale
    if abs(b)>64*eps(1)*scale,roots=-c/b;end
    return
end
discriminant=b*b-4*a*c;
if discriminant<0,return;end
signed=1;if b<0,signed=-1;end
v=-.5*(b+signed*sqrt(discriminant));
if v==0,roots=-b/(2*a);else,roots=[v/a,c/v];end
end

function [w,wx,wxx]=functional_rows(g,index,t)
n=numel(g.knots);H=g.h(index);
r=zeros(n,1);r(index:index+1)=[1-t;t];
c=zeros(n,1);c(index:index+1)=H^2/6*[(1-t)^3-(1-t);t^3-t];
w=(r+g.B'*(g.L'\(g.L\c(2:end-1))))'*g.E;
r(:)=0;r(index:index+1)=[-1/H;1/H];
c(:)=0;c(index:index+1)=H/6*[-3*(1-t)^2+1;3*t^2-1];
wx=(r+g.B'*(g.L'\(g.L\c(2:end-1))))'*g.E;
c(:)=0;c(index:index+1)=[1-t;t];
wxx=(g.B'*(g.L'\(g.L\c(2:end-1))))'*g.E;
end
