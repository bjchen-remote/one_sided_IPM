function centroid=centroid(values,ops)
%WALLGRADIENTCENTROID Positive-quadrature centroid of wall R_X squared.
% On fixed [-L,L], q is sampled by piecewise-linear interpolation at the
% clipped endpoints and unchanged interior nodes. Positive trapezoidal
% weights w give E=sum(w*q^2), C=sum(w*x*q^2)/E, and the exact derivative of
% this discrete observable Cdot=sum(2*w*(x-C)*q*qdot)/E. Endpoint sampling
% is included in the returned global-node functional. No peak curvature,
% local branch, time-step, or basin condition is used.
x=ops.x(:)';q=values(:)';L=1;
if isfield(ops,'rescaling') && isfield(ops.rescaling,'wallIntegralRadius')
    L=ops.rescaling.wallIntegralRadius;
end
assert(numel(q)==numel(x) && numel(x)>=2 && all(isfinite(x)) && ...
    all(isfinite(q)) && all(diff(x)>0), ...
    'ipm:WallCentroidInput','Centroid samples need a finite increasing grid.');
assert(isnumeric(L) && isscalar(L) && isfinite(L) && L>0 && ...
    x(1)<=-L && x(end)>=L, ...
    'ipm:WallCentroidWindow','The fixed centroid window must fit on the wall.');
inside=find(x>-L & x<L);z=[-L,x(inside),L];n=numel(z);
B=sparse(n,numel(x));B(1,:)=sample_row(x,-L);B(end,:)=sample_row(x,L);
for j=1:numel(inside),B(j+1,inside(j))=1;end
h=diff(z);w=[h(1),h(1:end-1)+h(2:end),h(end)]/2;
sampled=(B*q')';qScale=max(abs(sampled));
value=0;weights=zeros(1,numel(x));energyWeights=weights;normalizedEnergy=0;
if qScale>0
    scaled=sampled/qScale;
    normalizedEnergy=sum(w.*scaled.^2);
    value=sum(w.*z.*scaled.^2)/normalizedEnergy;
    weights=(2*(w.*(z-value).*scaled)/(qScale*normalizedEnergy))*B;
    energyWeights=(2*qScale*(w.*scaled))*B;
end
indices=find(weights~=0 | energyWeights~=0);
reference=value;referencePersisted=false;
if isfield(ops,'rescaling') && isfield(ops.rescaling,'referenceWallGradientCentroid')
    reference=ops.rescaling.referenceWallGradientCentroid;
    referencePersisted=true;
end
centroid=struct('definition','fixed_window_wall_R_X_squared_centroid', ...
    'value',value,'energy',(qScale*qScale)*normalizedEnergy, ...
    'normalizedEnergy',normalizedEnergy,'normalizationScale',qScale, ...
    'degenerateEnergy',qScale==0,'radius',L,'window',[-L,L], ...
    'indices',indices,'weights',weights(indices), ...
    'energyWeights',energyWeights(indices), ...
    'samplePoints',z,'sampleValues',sampled, ...
    'quadratureWeights',w,'sampleMap',B, ...
    'referenceValue',reference,'referencePersisted',referencePersisted, ...
    'quadratureOrder',2);
end

function row=sample_row(x,query)
row=sparse(1,numel(x));index=find(x==query,1);
if ~isempty(index),row(index)=1;return;end
left=find(x<query,1,'last');
t=(query-x(left))/(x(left+1)-x(left));
row(left:left+1)=[1-t,t];
end
