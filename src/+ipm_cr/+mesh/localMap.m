function M=localMap(nodes,queries,width)
%LOCALMAP Sparse local polynomial value map on arbitrary increasing nodes.
% Barycentric weights use the actual local nodes after translation/scaling;
% Chebyshev endpoint weights are inapplicable to the stretched carrier axis.
nodes=nodes(:);queries=queries(:);n=numel(nodes);nq=numel(queries);
assert(all(isfinite(nodes))&&all(diff(nodes)>0)&&all(isfinite(queries))&& ...
    width<=n,'ipm_cr:CarrierInterpolation','Finite increasing interpolation axes are required.');
lowerTolerance=64*eps(max(1,abs(nodes(1))));
upperTolerance=64*eps(max(1,abs(nodes(end))));
assert(all(queries>=nodes(1)-lowerTolerance&queries<=nodes(end)+upperTolerance), ...
    'ipm_cr:CarrierInterpolation','Local source interpolation cannot extrapolate.');
queries=max(nodes(1),min(nodes(end),queries));
nearest=interp1(nodes,(1:n)',queries,'nearest');
first=max(1,min(n-width+1,nearest-floor(width/2)));
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
