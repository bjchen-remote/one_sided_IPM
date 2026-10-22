function B=bary(nodes,queries)
%BARY Stable second-form interpolation on ascending Lobatto nodes.
% Only called on affine Chebyshev-Lobatto axes, including panel log-r axes.
x=nodes(:);q=queries(:);n=numel(x)-1;
w=(-1).^(0:n);w([1,end])=w([1,end])/2;
d=q-x';hit=abs(d)<=8*eps(max(1,max(abs(x))));
d(hit)=1;B=w./d;B=B./sum(B,2);
rows=find(any(hit,2));
for k=rows',B(k,:)=0;B(k,find(hit(k,:),1))=1;end
end
