function M=tensorMap(r,theta,rq,tq,width)
%TENSORMAP Sparse rowwise local interpolation of an r-by-theta field.
assert(numel(rq)==numel(tq),'ipm_cr:CarrierInterpolation','Query axes must have matching lengths.');
R=ipm_cr.carrier.localMap(r,rq,width);
T=ipm_cr.carrier.localMap(theta,tq,width);
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
    rows{j}=row(use);cols{j}=col(use);vals{j}=value(use);
end
M=sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),nq,nr*nt);
end
