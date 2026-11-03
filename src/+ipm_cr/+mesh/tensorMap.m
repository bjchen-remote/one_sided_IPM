function M=tensorMap(r,theta,rq,tq,width,extrapolate)
%TENSORMAP Sparse rowwise local interpolation of an r-by-theta field.
% Optional endpoint polynomial continuation is explicit; sampling and source
% maps keep the default no-extrapolation contract.
if nargin<6,extrapolate=false;end
assert(numel(rq)==numel(tq),'ipm_cr:CarrierInterpolation','Query axes must have matching lengths.');
R=ipm_cr.mesh.localMap(r,rq,width,extrapolate);
T=ipm_cr.mesh.localMap(theta,tq,width,extrapolate);
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
