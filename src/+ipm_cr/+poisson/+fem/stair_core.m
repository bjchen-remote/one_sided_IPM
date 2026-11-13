function c=stair_core(c,geometry)
%STAIR_CORE Assemble native tensor elements retained by the staircase.
s=geometry.stair;x=c.xAxis;y=c.yAxis;
n=numel(c.X);A=sparse(n,n);Dx=A;Dy=A;counts=zeros(n,1);
trace=false(size(c.X));indices=reshape(1:n,size(c.X));
qmask=false(size(c.qX));mask=s.elements;
assert(isequal(size(mask),[y.elementCount,x.elementCount])&& ...
    x.quadratureOrder==s.quadOrder&&y.quadratureOrder==s.quadOrder);
[ey,ex]=find(mask);
for k=1:numel(ey)
    i=ey(k);j=ex(k);
    rows=y.connectivity{i};cols=x.connectivity{j};
    ids=indices(rows,cols);ids=ids(:);
    [mx,kx,dx,qx]=element(x,j,cols);[my,ky,dy,qy]=element(y,i,rows);
    A(ids,ids)=A(ids,ids)+kron(mx,ky)+kron(kx,my);
    Dx(ids,ids)=Dx(ids,ids)+kron(dx,speye(numel(rows)));
    Dy(ids,ids)=Dy(ids,ids)+kron(speye(numel(cols)),dy);
    counts(ids)=counts(ids)+1;qmask(qy,qx)=true;
    if i==1||~mask(i-1,j),trace(rows(1),cols)=true;end
    if i==size(mask,1)||~mask(i+1,j),trace(rows(end),cols)=true;end
    if j==1||~mask(i,j-1),trace(rows,cols(1))=true;end
    if j==size(mask,2)||~mask(i,j+1),trace(rows,cols(end))=true;end
end
assert(isequal(reshape(counts>0,size(c.X)),s.support));
trace(c.Y==0)=false;free=s.support&~trace&c.Y>0;
assert(all(s.support(geometry.active)));
c.A=(A+A')/2;c.free=find(free);c.traceIndex=find(trace);c.stairQuadrature=qmask;
W=spdiags(1./max(counts,1),0,n,n);
c.stairReturns=struct('native',true,'chart','core','stair',true, ...
    'Da',W*Dx,'Db',W*Dy);
end

function [M,K,D,q]=element(a,e,nodes)
q=(e-1)*a.quadratureOrder+(1:a.quadratureOrder);
B=a.Bq(q,nodes);d=a.Dq(q,nodes);w=a.qw(q);
M=B'*(w.*B);K=d'*(w.*d);
[~,D]=ipm_cr.poisson.fem.assembly('basis',a.referenceNodes{e},a.referenceNodes{e});
D=D/a.elementJacobian(e);
end
