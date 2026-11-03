function a=axis(input,p,qOrder,nodalLayout,outwardBoundary)
%AXIS Genuine C0 Galerkin forms on native nodes or an explicit GLL reference.
% Native INPUT is the complete transport coordinate axis. Consecutive
% p/p+1-interval elements share endpoints and introduce no new unknowns.
% Independent Gauss points are integration points only. GLL_REFERENCE takes
% element breaks instead and is retained solely for component comparisons.
% On an outward native radial boundary, target P5 uses a complete P6 last
% element. This is a source-independent FE space choice, not a flux patch.
if nargin<2||isempty(p),p=5;end
if nargin<3||isempty(qOrder),qOrder=max(8,p+2);end
if nargin<4,nodalLayout='native';end
if nargin<5,outwardBoundary=false;end
validateattributes(outwardBoundary,{'logical'},{'scalar'});
validateattributes(input,{'numeric'},{'vector','real','finite','nonempty'});
validateattributes(p,{'numeric'},{'scalar','real','finite','integer','>=',1});
validateattributes(qOrder,{'numeric'},{'scalar','real','finite','integer','>=',p+1});
nodalLayout=char(validatestring(nodalLayout,{'native','gll_reference'}));
x=double(input(:));p=double(p);qOrder=double(qOrder);
assert(numel(x)>=2&&all(diff(x)>0)&&all(isfinite(diff(x))), ...
    'ipm_cr:FEMAxisBreaks','Use a finite strictly increasing coordinate axis.');
native=strcmp(nodalLayout,'native');
[gq,wq]=gauss_legendre(qOrder);
if native
    n=numel(x)-1;ne=floor(n/p);remainder=n-ne*p;
    originalElementCount=ne;boundaryElevated=outwardBoundary&&p==5;
    if boundaryElevated&&remainder==0
        % n=5*ne can also be represented with one fewer element and five
        % P6 elements when ne>=6. No interpolation node is introduced.
        ne=ne-1;remainder=n-ne*p;
    end
    assert(ne>=1&&remainder<=ne,'ipm_cr:FEMNativeGrid', ...
        'The native axis needs enough intervals for a partition of degree %d/%d without adding nodes.',p,p+1);
    degrees=p*ones(ne,1);
    if remainder>0
        % Distribute the higher-degree elements across the whole axis.
        extra=floor(((1:remainder)-.5)*ne/remainder)+1;degrees(extra)=p+1;
    end
    [degrees,partition]=stable_partition(x,p,gq,degrees,boundaryElevated);
    edges=[1;1+cumsum(degrees)];nodes=x;breaks=x(edges);
    boundaryPolicy=struct('outwardBoundaryRequested',outwardBoundary,'enabled',boundaryElevated, ...
        'strategy','target_P5_entire_outer_element_P6_on_existing_nodes', ...
        'lastElementDegree',degrees(end),'elementCountAdjusted',ne~=originalElementCount, ...
        'nativeUnknownsAdded',0,'sourceIndependent',true);
else
    breaks=x;ne=numel(x)-1;degrees=p*ones(ne,1);
    edges=(1:p:ne*p+1)';nodes=zeros(ne*p+1,1);
    partition=struct('strategy','explicit_GLL_reference_breaks','sourceIndependent',true);
    boundaryPolicy=struct('outwardBoundaryRequested',outwardBoundary,'enabled',false, ...
        'strategy','explicit_reference_space_unchanged','lastElementDegree',p, ...
        'elementCountAdjusted',false,'nativeUnknownsAdded',0,'sourceIndependent',true);
end
assert(qOrder>=max(degrees)+1,'ipm_cr:FEMQuadrature', ...
    'Gauss order must integrate the mass products of the highest actual element degree.');
width=diff(breaks);nq=ne*qOrder;
connections=cell(ne,1);references=cell(ne,1);
qx=zeros(nq,1);qw=qx;nentry=qOrder*sum(degrees+1);
ii=zeros(nentry,1);jj=ii;vb=ii;vd=ii;cursor=0;
for e=1:ne
    conn=edges(e):edges(e+1);h=width(e);J=h/2;
    if native
        reference=2*((nodes(conn)-breaks(e))/h)-1;
        reference([1,end])=[-1,1];
    else
        reference=gll_nodes(p);local=breaks(e)+h*((reference+1)/2);
        local([1,end])=breaks(e:e+1);nodes(conn)=local;
    end
    connections{e}=conn;references{e}=reference;
    [B,D]=ipm_cr.poisson.fem.assembly('basis',reference,gq);
    rows=(e-1)*qOrder+(1:qOrder);
    qx(rows)=breaks(e)+h*((gq+1)/2);qw(rows)=J*wq;
    [ri,ci]=ndgrid(rows,conn);slot=cursor+(1:numel(ri));cursor=cursor+numel(ri);
    ii(slot)=ri(:);jj(slot)=ci(:);vb(slot)=B(:);vd(slot)=D(:)/J;
end
assert(all(diff(nodes)>0),'ipm_cr:FEMAxisResolution','FE interpolation nodes must remain distinct.');
ndof=numel(nodes);Bq=sparse(ii,jj,vb,nq,ndof);Dq=sparse(ii,jj,vd,nq,ndof);
Wq=spdiags(qw,0,nq,nq);M=Bq'*Wq*Bq;K=Dq'*Wq*Dq;
symmetry=max(norm(M-M',inf)/max(1,norm(M,inf)),norm(K-K',inf)/max(1,norm(K,inf)));
M=(M+M')/2;K=(K+K')/2;
assert(all(isfinite(nonzeros(M)))&&all(isfinite(nonzeros(K)))&&all(qw>0), ...
    'ipm_cr:FEMAxisRange','Finite-element quadrature or assembly overflowed.');
if native
    connectivity=connections;referenceNodes=references;basis='C0 native-node piecewise Lagrange';
else
    connectivity=vertcat(connections{:});referenceNodes=references{1};basis='C0 affine-element GLL Lagrange reference';
end
a=struct('nodes',nodes,'breaks',breaks,'degree',p,'elementDegrees',degrees, ...
    'actualDegreeRange',[min(degrees),max(degrees)],'nodalLayout',nodalLayout, ...
    'partitionConditioning',partition, ...
    'boundaryDegreePolicy',boundaryPolicy, ...
    'sameNativeNodes',native,'M',M,'K',K,'Bq',Bq,'Dq',Dq,'qx',qx,'qw',qw, ...
    'connectivity',{connectivity},'referenceNodes',{referenceNodes}, ...
    'referenceQuadratureNodes',gq,'referenceQuadratureWeights',wq,'quadratureOrder',qOrder, ...
    'elementCount',ne,'dofCount',ndof,'quadratureCount',nq,'elementJacobian',width/2, ...
    'basis',basis,'massLumped',false,'derivativeCoordinate','supplied physical axis', ...
    'symmetryRoundingBeforeProjection',symmetry, ...
    'faceDerivative','arithmetic average; endpoint unique interior derivative');
end

function [degrees,info]=stable_partition(x,p,gq,initial,forceOuter)
% Moving only the p+1 elements changes neither the nodes nor their number.
% Minimize the worst Gauss-point Lebesgue factor by dynamic programming;
% ties retain the most evenly distributed placement. No source/solution
% value or measured PDE error participates in this geometric choice.
ne=numel(initial);extra=sum(initial-p);oldExtra=find(initial>p);
cost=inf(ne+1,extra+1);distance=cost;choice=zeros(ne+1,extra+1);
cost(1,1)=0;distance(1,1)=0;
for e=0:ne-1
    for used=0:min(e,extra)
        if ~isfinite(cost(e+1,used+1)),continue;end
        first=1+e*p+used;
        for take=0:1
            if forceOuter&&e==ne-1&&take~=1,continue;end
            after=used+take;
            if after>extra||extra-after>ne-e-1,continue;end
            value=local_condition(x,first,p+take,gq);
            worst=max(cost(e+1,used+1),value);
            displacement=distance(e+1,used+1);
            if take,displacement=displacement+abs((e+1)-oldExtra(after));end
            previous=cost(e+2,after+1);tolerance=1e-12*max(1,worst);
            if ~isfinite(previous)||worst<previous-tolerance|| ...
                    (abs(worst-previous)<=tolerance&&displacement<distance(e+2,after+1))
                cost(e+2,after+1)=worst;distance(e+2,after+1)=displacement;
                choice(e+2,after+1)=take;
            end
        end
    end
end
degrees=p*ones(ne,1);used=extra;
for e=ne:-1:1
    take=choice(e+1,used+1);degrees(e)=p+take;used=used-take;
end
assert(used==0&&sum(degrees)==numel(x)-1,'ipm_cr:FEMNativeGrid','Native element partition failed.');
old=0;worst=0;gapRatio=0;first=1;
for e=1:ne
    [v,gap]=local_condition(x,first,degrees(e),gq);worst=max(worst,v);gapRatio=max(gapRatio,gap);
    first=first+degrees(e);
end
first=1;
for e=1:ne,old=max(old,local_condition(x,first,initial(e),gq));first=first+initial(e);end
info=struct('strategy','geometry_only_minimax_Gauss_Lebesgue_p_or_pplus1', ...
    'sourceIndependent',true,'nativeNodesMoved',false,'originalEvenDegrees',initial, ...
    'maximumQuadratureLebesgueBefore',old,'maximumQuadratureLebesgue',worst, ...
    'maximumLocalGapRatio',gapRatio,'degreePlacementChanged',~isequal(degrees,initial));
info.outerElementElevated=forceOuter;
end
function [score,gapRatio]=local_condition(x,first,degree,gq)
values=x(first:first+degree);h=values(end)-values(1);
reference=2*(values-values(1))/h-1;reference([1,end])=[-1,1];
B=ipm_cr.poisson.fem.assembly('basis',reference,gq);
score=norm(B,inf);gap=diff(reference);gapRatio=max(gap)/min(gap);
end

function x=gll_nodes(p)
n=p-1;
if n==0,x=[-1;1];return;end
if n==1,interior=0;
else
    k=(1:n-1)';off=sqrt(k.*(k+2)./((2*k+1).*(2*k+3)));
    interior=sort(eig(diag(off,1)+diag(off,-1)));interior=(interior-flipud(interior))/2;
end
x=[-1;interior;1];
end
function [x,w]=gauss_legendre(n)
k=(1:n-1)';off=k./sqrt(4*k.^2-1);[V,L]=eig(diag(off,1)+diag(off,-1));
[x,order]=sort(diag(L));V=V(:,order);w=2*(V(1,:).^2)';
end
