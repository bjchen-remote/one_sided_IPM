function [E,D,metadata]=assembly(action,a,queries,side)
%ASSEMBLY Finite-element local basis and physical point-evaluation maps.
% [E,D]=assembly('evaluate',AXIS,QUERIES,SIDE) flattens QUERIES columnwise.
% E is continuous C0 point evaluation; D differentiates the SAME FE field.
% SIDE='average' (default), 'left', or 'right' selects the derivative at
% exact internal element faces. Domain endpoints use their unique interior
% derivative for all SIDE choices. Non-face queries use their owning element.
% No extrapolation is permitted. metadata.queryShape preserves caller shape.
%
% Axis connectivity/referenceNodes may be cells for variable-degree native
% elements. [B,D]=assembly('basis',REFERENCE_NODES,XI) returns values and
% exact first derivatives of that element's actual interpolation polynomial.
if nargin<4||isempty(side),side='average';end
assert((ischar(action)&&isrow(action))||(isstring(action)&&isscalar(action)), ...
    'ipm_cr:FEMAssemblyAction','Supply a scalar assembly action.');
action=char(action);
if strcmp(action,'basis')
    [E,D]=local_basis(a,queries);
    metadata=struct('derivativeCoordinate','reference','queryShape',size(queries));return
end
assert(strcmp(action,'evaluate'),'ipm_cr:FEMAssemblyAction','Use evaluate or basis.');
assert(isstruct(a)&&isscalar(a)&&all(isfield(a, ...
    {'nodes','breaks','degree','referenceNodes','connectivity'})), ...
    'ipm_cr:FEMAxisContract','Supply an axis assembled by ipm_cr.poisson.fem.axis.');
assert(((ischar(side)&&isrow(side))||(isstring(side)&&isscalar(side)))&& ...
    any(strcmp(side,{'average','left','right'})), ...
    'ipm_cr:FEMDerivativeSide','Use average, left, or right.');
side=char(side);
validateattributes(queries,{'numeric'},{'real','finite'},mfilename,'queries');
shape=size(queries);q=double(queries(:));breaks=a.breaks(:);
ne=numel(breaks)-1;ndof=numel(a.nodes);nq=numel(q);nlocal=a.degree+1;
if isfield(a,'elementDegrees'),nlocal=max(a.elementDegrees)+1;end
assert(all(q>=breaks(1)&q<=breaks(end)),'ipm_cr:FEMQueryOutside', ...
    'Finite-element queries must lie inside the closed axis interval.');
if isempty(q)
    E=sparse(0,ndof);D=E;
    metadata=struct('queryShape',shape,'internalFace',false(0,1), ...
        'valueElement',zeros(0,1),'derivativeSide',side);return
end
element=discretize(q,breaks);
assert(all(isfinite(element)),'ipm_cr:FEMQueryOwner','Every query needs an element owner.');
% discretize assigns internal faces to their right element, the final
% domain endpoint to the last element. Faces mean exact physical breaks.
face=(element>1)&(q==breaks(element));
baseWeight=ones(nq,1);leftWeight=0;
if strcmp(side,'average'),baseWeight(face)=.5;leftWeight=.5;
elseif strcmp(side,'left'),baseWeight(face)=0;leftWeight=1;
end
nentry=nq*nlocal;ii=zeros(nentry,1);jj=ii;vb=ii;vd=ii;
cursor=0;
[~,order]=sort(element);offset=[0;cumsum(accumarray(element,1,[ne,1]))];
for e=1:ne
    rows=order(offset(e)+1:offset(e+1));if isempty(rows),continue;end
    [conn,reference]=element_data(a,e);h=breaks(e+1)-breaks(e);
    xi=2*((q(rows)-breaks(e))/h)-1;
    % Exact stored FE nodes must give exact nodal evaluation, even when
    % the inverse affine map accumulates a final rounding error.
    for k=1:numel(conn)
        xi(q(rows)==a.nodes(conn(k)))=reference(k);
    end
    [B,Dr]=local_basis(reference,xi);
    [ri,ci]=ndgrid(rows,conn);count=numel(ri);slot=cursor+(1:count);
    ii(slot)=ri(:);jj(slot)=ci(:);vb(slot)=B(:);
    derivative=bsxfun(@times,Dr*(2/h),baseWeight(rows));vd(slot)=derivative(:);
    cursor=cursor+count;
end
E=sparse(ii(1:cursor),jj(1:cursor),vb(1:cursor),nq,ndof);
D=sparse(ii(1:cursor),jj(1:cursor),vd(1:cursor),nq,ndof);
if leftWeight>0&&any(face)
    rows=find(face);leftElement=element(rows)-1;
    ni=numel(rows)*nlocal;li=zeros(ni,1);lj=li;lv=li;cursor=0;
    [~,order]=sort(leftElement);offset=[0;cumsum(accumarray(leftElement,1,[ne,1]))];
    for e=1:ne
        r=rows(order(offset(e)+1:offset(e+1)));if isempty(r),continue;end
        h=breaks(e+1)-breaks(e);[conn,reference]=element_data(a,e);
        [~,leftDerivative]=local_basis(reference,1);
        [ri,ci]=ndgrid(r,conn);count=numel(ri);slot=cursor+(1:count);
        values=repmat(leftWeight*(2/h)*leftDerivative,numel(r),1);
        li(slot)=ri(:);lj(slot)=ci(:);lv(slot)=values(:);cursor=cursor+count;
    end
    D=D+sparse(li(1:cursor),lj(1:cursor),lv(1:cursor),nq,ndof);
end
metadata=struct('queryShape',shape,'internalFace',face, ...
    'valueElement',element,'derivativeSide',side, ...
    'faceDetection','exact supplied physical break; no off-node snapping', ...
    'derivativeCoordinate','supplied physical axis');
end

function [conn,reference]=element_data(a,e)
if iscell(a.connectivity),conn=a.connectivity{e};else,conn=a.connectivity(e,:);end
if iscell(a.referenceNodes),reference=a.referenceNodes{e};else,reference=a.referenceNodes;end
conn=reshape(conn,1,[]);reference=reference(:);
end

function [B,D]=local_basis(nodes,x)
validateattributes(nodes,{'numeric'},{'vector','real','finite','nonempty'});
validateattributes(x,{'numeric'},{'real','finite'});
nodes=double(nodes(:));x=double(x(:));n=numel(nodes);
assert(n>=2&&all(diff(nodes)>0),'ipm_cr:FEMLocalNodes', ...
    'Use distinct increasing local interpolation nodes.');
B=ones(numel(x),n);D=zeros(numel(x),n);
% Product-rule recurrence differentiates the Lagrange polynomial itself.
% It remains regular at interpolation nodes and never forms a finite-
% difference derivative or a pointwise second-derivative operator.
for i=1:n
    value=ones(size(x));derivative=zeros(size(x));
    for j=1:n
        if j==i,continue;end
        denominator=nodes(i)-nodes(j);factor=(x-nodes(j))/denominator;
        derivative=derivative.*factor+value/denominator;
        value=value.*factor;
    end
    B(:,i)=value;D(:,i)=derivative;
end
end
