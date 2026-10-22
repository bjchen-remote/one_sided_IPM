function D = quadrantDerivative(x,order,parity,width)
%IPM.MESH.QUADRANTDERIVATIVE Differentiate on x>=0 with virtual parity nodes.
% Only an n-by-n sparse operator is stored. Negative coordinates exist only
% inside a local stencil; no mirrored field or negative-x grid is allocated.
validateattributes(x,{'numeric'},{'vector','real','finite','increasing'});
validateattributes(order,{'numeric'},{'scalar','integer','>=',1,'<=',2});
validateattributes(width,{'numeric'},{'scalar','integer','>=',order+1});
assert(mod(width,2)==1,'ipm:QuadrantStencilWidth', ...
    'The parity stencil width must be odd.');
assert(x(1)==0 && numel(x)>=width,'ipm:QuadrantDerivativeAxis', ...
    'The quadrant axis must start at zero and contain one full stencil.');
assert(any(strcmp(parity,{'even','odd'})),'ipm:QuadrantParity', ...
    'Parity must be even or odd.');
x=x(:);n=numel(x);virtualCount=2*n-1;half=floor(width/2);
rows=zeros(n*width,1);cols=rows;values=rows;cursor=0;
for i=1:n
    center=n+i-1;
    first=min(max(center-half,1),virtualCount-width+1);
    stencil=first:first+width-1;
    indices=abs(stencil-n)+1;
    coordinates=x(indices);
    negative=stencil<n;
    coordinates(negative)=-coordinates(negative);
    offsets=coordinates-x(i);
    scale=max(abs(offsets));
    if scale==0,error('ipm:FiniteDifferenceScale','A quadrant stencil collapsed.');end
    moments=(offsets(:)'/scale).^((0:width-1)');
    target=zeros(width,1);target(order+1)=factorial(order)/scale^order;
    coefficients=moments\target;
    if strcmp(parity,'odd')
        coefficients(negative)=-coefficients(negative);
        coefficients(indices==1)=0;
    end
    slots=cursor+(1:width);rows(slots)=i;cols(slots)=indices;
    values(slots)=coefficients;cursor=cursor+width;
end
D=sparse(rows,cols,values,n,n);
if strcmp(parity,'even') && mod(order,2)==1
    D(1,:)=0;
end
end
