function [coeff,info]=corner(source,g,analyticSource)
%CORNER_POLYNOMIAL Cartesian Hermite jets for the two outer-wall corners.
% qP=sum coeff(a+1,b+1)*(x/R)^a*(y/R)^b. Its 20 coefficients
% match the interpolated total-order <=3 Taylor jets at (+/-R,0).
% Source interpolation is degree6 in normalized radial coordinate and angle;
% no lower-order wall derivative or clearing of a physical wall source.
if nargin<3,analyticSource=[];end
R=g.polar.r(end);ri=numel(g.polar.r)-6:numel(g.polar.r);
assert(numel(g.polar.r)>=7&&numel(g.polar.theta)>=7,'ipm_cr:FEMCornerNodes');
rho=g.polar.r(ri)/R-1;rs=max(abs(rho));
Vr=vander7(rho/rs);coeff=zeros(8,4);jets=zeros(4,4,2);
for side=1:2
    e=3-2*side;
    if e==1,ti=1:7;phi=g.polar.theta(ti);else,ti=numel(g.polar.theta):-1:numel(g.polar.theta)-6;phi=pi-g.polar.theta(ti);end
    ts=max(abs(phi));Vt=vander7(phi/ts);
    x=g.polar.X(ri,ti);y=g.polar.Y(ri,ti);
    if isa(source,'function_handle'),q=source(x,y);else,q=source.polar(ri,ti);end
    if ~isempty(analyticSource),q=q+analyticSource(x,y);end
    assert(isequal(size(q),size(x))&&isreal(q),'ipm_cr:FEMCornerSource', ...
        'Corner source values must be real and match the supplied stencil shape.');
    anchor=q(end,1);a=(Vr\(q-anchor))/Vt';a(1,1)=a(1,1)+anchor;
    a=a(1:4,1:4)./(rs.^(0:3)')./(ts.^(0:3));
    % rho=e*dx+dy^2/2-e*dx*dy^2/2, phi=dy-e*dx*dy+dx^2*dy-dy^3/3.
    rp=zeros(4);rp(2,1)=e;rp(1,3)=.5;rp(2,3)=-e/2;
    tp=zeros(4);tp(1,2)=1;tp(2,2)=-e;tp(3,2)=1;tp(1,4)=-1/3;
    jet=zeros(4);rpow=zeros(4);rpow(1)=1;
    for i=0:3
        tpow=zeros(4);tpow(1)=1;
        for j=0:3-i
            jet=jet+a(i+1,j+1)*mul(rpow,tpow);tpow=mul(tpow,tp);
        end
        rpow=mul(rpow,rp);
    end
    jets(:,:,side)=jet;
end
for b=0:3
    order=3-b;degree=2*order+1;H=zeros(degree+1);rhs=zeros(degree+1,1);row=0;
    for side=1:2
        e=3-2*side;
        for i=0:order
            row=row+1;rhs(row)=jets(i+1,b+1,side);
            for a=i:degree,H(row,a+1)=nchoosek(a,i)*e^(a-i);end
        end
    end
    coeff(1:degree+1,b+1)=H\rhs;
end
info=struct('R',R,'fit','degree6 normalized-radius/angle interpolation; Cartesian total-order3 corner jets', ...
    'sourceCoefficientCount',20,'cornerTaylorJets',jets,'coefficientInfinity',max(abs(coeff),[],'all'), ...
    'cornerSourceValues',[jets(1,1,1),jets(1,1,2)],'jetAccuracyClaim','formal fourth-order for third jets of a smooth source; measured separately');
end
function V=vander7(x)
x=x(:);V=ones(numel(x),7);for j=2:7,V(:,j)=V(:,j-1).*x;end
end
function C=mul(A,B)
C=zeros(4);for i=0:3,for j=0:3-i,for k=0:3-i,for l=0:3-i-j-k
    C(i+k+1,j+l+1)=C(i+k+1,j+l+1)+A(i+1,j+1)*B(k+1,l+1);
end,end,end,end
end
