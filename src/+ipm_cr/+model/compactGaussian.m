function v=compactGaussian(m,x,y,second)
%COMPACTGAUSSIAN Smoothly continue the primitive to its horizontal limits.
% The density can have a plateau. Its x derivative has compact support.
[wx,wx1,wx2]=cutoff(x,m.cutoff,second);
[wy,wy1]=cutoff(y,m.cutoff,false);
sx=m.widthX;sy=m.widthY;w=m.skewWeight;A=m.amplitude;
z=(x-m.centerX)/sx;z2=(z-.8)/1.25;
amplitude=A*sx*sqrt(pi)/2;
primitive=amplitude*((1-w)*erfc(-z)+w*erfc(-z2));
tail=primitive;
right=x>m.cutoff(1);
tail(right)=-amplitude*((1-w)*erfc(z(right))+w*erfc(z2(right)));
horizontal=wx.*primitive;
horizontal(right)=2*amplitude+wx(right).*tail(right);
E=exp(-(y/sy).^2);vertical=E.*wy;
v=struct('rho',horizontal.*vertical, ...
    'correction',E.*((wy-1).*primitive+wy.*((wx-1).*tail)));
e1=exp(-z.^2);e2=exp(-z2.^2);
gx=A*((1-w)*e1+(w/1.25)*e2);
horizontalX=wx.*gx+wx1.*tail;
verticalY=E.*(wy1-2*(y/sy^2).*wy);
v.x=horizontalX.*vertical;
v.y=horizontal.*verticalY;
if second
    gxx=(-2*A/sx)*((1-w)*(z.*e1)+(w/1.25^2)*(z2.*e2));
    v.xx=(wx.*gxx+2*wx1.*gx+wx2.*tail).*vertical;
    v.xy=horizontalX.*verticalY;
end
end

function [w,d1,d2]=cutoff(x,radii,second)
start=radii(1);stop=radii(2);width=stop-start;
z=(abs(x)-start)/width;
w=double(z<=0);d1=zeros(size(z));d2=[];
inside=z>0&z<1;
u=z(inside);v=1-u;
rise=exp(-1./u);fall=exp(-1./v);total=rise+fall;
a=fall./total;b=rise./total;
L1=1./u.^2+1./v.^2;
w(inside)=a;
d1(inside)=-sign(x(inside)).*(a.*b).*L1/width;
if second
    d2=zeros(size(z));
    L2=-2./u.^3+2./v.^3;
    d2(inside)=(a.*b).*((1-2*a).*L1.^2-L2)/width^2;
end
end
