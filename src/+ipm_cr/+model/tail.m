function d=tail(model,scale,X,Y,wantSecondDerivatives,sourceOnly,shared)
%DIPOLETAIL Scaled analytic vertical enrichment, before physical-time weight.
% A=K*x*y^2*rho0(x,y)/(1+x^2+y^2)^2, x=(X-S)/L,y=Y/L,
% K=-4*C/widthY^2=-2*amplitude*widthX/sqrt(pi),
% C=amplitude*widthX*widthY^2/(2*sqrt(pi)).
% H=exp(omega)*A and x/y/xx/xy are its Cartesian X/Y derivatives.
% timeH/timeX ALREADY include t'=exp(omega-ell) but not physicalTime.
% A=A_x=0 at y=0, hence centered=A and both origin corrections are zero.
% No finite differences, source truncation or tail mask is used. Gaussian
% derivatives and rational product rules are analytic. Ordinary trajectories
% use scalar exponentials; exceptional scale/Y products use combined logs.
% This combined-log guarantee covers scale, Y powers and rational factors.
% Gaussian primitive factors themselves use ordinary arithmetic; arbitrary
% rescaling of already-underflowed Gaussian exponentials is not certified.
% The unchanged no-tail gaussianBackground retains its wider log arithmetic.
% Optional sourceOnly=true returns ONLY d.x, using the identical Tx and
% scaled_field arithmetic while skipping unused background/time fields.
% Optional seventh input is the raw second output of gaussianBackground for
% the IDENTICAL model, ell, shift and X/Y queries. It is local input reuse,
% never a persistent state cache; all raw field arithmetic stays in order.
if nargin<5,wantSecondDerivatives=true;end
if nargin<6,sourceOnly=false;end
if nargin<7,shared=[];end
validateattributes(wantSecondDerivatives,{'logical'},{'scalar'});
validateattributes(sourceOnly,{'logical'},{'scalar'});
wantSecondDerivatives=wantSecondDerivatives&&~sourceOnly;
assert(isnumeric(X)&&isreal(X)&&isnumeric(Y)&&isreal(Y)&& ...
    isequal(size(X),size(Y))&&all(isfinite(X),'all')&&all(isfinite(Y),'all'), ...
    'ipm_cr:TailQueries','Use matching finite real Cartesian queries.');
assert(isstruct(model)&&strcmp(model.kind,'gaussian_primitive')&& ...
    all(isfield(model,{'amplitude','widthX','widthY','centerX','skewWeight'})), ...
    'ipm_cr:TailModel','Use the original Gaussian primitive model.');
sx=model.widthX;sy=model.widthY;amplitude=model.amplitude;w=model.skewWeight;
assert(all(isfinite([sx,sy,amplitude,model.centerX,w]))&&sx>0&&sy>0&&amplitude>0&&w>=0&&w<=1, ...
    'ipm_cr:TailModel','Invalid Gaussian primitive parameters.');
ell=scale.logC_l;omega=scale.logC_omega;
assert(all(isfinite([ell,omega,scale.X_shift])),'ipm_cr:TailScale');
% This is a representability guard on actual initial coordinates, not an
% artificial bound on log scales. In particular ell=400, omega=600 is valid.
inverse=exp(-ell);
if isfinite(inverse)&&inverse>0
    difference=X-scale.X_shift;x=difference*inverse;y=Y*inverse;
else
    difference=X-scale.X_shift;
    x=scaled_number(difference,-ell);y=scaled_number(Y,-ell);
end
overflow=~isfinite(difference);
if any(overflow,'all')
    x(overflow)=scaled_number(X(overflow),-ell)-scaled_number(scale.X_shift,-ell);
end
assert(all(isfinite(x),'all')&&all(isfinite(y),'all')&& ...
    max([abs(x(:));abs(y(:));0])<=1e150, ...
    'ipm_cr:TailCoordinateRange','Tail initial coordinates must be representable with magnitude <= 1e150.');
eta=y/sy;eta2=eta.^2;
assert(all(isfinite(eta2),'all'),'ipm_cr:TailCoordinateRange');
useShared=~isempty(shared);
if useShared
    required={'kind','modelParameters','scaleKey','X','Y','zx','z2','E1','E2','EY','erfc1','erfc2'};
    assert(isstruct(shared)&&isscalar(shared)&&all(isfield(shared,required))&& ...
        strcmp(shared.kind,'common_gaussian_raw_v1')&& ...
        isequaln(shared.modelParameters,[amplitude,model.centerX,sx,sy,w])&& ...
        isequaln(shared.scaleKey,[ell,scale.X_shift])&& ...
        isequaln(shared.X,X)&&isequaln(shared.Y,Y), ...
        'ipm_cr:TailSharedIdentity','Raw Gaussian reuse requires the identical model/scale/query.');
    for name={'zx','z2','E1','E2','EY','erfc1','erfc2'}
        assert(isequal(size(shared.(name{1})),size(X)), ...
            'ipm_cr:TailSharedShape','Raw Gaussian fields must match the query shape.');
    end
    EY=shared.EY;
else
    EY=exp(-eta2);
end
R=zeros(size(X));Rx=R;Rxx=R;
weights=[1-w,w];stretch=[1,1.25];shift=[0,.8];
for k=1:2
    if weights(k)==0,continue;end
    if useShared
        if k==1,z=shared.zx;rawE=shared.E1;rawErfc=shared.erfc1;
        else,z=shared.z2;rawE=shared.E2;rawErfc=shared.erfc2;end
        E=rawE.*EY;
        R=R+(amplitude*sx*sqrt(pi)/2*weights(k))*(rawErfc.*EY);
    else
        z=((x-model.centerX)/sx-shift(k))/stretch(k);
        E=exp(-z.^2).*EY;
        R=R+(amplitude*sx*sqrt(pi)/2*weights(k))*(erfc(-z).*EY);
    end
    Rx=Rx+(amplitude*weights(k)/stretch(k))*E;
    if wantSecondDerivatives
        ze=z.*E;ze(~isfinite(z)&E==0)=0;
        Rxx=Rxx+(-2*amplitude/sx*weights(k)/stretch(k)^2)*ze;
    end
end
% The admitted finite-coordinate range keeps D representable. Cache these
% ratios once; the ordinary trajectory performs no per-query logs here.
D=x.^2+y.^2+1;inverseD=1./D;
a2=(x.^2).*inverseD;Y2=Y.^2;
K=-2*amplitude*sx/sqrt(pi);
Tx=(1-4*a2).*R+x.*Rx;
if sourceOnly
    d=struct('x',scaled_field(K,Tx,omega-3*ell,Y,Y2,2,D,inverseD));
    d.x(Y==0)=0;
    assert(all(isfinite(d.x),'all'),'ipm_cr:TailRange', ...
        'The requested analytic tail source is not representable.');
    return
end
b2=(y.^2).*inverseD;xi=x.*inverseD;
H=x.*R;
Ty=2*x.*((1-2*b2).*R-eta2.*R);
Txy=[];Txx=[];
if wantSecondDerivatives
    Txx=12*xi.*(2*a2-1).*R+2*(1-4*a2).*Rx+x.*Rxx;
    Txy=2*(((1-4*a2).*(1-eta2)+2*b2.*(6*a2-1)).*R+ ...
        x.*((1-2*b2).*Rx-eta2.*Rx));
end
d=struct('H',scaled_field(K,H,omega-2*ell,Y,Y2,2,D,inverseD), ...
    'x',scaled_field(K,Tx,omega-3*ell,Y,Y2,2,D,inverseD), ...
    'y',scaled_field(K,Ty,omega-2*ell,Y,Y2,1,D,inverseD), ...
    'xx',[],'xy',[],'originValue',0,'originX',0, ...
    'timeH',scaled_field(K,H,2*omega-3*ell,Y,Y2,2,D,inverseD), ...
    'timeX',scaled_field(K,Tx,2*omega-4*ell,Y,Y2,2,D,inverseD));
if wantSecondDerivatives
    d.xx=scaled_field(K,Txx,omega-4*ell,Y,Y2,2,D,inverseD);
    d.xy=scaled_field(K,Txy,omega-3*ell,Y,Y2,1,D,inverseD);
end
wall=Y==0;d.H(wall)=0;d.x(wall)=0;d.y(wall)=0;
d.timeH(wall)=0;d.timeX(wall)=0;
if wantSecondDerivatives,d.xx(wall)=0;d.xy(wall)=0;end
for name={'H','x','y','xx','xy','timeH','timeX'}
    assert(all(isfinite(d.(name{1})),'all'),'ipm_cr:TailRange', ...
        'The requested analytic tail is not representable.');
end
end

function value=scaled_field(K,bracket,exponent,Y,Y2,power,D,inverseD)
% Common path: no per-point logarithm/exponential beyond the already cached
% geometry denominator. Exceptional Y^p or scale products use a joint log.
coefficient=exp(exponent);if power==2,yp=Y2;else,yp=Y;end
factor=coefficient*yp;
ordinary=isfinite(coefficient)&&coefficient>0&&all(isfinite(yp),'all')&& ...
    coefficient>=realmin&&all(Y==0|abs(yp)>=realmin,'all')&& ...
    all(isfinite(factor),'all')&&all(Y==0|abs(factor)>=realmin,'all');
if ordinary
    first=bracket.*inverseD;reduced=first.*inverseD;
    value=K*(factor.*reduced);
    % Subnormal intermediates can lose a product that later scale factors
    % would restore. Recompute only exceptional points in the joint form.
    exceptional=Y~=0&bracket~=0&(~isfinite(value)|value==0|abs(value)<realmin| ...
        abs(first)<realmin|abs(reduced)<realmin);
    if ~any(exceptional,'all'),return;end
else
    value=zeros(size(Y));exceptional=Y~=0&bracket~=0;
end
if any(exceptional,'all')
    logarithm=log(abs(K))+exponent+power*log(abs(Y(exceptional)))+ ...
        log(abs(bracket(exceptional)))-2*log(D(exceptional));
    value(exceptional)=sign(K)*sign(Y(exceptional)).^power.* ...
        sign(bracket(exceptional)).*exp(logarithm);
end
end

function value=scaled_number(number,exponent)
value=sign(number).*exp(log(abs(number))+exponent);value(number==0)=0;
end
