function b=background(model,scale,X,Y,wantSecondDerivatives)
%BACKGROUND Gaussian primitive plus an optional analytic physical-time tail.
% The unchanged Gaussian implementation owns its original extreme-scale
% arithmetic. With tailCorrection=true, G=exp(omega)*(rho0+t*A), where
% A is a smooth vertical dipole-tail enrichment of rho0. This is a
% representation change; timeH/timeX explicitly expose its time dependence.
if nargin<5,wantSecondDerivatives=true;end
tailEnabled=isfield(model,'tailCorrection')&&model.tailCorrection;
shared=[];
if strcmp(model.kind,'compact_gaussian_source')
    b=compact_background(model,scale,X,Y,wantSecondDerivatives);
elseif tailEnabled
    [b,shared]=gaussian(model,scale,X,Y,wantSecondDerivatives);
else
    b=gaussian(model,scale,X,Y,wantSecondDerivatives);
end
b.timeH=zeros(size(X));b.timeX=zeros(size(X));
if ~tailEnabled,return;end
assert(strcmp(model.kind,'gaussian_primitive'),'ipm_cr:TailBackgroundModel', ...
    'Tail correction requires the original Gaussian primitive.');
[high,low]=physical_time(scale);
d=tail(model,scale,X,Y,wantSecondDerivatives,shared);
b.H=b.H+high*d.H+low*d.H;
b.x=b.x+high*d.x+low*d.x;
b.y=b.y+high*d.y+low*d.y;
if isfield(scale,'physicalTimeTail')
    for word=reshape(scale.physicalTimeTail,1,[])
        if word~=0,b.H=b.H+word*d.H;b.x=b.x+word*d.x;b.y=b.y+word*d.y;end
    end
end
% A and A_x vanish on the physical wall, so neither origin field changes.
if wantSecondDerivatives
    b.xx=b.xx+high*d.xx+low*d.xx;
    b.xy=b.xy+high*d.xy+low*d.xy;
    if isfield(scale,'physicalTimeTail')
        for word=reshape(scale.physicalTimeTail,1,[])
            if word~=0,b.xx=b.xx+word*d.xx;b.xy=b.xy+word*d.xy;end
        end
    end
end
% Already includes t'=exp(omega-ell), the canonical-to-physical clock rate.
b.timeH=d.timeH;b.timeX=d.timeX;
origin=X==0&Y==0;b.H(origin)=0;b.timeH(origin)=0;b.y(Y==0)=0;
if wantSecondDerivatives,b.xy(Y==0)=0;end
end

function [high,low]=physical_time(scale)
assert(isfield(scale,'physicalTime'),'ipm_cr:TailBackgroundClock', ...
    'Tail correction requires the physical-time clock.');
high=scale.physicalTime;low=0;
if isfield(scale,'physicalTimeLow'),low=scale.physicalTimeLow;end
validateattributes(high,{'numeric'},{'scalar','real'});
validateattributes(low,{'numeric'},{'scalar','real'});
end

function result=compact_background(model,scale,X,Y,second)
%COMPACTBACKGROUND Keep the original stable Gaussian arithmetic inside.
base=model;base.kind='gaussian_primitive';
ell=scale.logC_l;omega=scale.logC_omega;inverse=exp(-ell);
x=(X-scale.X_shift)*inverse;y=Y*inverse;
outside=abs(x)>model.cutoff(1)|abs(y)>model.cutoff(1);
result=gaussian(base,scale,X,Y,second);
origin=-scale.X_shift*inverse;
if any(outside,'all')
    value=ipm_cr.model.initial.compactGaussian(model,x(outside),y(outside),second);
    result.H(outside)=result.H(outside)+exp(omega)*value.correction;
    result.x(outside)=exp(omega-ell)*value.x;
    result.y(outside)=exp(omega-ell)*value.y;
    if second
        result.xx(outside)=exp(omega-2*ell)*value.xx;
        result.xy(outside)=exp(omega-2*ell)*value.xy;
    end
end
if abs(origin)>model.cutoff(1)
    v=ipm_cr.model.initial.compactGaussian(model,origin,0,false);
    result.H=result.H-exp(omega)*v.correction;
    result.originValue=exp(omega)*v.rho;
    result.originX=exp(omega-ell)*v.x;
end
result.H(X==0&Y==0)=0;
end

function [b,shared]=gaussian(model,scale,X,Y,wantSecondDerivatives)
%BACKGROUND Analytic gauge-advected Gaussian primitive and centered values.
% G=exp(logC_omega)*rho0((X-X_shift)/exp(logC_l),Y/exp(logC_l)).
% H=G-G(0,0) is evaluated as an increment, never by subtracting two G's.
% X/Y can have any matching finite real shape. Derivatives are Cartesian.
% Optional WANTSECONDDERIVATIVES defaults to true; false returns empty xx/xy
% and leaves all other outputs bitwise identical to the complete path.
% This routine changes no state and makes no spatial-order assertion.
% Optional second output shares raw ordinary-range transcendental results
% only for this model/ell/shift/query. Empty means the original log fallback.
shared=[];
if nargin<5,wantSecondDerivatives=true;end
validateattributes(wantSecondDerivatives,{'logical'},{'scalar'});
assert(isnumeric(X)&&isreal(X)&&isnumeric(Y)&&isreal(Y)&& ...
    isequal(size(X),size(Y))&&all(isfinite(X),'all')&&all(isfinite(Y),'all'), ...
    'ipm_cr:BackgroundQueries','Use matching finite real Cartesian queries.');
assert(isstruct(model)&&isscalar(model)&&isfield(model,'kind'), ...
    'ipm_cr:BackgroundModel','Supply one immutable background model.');
kind=char(model.kind);
z=zeros(size(X));
b=struct('H',z,'x',z,'y',z,'xx',z,'xy',z,'originValue',0,'originX',0);
if ~wantSecondDerivatives,b.xx=[];b.xy=[];end
if strcmp(kind,'none'),return;end
assert(strcmp(kind,'gaussian_primitive'),'ipm_cr:BackgroundModel', ...
    'The background kind must be none or gaussian_primitive.');
required={'amplitude','centerX','widthX','widthY','skewWeight'};
assert(all(isfield(model,required)),'ipm_cr:BackgroundModel', ...
    'The Gaussian background requires all original analytic parameters.');
for k=1:numel(required)
    validateattributes(model.(required{k}),{'numeric'},{'scalar','real','finite'});
end
assert(model.amplitude>0&&model.widthX>0&&model.widthY>0&& ...
    model.skewWeight>=0&&model.skewWeight<=1,'ipm_cr:BackgroundModel', ...
    'Use positive amplitude/widths and a skew weight in [0,1].');
assert(isstruct(scale)&&isscalar(scale)&& ...
    all(isfield(scale,{'logC_l','logC_omega','X_shift'})), ...
    'ipm_cr:BackgroundScale','Supply the current complete RK stage scale.');
for name={'logC_l','logC_omega','X_shift'}
    validateattributes(scale.(name{1}),{'numeric'},{'scalar','real'});
end

ell=scale.logC_l;omega=scale.logC_omega;
if abs(ell)<=80&&abs(omega)<=80
    if nargout>1
        [ready,b,shared]=common_background(model,scale,X,Y,b,wantSecondDerivatives,true);
    else
        [ready,b]=common_background(model,scale,X,Y,b,wantSecondDerivatives);
    end
    if ready,return;end
end
logSx=log(model.widthX);logSy=log(model.widthY);
logA=log(model.amplitude);logC=logA+logSx+log(pi)/2-log(2);
% The increment has a separate logarithm even when X/exp(ell) underflows.
% Keeping it separate also avoids losing X before subtraction from a large S.
logDx=log(abs(X))-ell-logSx;signDx=sign(X);
baseX=-scaled_number(scale.X_shift,-ell);
z0=(baseX-model.centerX)/model.widthX;
% For ordinary queries z0+delta is sufficient. The direct coordinate also
% handles X near X_shift when the two individually scaled terms overflow.
physicalDifference=X-scale.X_shift;
initialX=scaled_number(physicalDifference,-ell);
overflow=~isfinite(physicalDifference);
if any(overflow,'all')
    initialX(overflow)=scaled_number(X(overflow),-ell) ...
        -scaled_number(scale.X_shift,-ell);
end
zx=(initialX-model.centerX)/model.widthX;
logY=log(abs(Y));logEta2=2*(logY-ell-logSy);
eta2=exp(logEta2);
% log(1-exp(-eta^2)) retains eta^2 even below the subnormal range.
logVerticalDifference=logEta2;
ordinary=logEta2>-36;
logVerticalDifference(ordinary)=log(-expm1(-eta2(ordinary)));

weights=[1-model.skewWeight,model.skewWeight];
shifts=[0,.8];stretches=[1,1.25];
logP=-Inf(size(X));logP0=-Inf;
for k=1:2
    if weights(k)==0,continue;end
    w=log(weights(k));s=stretches(k);ls=log(s);
    t=(zx-shifts(k))/s;t0=(z0-shifts(k))/s;
    logP=log_add(logP,w+log_erfc_minus(t));
    logP0=log_add(logP0,w+log_erfc_minus(t0));
    % C*exp(-eta^2)*[erfc(-t)-erfc(-t0)], with the independent
    % signed increment X/(L*sx*s). All products use combined exponents.
    b.H=b.H+scaled_erf_increment(t0,t,signDx,logDx-ls, ...
        logC+omega+w-eta2);
    squared=t.^2;
    exponent=logA+omega-ell+w-ls-squared-eta2;
    b.x=b.x+exp(exponent);
    originExponent=logA+omega-ell+w-ls-t0^2;
    b.originX=b.originX+exp(originExponent);
    if wantSecondDerivatives
        finite=isfinite(t)&t~=0;
        term=zeros(size(X));
        term(finite)=-sign(t(finite)).*exp(log(2)+logA+omega-2*ell ...
            +w-logSx-2*ls+log(abs(t(finite)))-squared(finite)-eta2(finite));
        b.xx=b.xx+term;
        % Form G_xy directly instead of multiplying an underflowed G_x.
        b.xy=b.xy-sign(Y).*exp(log(2)+logY+logA+omega-3*ell ...
            +w-ls-2*logSy-squared-eta2);
    end
end
b.H=b.H-exp(logC+omega+logP0+logVerticalDifference);
b.y=-sign(Y).*exp(log(2)+logY+logC+omega-2*ell-2*logSy+logP-eta2);
b.originValue=exp(logC+omega+logP0);
% These identities hold independently of subtraction/underflow arithmetic.
origin=X==0&Y==0;b.H(origin)=0;b.y(Y==0)=0;
if wantSecondDerivatives,b.xy(Y==0)=0;end
end

function [ready,b,shared]=common_background(model,scale,X,Y,b,wantSecondDerivatives,wantShared)
% Ordinary trajectory range: reuse Gaussian exponentials and scalar scales.
% H still uses an independent coordinate increment and expm1; the fallback
% below remains responsible for extreme combined exponents/coordinates.
if nargin<7,wantShared=false;end
shared=[];
ready=false;ell=scale.logC_l;omega=scale.logC_omega;
inverse=exp(-ell);sx=model.widthX;sy=model.widthY;w=model.skewWeight;
z0=(-scale.X_shift*inverse-model.centerX)/sx;
delta=(X*inverse)/sx;
zx=((X-scale.X_shift)*inverse-model.centerX)/sx;
eta=(Y*inverse)/sy;
C0=model.amplitude*sx*sqrt(pi)/2;
C=C0*exp(omega);
qx=model.amplitude*exp(omega-ell);
qxx=(-2*model.amplitude/sx)*exp(omega-2*ell);
qy=(-2*C0/sy^2)*exp(omega-2*ell);
qxy=(-2*model.amplitude/sy^2)*exp(omega-3*ell);
coefficients=[C,qx,qxx,qy,qxy];
if ~isfinite(z0)||abs(z0)>1e150||any(~isfinite(coefficients))|| ...
        any(coefficients==0)||any(abs(coefficients)>1e150)||abs(C)>1e60|| ...
        any(abs(X)>1e150|abs(Y)>1e150,'all')||any(~isfinite(zx),'all')|| ...
        any(abs(zx)>1e150,'all')||any(~isfinite(delta),'all')|| ...
        any(abs(delta)>1e150,'all')||any(~isfinite(eta),'all')|| ...
        any(abs(eta)>1e150,'all')||any(X~=0&delta==0,'all')|| ...
        any(delta~=0&abs(delta)<realmin,'all')|| ...
        any(Y~=0&eta.^2==0,'all')
    return;
end
z2=(zx-.8)/1.25;z02=(z0-.8)/1.25;
E1=exp(-zx.^2);E2=exp(-z2.^2);EY=exp(-eta.^2);
if wantShared
    erfc1=erfc(-zx);erfc2=erfc(-z2);
    P=(1-w)*erfc1+w*erfc2;
else
    P=(1-w)*erfc(-zx)+w*erfc(-z2);
end
P0=(1-w)*erfc(-z0)+w*erfc(-z02);
increment=(1-w)*erf_increment(z0,zx,delta) ...
    +w*erf_increment(z02,z2,delta/1.25);
gaussian=(1-w)*E1+(w/1.25)*E2;
b.H=C*(increment.*EY+P0*expm1(-eta.^2));
b.x=qx*(gaussian.*EY);
b.y=qy*(Y.*P.*EY);
if wantSecondDerivatives
    b.xx=qxx*(((1-w)*zx.*E1+(w/1.25^2)*z2.*E2).*EY);
    b.xy=qxy*(Y.*gaussian.*EY);
end
b.originValue=C*P0;
b.originX=qx*((1-w)*exp(-z0^2)+(w/1.25)*exp(-z02^2));
origin=X==0&Y==0;b.H(origin)=0;b.y(Y==0)=0;
if wantSecondDerivatives,b.xy(Y==0)=0;end
if wantShared
    shared=struct('kind','common_gaussian_raw_v1', ...
        'modelParameters',[model.amplitude,model.centerX,model.widthX,model.widthY,model.skewWeight], ...
        'scaleKey',[ell,scale.X_shift],'X',X,'Y',Y, ...
        'zx',zx,'z2',z2,'E1',E1,'E2',E2,'EY',EY,'erfc1',erfc1,'erfc2',erfc2);
end
ready=true;
end

function value=erf_increment(a,endpoint,delta)
% Twelve scalar integral coefficients evaluated by Horner on a small delta.
% The guard makes every product small. For abs(a)>30 the whole small-step
% Gaussian tail is below double range after any admitted common-range scale.
value=zeros(size(delta));small=abs(delta)*(1+abs(a))<=.075;
if any(small,'all')&&abs(a)<=30
    coefficients=zeros(1,12);coefficients(1)=1;coefficients(2)=-2*a;
    for n=1:10
        coefficients(n+2)=(-2*a*coefficients(n+1)-2*coefficients(n))/(n+1);
    end
    d=delta(small);polynomial=coefficients(12)/12+zeros(size(d));
    for n=11:-1:1,polynomial=polynomial.*d+coefficients(n)/n;end
    value(small)=(2/sqrt(pi)*exp(-a^2))*(d.*polynomial);
end
remaining=~small;
positive=remaining&(a>=0)&(endpoint>=0);
negative=remaining&(a<=0)&(endpoint<=0);
cross=remaining&~positive&~negative;
value(positive)=erfc(a)-erfc(endpoint(positive));
value(negative)=erfc(-endpoint(negative))-erfc(-a);
value(cross)=erf(endpoint(cross))-erf(a);
end

function value=scaled_number(x,exponent)
% Signed x*exp(exponent), without an overflowing/underflowing scale factor.
value=sign(x).*exp(log(abs(x))+exponent);
value(x==0)=0;
end

function value=log_erfc_minus(x)
% log(erfc(-x)); erfcx protects the Gaussian primitive's negative tail.
value=zeros(size(x));negative=x<0;
value(negative)=log(erfcx(-x(negative)))-x(negative).^2;
value(~negative)=log(erfc(-x(~negative)));
end

function value=log_add(a,b)
% Positive log-sum-exp, including the exact zero represented by -Inf.
high=max(a,b);low=min(a,b);
value=high+log1p(exp(low-high));
value(isinf(high)&high<0)=-Inf;
end

function value=scaled_erf_increment(a,endpoint,sgn,logDelta,exponent)
% exp(exponent)*(erf(a+delta)-erf(a)), with delta kept independently.
% Twelve integral-Taylor terms give a tiny remainder for
% abs(delta)*(1+abs(a))<=.075. Scaled coefficient recurrence avoids powers
% of a overflowing when delta is very small. No per-query quadrature is used.
value=zeros(size(endpoint));nonzero=sgn~=0;
small=nonzero&logDelta+log1p(abs(a))<=log(.075);
if any(small,'all')
    ld=logDelta(small);
    ad=sign(a).*sgn(small).*exp(log(abs(a))+ld);
    d2=exp(2*ld);previous=ones(size(ld));current=-2*ad;
    total=previous+current/2;
    for n=1:10
        next=(-2*ad.*current-2*d2.*previous)/(n+1);
        total=total+next/(n+2);previous=current;current=next;
    end
    value(small)=sgn(small).*exp((exponent(small)-a^2) ...
        +log(2/sqrt(pi))+ld+log(total));
end
remaining=nonzero&~small;
if ~any(remaining,'all'),return;end
same=remaining&((a>=0&endpoint>=0)|(a<=0&endpoint<=0));
if any(same,'all')
    ld=logDelta(same);e=endpoint(same);
    delta=sgn(same).*exp(ld);
    squareDifference=2*a.*delta+delta.^2;
    first=log(erfcx(abs(a)));second=log(erfcx(abs(e)))-squareDifference;
    correction=max(first,second)+log(-expm1(-abs(first-second)));
    candidate=(exponent(same)-a^2)+correction;
    % Extreme coordinates may overflow squared terms even though the
    % primitive difference itself is bounded. Fall back to endpoint tails.
    extreme=~isfinite(squareDifference)|~isfinite(a^2)|isnan(candidate);
    if any(extreme)
        lo=min(abs(a),abs(e(extreme)));hi=max(abs(a),abs(e(extreme)));
        logLo=log(erfcx(lo))-lo.^2;logHi=log(erfcx(hi))-hi.^2;
        difference=logLo+log(-expm1(logHi-logLo));
        difference(isinf(logLo)&logLo<0)=-Inf;
        original=exponent(same);candidate(extreme)=original(extreme)+difference;
    end
    value(same)=sgn(same).*exp(candidate);
end
cross=remaining&~same;
if any(cross,'all')
    magnitude=erf(abs(a))+erf(abs(endpoint(cross)));
    value(cross)=sgn(cross).*exp(exponent(cross)+log(magnitude));
end
end

function d=tail(model,scale,X,Y,wantSecondDerivatives,shared)
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
% Optional sixth input is the raw second output of gaussianBackground for
% the IDENTICAL model, ell, shift and X/Y queries. It is local input reuse,
% never a persistent state cache; all raw field arithmetic stays in order.
if nargin<5,wantSecondDerivatives=true;end
if nargin<6,shared=[];end
validateattributes(wantSecondDerivatives,{'logical'},{'scalar'});
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
inverse=exp(-ell);
if isfinite(inverse)&&inverse>0
    difference=X-scale.X_shift;x=difference*inverse;y=Y*inverse;
else
    difference=X-scale.X_shift;
    x=tail_scaled_number(difference,-ell);y=tail_scaled_number(Y,-ell);
end
overflow=~isfinite(difference);
if any(overflow,'all')
    x(overflow)=tail_scaled_number(X(overflow),-ell)-tail_scaled_number(scale.X_shift,-ell);
end
eta=y/sy;eta2=eta.^2;
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
% Cache the geometry ratios once; the ordinary trajectory performs no
% per-query logs here.
D=x.^2+y.^2+1;inverseD=1./D;
a2=(x.^2).*inverseD;Y2=Y.^2;
K=-2*amplitude*sx/sqrt(pi);
Tx=(1-4*a2).*R+x.*Rx;
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

function value=tail_scaled_number(number,exponent)
value=sign(number).*exp(log(abs(number))+exponent);value(number==0)=0;
end
