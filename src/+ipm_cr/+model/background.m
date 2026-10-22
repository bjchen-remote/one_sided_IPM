function b=background(model,scale,X,Y,wantSecondDerivatives)
%BACKGROUND Gaussian primitive plus an optional analytic physical-time tail.
% The unchanged Gaussian implementation owns its original extreme-scale
% arithmetic. With tailCorrection=true, G=exp(omega)*(rho0+t*A), where
% A is a smooth vertical dipole-tail enrichment of rho0. This is a
% representation change; timeH/timeX explicitly expose its time dependence.
if nargin<5,wantSecondDerivatives=true;end
tailEnabled=isfield(model,'tailCorrection')&&model.tailCorrection;
shared=[];
if tailEnabled
    [b,shared]=ipm_cr.model.gaussian(model,scale,X,Y,wantSecondDerivatives);
else
    b=ipm_cr.model.gaussian(model,scale,X,Y,wantSecondDerivatives);
end
b.timeH=zeros(size(X));b.timeX=zeros(size(X));
if ~tailEnabled,return;end
assert(strcmp(model.kind,'gaussian_primitive'),'ipm_cr:TailBackgroundModel', ...
    'Tail correction requires the original Gaussian primitive.');
[high,low]=physical_time(scale);
d=ipm_cr.model.tail(model,scale,X,Y,wantSecondDerivatives,false,shared);
b.H=b.H+high*d.H+low*d.H;
b.x=b.x+high*d.x+low*d.x;
b.y=b.y+high*d.y+low*d.y;
% A and A_x vanish on the physical wall, so neither origin field changes.
if wantSecondDerivatives
    b.xx=b.xx+high*d.xx+low*d.xx;
    b.xy=b.xy+high*d.xy+low*d.xy;
end
% Already includes t'=exp(omega-ell), the canonical-to-physical clock rate.
b.timeH=d.timeH;b.timeX=d.timeX;
origin=X==0&Y==0;b.H(origin)=0;b.timeH(origin)=0;b.y(Y==0)=0;
if wantSecondDerivatives,b.xy(Y==0)=0;end
for name={'H','x','y','xx','xy','originValue','originX','timeH','timeX'}
    assert(all(isfinite(b.(name{1})),'all'),'ipm_cr:TailBackgroundRange', ...
        'The tail background is not representable in the admitted range.');
end
end

function [high,low]=physical_time(scale)
assert(isfield(scale,'physicalTime'),'ipm_cr:TailBackgroundClock', ...
    'Tail correction requires the physical-time clock.');
high=scale.physicalTime;low=0;
if isfield(scale,'physicalTimeLow'),low=scale.physicalTimeLow;end
validateattributes(high,{'numeric'},{'scalar','real','finite'});
validateattributes(low,{'numeric'},{'scalar','real','finite'});
end
