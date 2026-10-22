function report=test_polar_cartesian()
%TEST_POLAR_CARTESIAN Manufactured FE reproduction, masks and physical units.
r=linspace(0,2,35)';theta=linspace(0,pi,32);X=r*cos(theta);Y=r*sin(theta);
X(:,end)=-r;Y(:,[1,end])=0;
d=struct('r',r,'theta',theta,'scale',struct('logC_l',0,'logC_omega',0,'X_shift',0));
errors=zeros(1,4);
for degree=0:3
    [d.rho,d.rho_X,d.rho_Y]=polynomial(X,Y,degree);
    c=ipm_c.cartesian(d,[-.8,.8,0,1.6],101,struct('fields',{{'rho'}}));
    exact=polynomial(c.X,c.Y,degree);errors(degree+1)=max(abs(c.rho-exact),[],'all');
    assert(all(c.mask,'all')&&errors(degree+1)<3e-10,'ipm_c:TestCubicReproduction', ...
        'P3 reproduction failed at degree %d: %.6g.',degree,errors(degree+1));
    assert(c.metadata.centroidReconstruction.quadraticTriangles==0&& ...
        c.metadata.centroidReconstruction.constantFallbackTriangles==0);
end
[d.rho,d.rho_X,d.rho_Y]=polynomial(X,Y,1);
c=ipm_c.cartesian(rmfield(d,{'rho_X','rho_Y'}),[-.8,.8,0,1.6],81,struct('fields',{{'rho'}}));
exact=polynomial(c.X,c.Y,1);inferredGradientError=max(abs(c.rho-exact),[],'all');
assert(inferredGradientError<1e-7,'ipm_c:TestInferredNodalGradients');
c=ipm_c.cartesian(d,[-.8,.8,0,1.6],[85,73],struct('fields',{{'rho'}},'method','linear'));
exact=polynomial(c.X,c.Y,1);linearError=max(abs(c.rho-exact),[],'all');
assert(linearError<2e-14&&isequal(size(c.rho),[73,85]));
c=ipm_c.cartesian(d,[-2.2,2.2,-.2,4.2],111,struct('fields',{{'rho'}}));
strictOutside=c.Y<0|hypot(c.X,c.Y)>r(end);
safeInside=c.Y>=0&hypot(c.X,c.Y)<r(end)*cos(max(diff(theta))/2)*(1-1e-12);
assert(~any(c.mask(strictOutside))&&all(c.mask(safeInside))&& ...
    all(isnan(c.rho(~c.mask))),'ipm_c:TestDomainMask','Mask extrapolates or has interior holes.');
c=ipm_c.cartesian(d,[-1e-8,1e-8,0,2e-8],5,struct('fields',{{'rho'}}));
originError=abs(c.rho(1,3)-d.rho(1,1));
originVariation=max(abs(c.rho-d.rho(1,1)),[],'all');
assert(originError==0&&originVariation<5e-8,'ipm_c:TestOriginContinuity');
d.scale=struct('logC_l',log(2),'logC_omega',log(3),'X_shift',.4);
d.u1=2+0*X;d.u2=-1+0*X;d.psi=.5+0*X;d.rho_tau=.8+.2*X-.1*Y;
d.selected=struct('tau',1.2,'t',.7,'scale',d.scale,'rates',struct('c_l',.2,'c_omega',-.1,'c_r',.7));
% Trajectory vectors must never replace selected scalar gauge rates.
d.c_l=[.1,.2];d.c_omega=[-.05,-.1];d.c_r=[.6,.7];
c=ipm_c.cartesian(d,[-.3,.3,0,.6],71,struct('coordinates','physical', ...
    'method','linear','fields',{{'rho','rho_X','rho_Y','u1','u2','psi','rho_t','rho_tau'}}));
nX=2*c.X+.4;nY=2*c.Y;[exact,ex,ey]=polynomial(nX,nY,1);
physicalRhoError=max(abs(c.rho-exact/3),[],'all');
physicalGradientError=max(abs(c.rho_X-(2/3)*ex),[],'all');
expectedT=(2/9)*(.8+.2*nX-.1*nY+(.2*nX+.7).*ex+.2*nY.*ey+.1*exact);
physicalTimeError=max(abs(c.rho_t-expectedT),[],'all');
assert(all(c.mask,'all')&&physicalRhoError<2e-14&&physicalGradientError<2e-14&& ...
    physicalTimeError<2e-14&&max(abs(c.u1-2/3),[],'all')<2e-14&& ...
    max(abs(c.psi-1/12),[],'all')<2e-14,'ipm_c:TestPhysicalUnits');
assert(contains(c.metadata.fieldUnits.rho_tau,'canonical'));
report=struct('status','passed','cubicPolynomialDegree',0:3,'cubicPolynomialErrors',errors, ...
    'linearAffineError',linearError,'originError',originError,'originVariation',originVariation, ...
    'inferredGradientAffineError',inferredGradientError, ...
    'physicalRhoError',physicalRhoError,'physicalGradientError',physicalGradientError, ...
    'physicalEulerianTimeDerivativeError',physicalTimeError,'strictOutsideMasked',true, ...
    'interiorHasNoHoles',true,'rectangularResolutionPreserved',true,'sharedOriginVertex',true);
end

function [f,fx,fy]=polynomial(x,y,degree)
f=2.5+0*x;fx=0*x;fy=0*x;
if degree>=1,f=f+.3*x-.7*y;fx=fx+.3;fy=fy-.7;end
if degree>=2,f=f+.4*x.^2-.1*x.*y+.2*y.^2;fx=fx+.8*x-.1*y;fy=fy-.1*x+.4*y;end
if degree>=3
    f=f+.15*x.^3-.2*x.^2.*y+.13*x.*y.^2-.11*y.^3;
    fx=fx+.45*x.^2-.4*x.*y+.13*y.^2;
    fy=fy-.2*x.^2+.26*x.*y-.33*y.^2;
end
end
