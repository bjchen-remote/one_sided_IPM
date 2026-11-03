function report=test_central_transport(outputDir)
%TEST_CENTRAL_TRANSPORT Order, constant nullspace and periodic RK4 stability.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<1,outputDir=fullfile(root,'result','test_central_transport');end
if ~isfolder(outputDir),mkdir(outputDir);end
error=zeros(1,3);levels=[32,64,128];
for j=1:3
    h=1/levels(j);x=0:h:1;
    [p,m]=ipm_cr.op.central(sin(2*pi*x),h,.2);
    % The complete positive-speed split includes the explicitly chosen Q6.
    error(j)=max(abs(p-2*pi*cos(2*pi*x)));
    [p0,m0]=ipm_cr.op.central(1e12*ones(3,numel(x)),h,.2);
    assert(~any(p0,'all')&&~any(m0,'all'),'ipm_cr:CentralConstant');
    [p1,m1]=ipm_cr.op.central(x.^4,h,.2);
    assert(max(abs(p1-4*x.^3))<1e-10&&max(abs(m1-4*x.^3))<1e-10, ...
        'ipm_cr:CentralPolynomial','Boundary and interior must reproduce quartics.');
end
order=log2(error(1:2)./error(2:3));
assert(all(order>3.8),'ipm_cr:CentralOrder','Complete-array centered convergence fell below fourth order.');
n=48;h=1/n;field=repmat(eye(n),1,3);
[p,m]=ipm_cr.op.central(field,h,.2);A=-p(:,n+(1:n))';
lambda=eig(A);dt=.3*h/1.2;z=dt*lambda;
amplification=abs(1+z+z.^2/2+z.^3/6+z.^4/24);
assert(max(real(lambda))<1e-11&&max(amplification)<=1+1e-12, ...
    'ipm_cr:CentralPeriodicStability','Periodic constant-speed RK4 control is unstable.');
assert(max(abs((p(:,n+(1:n))-m(:,n+(1:n)))/2),[],'all')>0, ...
    'ipm_cr:CentralDamping','The positive sixth-order damping must be present.');
report=struct('passed',true,'levels',levels,'fullDerivativeError',error,'observedOrders',order, ...
    'constantExactlyZero',true,'quarticReproducedIncludingBoundaries',true, ...
    'maximumPeriodicRealEigenvalue',max(real(lambda)), ...
    'maximumRk4Amplification',max(amplification),'dissipationCoefficient',.2, ...
    'scope','fourth-order one-sided-smooth lines and periodic constant-speed stability; coupled patch stability requires trajectory tests');
test_report(outputDir,'central_transport',report);
end
