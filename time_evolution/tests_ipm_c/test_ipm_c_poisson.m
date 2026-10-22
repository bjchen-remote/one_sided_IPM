function result=test_ipm_c_poisson()
% Exact finite-disk manufactured solution psi=y*(1-r^2), -Delta psi=8*y.
g=ipm_c.mesh.build(struct('nr',65,'ntheta',65,'rmax',1,'firstSpacing',.01,'domain','semicircle'));
exact=g.Y.*(1-g.X.^2-g.Y.^2);source=8*g.Y;
[f,info]=ipm_c_legacy_fixture('poisson',source,g,struct('manufactured',true,'outerPsi',zeros(1,65)));
ex=2*g.X.*g.Y;ey=1-g.X.^2-3*g.Y.^2;
psiError=max(abs(f.psi-exact),[],'all');
uError=max([max(abs(f.u1+ey),[],'all'),max(abs(f.u2+ex),[],'all')]);
% u=(-psi_y,psi_x); psi_x=-2xy, hence expected u2=-2xy.
assert(psiError<1e-8 && uError<1e-6,'polar:PoissonManufacture','Polar Poisson manufacture failed.');
[zero,zinfo]=ipm_c_legacy_fixture('poisson',zeros(size(source)),g,struct());
assert(all(zero.psi==0,'all') && all(zero.u1==0,'all') && all(zero.u2==0,'all'));
single=zeros(size(source));single(22,19)=1;
[point,pinfo]=ipm_c_legacy_fixture('poisson',single,g,struct('greenSourceTolerance',1e-10,'greenMaxSources',2048));
wr=(g.r(23)-g.r(21))/2;wt=(g.theta(20)-g.theta(18))/2;
strength=g.r(22)*wr*wt;xs=g.X(22,19);ys=g.Y(22,19);
expected=log1p(4*g.Y(end,:)*ys./((g.X(end,:)-xs).^2+(g.Y(end,:)-ys).^2))/(4*pi)*strength;
greenError=max(abs(point.psi(end,:)-expected));
assert(greenError<1e-12,'polar:GreenManufacture','Nonzero image-Green boundary failed.');
forcing=-2+2*g.X.^2+14*g.Y.^2;
exact2=g.Y.^2.*(1-g.X.^2-g.Y.^2);
[nonzeroWall,winfo]=ipm_c_legacy_fixture('poisson',forcing,g,struct('manufactured',true,'outerPsi',zeros(1,65)));
wallSourceError=max(abs(nonzeroWall.psi-exact2),[],'all');
assert(wallSourceError<1e-8,'polar:WallSourceLiftManufacture', ...
    'Poisson manufacture with nonzero forcing on the physical wall failed.');
result=struct('passed',true,'psiError',psiError,'velocityError',uError,'poisson',info, ...
    'zeroGreen',zinfo,'singleSourceGreenError',greenError,'singleSourceGreen',pinfo, ...
    'nonzeroWallSourceManufactureError',wallSourceError,'nonzeroWallSource',winfo, ...
    'nativePdeCalls',0,'actualTrajectoryAdvanced',false);
end
