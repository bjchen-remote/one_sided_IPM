function c=build(cfg)
%BUILD Cached separable Cartesian core, no two-dimensional factorization.
[c.x,c.Dx]=ipm_cr.internal.cheb(cfg.coreDegree(1),-cfg.coreHalfWidth,cfg.coreHalfWidth);
[c.y,c.Dy]=ipm_cr.internal.cheb(cfg.coreDegree(2),0,cfg.coreHeight);
c.Dxx=c.Dx*c.Dx;c.Dyy=c.Dy*c.Dy;
[c.X,c.Y]=meshgrid(c.x,c.y);c.nx=numel(c.x);c.ny=numel(c.y);
c.ix=2:c.nx-1;c.iy=2:c.ny-1;
mask=false(c.ny,c.nx);mask(:,[1,end])=true;mask(end,:)=true;mask(1,:)=false;
c.traceIndex=find(mask);c.traceX=c.X(mask);c.traceY=c.Y(mask);
[c.Vx,Lx]=eig(-c.Dxx(c.ix,c.ix));[c.Vy,Ly]=eig(-c.Dyy(c.iy,c.iy));
assert(isreal(Lx)&&isreal(Ly)&&all(diag(Lx)>0)&&all(diag(Ly)>0), ...
    'ipm_cr:CoreSpectrum','Unexpected Dirichlet core spectrum.');
c.Vxi=c.Vx\eye(size(c.Vx));c.Vyi=c.Vy\eye(size(c.Vy));
c.denominator=diag(Ly)+diag(Lx)';
c.eigenvectorCondition=max(cond(c.Vx),cond(c.Vy));
end
