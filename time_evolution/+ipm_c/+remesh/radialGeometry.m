function g=radialGeometry(old,r)
%RADIALGEOMETRY Rebuild every polar spatial operator after a radial proposal.
% Angular nodes, physical radius and node counts stay fixed. No field transfer
% or acceptance is performed by this geometry-only constructor.
r=r(:);nr=old.nr;nt=old.ntheta;
assert(numel(r)==nr&&all(isfinite(r))&&r(1)==0&&r(end)==old.r(end)&& ...
    all(diff(r)>0),'ipm_c:AreaRemeshAxis','A radial proposal must retain count/origin/domain and be strictly increasing.');
theta=old.theta(:)';
assert(max(abs(theta-linspace(0,pi,nt)))<=128*eps(pi), ...
    'ipm_c:AreaRemeshTheta','The angular axis must remain uniform.');
g=old;g.r=r;g.rmax=r(end);g.firstSpacing=r(2);
if ~isequaln(r,old.r)&&isfield(g,'uniformCoreInfo')
    g=rmfield(g,'uniformCoreInfo');
end
g.R=repmat(r,1,nt);g.Theta=repmat(theta,nr,1);
g.X=r*old.cosTheta;g.Y=r*old.sinTheta;
g.Dr=ipm.mesh.fdMatrix(r,1,old.stencilWidth);
g.Drr=ipm.mesh.fdMatrix(r,2,old.stencilWidth);
logical=linspace(0,1,nr)';g.radialMetric=ipm.mesh.fdMatrix(logical,1,7)*r;
g.signedWallX=[-flipud(r(2:end));r]';wallLogical=linspace(-1,1,2*nr-1);
g.signedWallMetric=(ipm.mesh.fdMatrix(wallLogical,1,7)*g.signedWallX')';
assert(all(isfinite(g.radialMetric))&&all(g.radialMetric>0)&& ...
    all(isfinite(g.signedWallMetric))&&all(g.signedWallMetric>0), ...
    'ipm_c:AreaRemeshMetric','The rebuilt radial and joined-wall mapped metrics must stay positive.');
rows=old.originFit.rows;fitScale=r(rows(end));
x=g.X(rows,:)/fitScale;y=g.Y(rows,:)/fitScale;x=x(:);y=y(:);
powers=old.originFit.coefficientPowers;design=zeros(numel(x),size(powers,1));
for k=1:size(powers,1),design(:,k)=x.^powers(k,1).*y.^powers(k,2);end
[Q,R]=qr(design,0);
assert(rcond(R)>1e-10,'ipm_c:AreaRemeshOriginFit','The physical-origin polynomial fit must remain well conditioned.');
g.originFit=old.originFit;g.originFit.scale=fitScale;g.originFit.map=R\Q';g.originFit.rankCondition=rcond(R);
g.radialMap='accepted_radial_transfer';g.radialMapParameter=NaN;
g.originIsSinglePhysicalPoint=true;
end
