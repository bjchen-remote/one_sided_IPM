function [rate,position]=transportRate(vx,vy,g)
%TRANSPORTRATE Complete mapped polar and joined-wall transport diagnostic.
assert(isequal(size(vx),size(vy),[g.nr,g.ntheta])&& ...
    all(isfinite(vx),'all')&&all(isfinite(vy),'all'), 'ipm_c:SpatialVelocity');
vr=vx.*g.cosTheta+vy.*g.sinTheta;vt=-vx.*g.sinTheta+vy.*g.cosTheta;
field=zeros(size(vx));positive=2:g.nr;
field(positive,:)=abs(vr(positive,:))./(g.radialMetric(positive)/(g.nr-1))+ ...
    abs(vt(positive,:))./(g.r(positive)*pi/(g.ntheta-1));
signedV=[flipud(vx(2:end,end));vx(:,1)]';
wallRate=abs(signedV)./(g.signedWallMetric/(g.nr-1));
field(:,1)=wallRate(g.nr:end)';field(:,end)=flipud(wallRate(1:g.nr)');
field(1,:)=wallRate(g.nr);
[rate,index]=max(field(:));[row,column]=ind2sub(size(field),index);
position=struct('r',g.r(row),'theta',g.theta(column),'X',g.X(row,column),'Y',g.Y(row,column));
assert(isfinite(rate)&&rate>=0,'ipm_c:SpatialTransportRate');
end
