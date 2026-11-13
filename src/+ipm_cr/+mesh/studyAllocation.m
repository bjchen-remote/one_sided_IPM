function [out,info]=studyAllocation(a,plan,halo)
%STUDYALLOCATION Reallocate a coarse thin-staircase budget on native axes.
% The first halo/ring coordinates and the box endpoints are fixed.
validateattributes(plan,{'double'},{'vector','real','finite'});
assert(ismember(numel(plan),[6,7,8]));
nearOffset=0;
if numel(plan)>=7
    nearOffset=plan(7);
    validateattributes(nearOffset,{'double'},{'integer','>=',-1,'<=',2});
end
plan=reshape(plan,1,[]);
% Reuse one allocation in physical coordinates across the paired presets.
% The 181-angle medium chart has twice the base coordinate intervals.
refinement=1+double(numel(a.theta)==181);
addedPairs=refinement*plan(1);removedRings=refinement*plan(2);center=plan(3);
halfWidth=plan(4);warpStrength=plan(5);nearRings=refinement*plan(6);
validateattributes([addedPairs,removedRings],{'double'},{'integer','>=',-8,'<=',30});
validateattributes(nearRings,{'double'},{'integer','nonnegative','<=',30});
assert(center>0&&halfWidth>0&&warpStrength>=0&&warpStrength<1);
a.r=a.r(:).';a.drMetric=a.drMetric(:).';
cap=min(a.x(end-halo),a.y(end-halo));
capIndex=find(a.r<=cap,1,'last');
assert(capIndex>halo,'ipm_cr:StudyGeometry','The inner ring must remain inside the donor box.');
second=a.r(capIndex);
baseCounts=[numel(a.x),numel(a.y),numel(a.r)-capIndex+halo+2,numel(a.theta)];
tailAdjustment=0;
if refinement==2
    % Count the coarse partner in the same physical frame. The three
    % ghost rows make the cropped sizes differ from the full-axis ratio.
    coarseR=a.r(end:-2:1);
    coarseCap=min(a.x(end-2*halo),a.y(end-2*halo));
    coarseNr=nnz(coarseR>coarseCap)+halo+2-plan(2);
    tailAdjustment=2*(coarseNr-1)+1-(baseCounts(3)-removedRings);
end
middle=(numel(a.x)+1)/2;x=a.x(middle:end);dx=a.dxMetric(middle:end);
start=max(x(2),(center-halfWidth)*second);
finish=min(x(end-halo-2),(center+halfWidth)*second);
assert(finish>start,'ipm_cr:StudyGeometry','The target allocation interval is empty.');
focus=min(finish,max(start,center*second));
if warpStrength>0
    width=finish-start;t=(x-start)/width;active=t>0&t<1;
    z=t(active);
    bump=64*(z.*(1-z)).^3;
    bumpPrime=192*(z.*(1-z)).^2.*(1-2*z)/width;
    offset=x(active)-focus;
    x(active)=x(active)-warpStrength*offset.*bump;
    dx(active)=dx(active).*(1-warpStrength*(bump+offset.*bumpPrime));
end
if addedPairs>0
    z=interp1(x,0:numel(x)-1,[start,finish],'linear');
    [x,dx]=redistribute(x,dx,[z(1),z(2),addedPairs],addedPairs);
elseif addedPairs<0
    % Keep the staircase and its active Cartesian columns inside the
    % unchanged prefix. Spend fewer columns on the unused tensor tail.
    tailStart=interp1(x,0:numel(x)-1,1.12*second,'linear');
    tailEnd=numel(x)-1;
    assert(tailEnd-tailStart>1.875*abs(addedPairs),'ipm_cr:StudyTailDensity', ...
        'The Cartesian tail must have positive coordinate density.');
    [x,dx]=redistribute(x,dx,[tailStart,tailEnd,addedPairs],addedPairs);
end
y=a.y;dy=a.dyMetric;
if numel(plan)==8
    % Refine the measured corner height with a C2 physical-coordinate map.
    strength=plan(8);
    validateattributes(strength,{'double'},{'scalar','nonnegative','<',1});
    lo=.06*second;hi=min(.44*second,y(end-halo-1));focusY=.24*second;
    t=(y-lo)/(hi-lo);inside=t>0&t<1;z=t(inside);
    bump=64*(z.*(1-z)).^3;
    bumpPrime=192*(z.*(1-z)).^2.*(1-2*z)/(hi-lo);
    offset=y(inside)-focusY;
    y(inside)=y(inside)-strength*offset.*bump;
    dy(inside)=dy(inside).*(1-strength*(bump+offset.*bumpPrime));
end
out=struct('x',[-fliplr(x(2:end)),x],'dxMetric',[fliplr(dx(2:end)),dx], ...
    'y',y,'dyMetric',dy,'r',a.r,'drMetric',a.drMetric,'theta',a.theta);
if removedRings~=0||nearRings>0||tailAdjustment~=0
    n=numel(a.r)-1;
    windows=[capIndex+refinement*(2+nearOffset), ...
        capIndex+refinement*(15+nearOffset),nearRings; ...
        capIndex+refinement*20,n,tailAdjustment-removedRings-nearRings];
    [out.r,out.drMetric]=redistribute(a.r,a.drMetric,windows,tailAdjustment-removedRings);
end
counts=baseCounts+[2*addedPairs,0,tailAdjustment-removedRings,0];
nominal=@(c)c(1)*c(2)+c(3)*c(4);
assert(all(diff(out.x)>0)&&all(out.dxMetric>0)&& ...
    all(diff(out.r)>0)&&all(out.drMetric>0)&&all(diff(out.y)>0)&&all(out.dyMetric>0));
newCap=min(out.x(end-halo),out.y(end-halo));
assert(isequal(newCap,cap)&&isequal(out.r(1:capIndex+1),a.r(1:capIndex+1)));
assert(isequal(out.x([1,end]),a.x([1,end]))&&out.r(end)==a.r(end));
measure=struct('plan',plan,'targetBand',[start,finish],'secondActiveRadius',second, ...
    'referenceSizes',baseCounts,'sizes',counts,'referenceNominalCount',nominal(baseCounts), ...
    'nominalCount',nominal(counts),'allocationRefinement',refinement);
info=struct('referenceX',a.x,'referenceDx',a.dxMetric,'measure',measure);
if numel(plan)==8
    info.referenceY=a.y;
    info.referenceDy=a.dyMetric;
end
end

function [v,d]=redistribute(v,d,windows,delta)
% Integrate a C2 quintic cumulative allocation and invert in index space.
n=numel(v)-1;target=0:n+delta;
lo=zeros(size(target));hi=n+lo;
for iteration=1:47
    mid=(lo+hi)/2;
    [mapped,~]=mass(mid,windows);
    below=mapped<target;
    lo(below)=mid(below);hi(~below)=mid(~below);
end
z=(lo+hi)/2;
% Restore every integer sample outside the windows exactly, including caps.
anchor=0:n;[position,~]=mass(anchor,windows);
fixed=position==round(position);
inside=false(size(anchor));
for j=1:size(windows,1)
    inside=inside|(anchor>windows(j,1)&anchor<windows(j,2));
end
keep=fixed&~inside;
z(position(keep)+1)=anchor(keep);
[~,density]=mass(z,windows);
j=min(floor(z),n-1)+1;t=z-j+1;
gap=v(j+1)-v(j);m0=d(j);m1=d(j+1);
c2=3*gap-2*m0-m1;c3=m0+m1-2*gap;
samples=v(j)+t.*(m0+t.*(c2+t.*c3));
samples(z==n)=v(end);
d=(m0+t.*(2*c2+3*t.*c3))./density;
v=samples;
end

function [value,derivative]=mass(z,windows)
value=z;derivative=ones(size(z));
for j=1:size(windows,1)
    row=windows(j,:);width=row(2)-row(1);
    t=min(1,max(0,(z-row(1))/width));
    value=value+row(3)*t.^3.*(10-15*t+6*t.^2);
    derivative=derivative+row(3)/width*30*t.^2.*(1-t).^2;
end
end
