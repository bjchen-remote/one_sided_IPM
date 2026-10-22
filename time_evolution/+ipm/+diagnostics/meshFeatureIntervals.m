function feature=meshFeatureIntervals(view)
%IPM.DIAGNOSTICS.MESHFEATUREINTERVALS Paired positive-wall core/front geometry.
% Uses an already accepted source and derivative matrix, never an elliptic
% solve or gauge-rate feedback. The caller retains the same field/axis units.
required={'rho','x','y','Dx','source','trusted'};
assert(isstruct(view)&&isscalar(view)&&all(isfield(view,required)), ...
    'ipm:AutonomousMeshView','A complete accepted rho/axes/Dx/source view is required.');
x=view.x(:)';y=view.y(:);rho=view.rho;omega=view.source;
validateattributes(x,{'numeric'},{'real','finite','increasing','numel',size(rho,2)});
validateattributes(y,{'numeric'},{'real','finite','increasing','numel',size(rho,1)});
assert(islogical(view.trusted)&&isscalar(view.trusted)&&view.trusted && ...
    isequal(size(omega),size(rho))&&isequal(size(view.Dx),[numel(x),numel(x)])&& ...
    numel(x)>=17&&numel(y)>=8&&y(1)==0&& ...
    (isequal(x,-fliplr(x)) || x(1)==0)&& ...
    all(isfinite(rho),'all')&&all(isfinite(omega),'all')&&isreal(rho)&&isreal(omega), ...
    'ipm:AutonomousMeshView','The accepted view must be finite and symmetric or quadrant-only.');
assert(isequal(omega,rho*view.Dx'), ...
    'ipm:AutonomousMeshSourcePair','Accepted source must equal rho*Dx'' exactly.');
if x(1)==0
    wall=smooth_signal(abs(omega(1,:)));
    envelope=smooth_signal(max(abs(omega),[],1));
else
    wall=smooth_signal(.5*(abs(omega(1,:))+fliplr(abs(omega(1,:)))));
    envelope=max(abs(omega),[],1);envelope=smooth_signal(.5*(envelope+fliplr(envelope)));
end
frontSignal=smooth_signal(abs((view.Dx*wall(:))'));
tolerance=100*eps(max(1,max(abs(x))));positive=x>=-tolerance;
xp=x(positive);xp(abs(xp)<=tolerance)=0;wallPositive=wall(positive);envelopePositive=envelope(positive);
[value,frontPeakIndex]=max(wallPositive);
if value<=100*eps(max(1,max(envelopePositive))),[value,frontPeakIndex]=max(envelopePositive);end
assert(value>0,'ipm:AutonomousMeshFeature','A nonzero positive-half source feature is required.');
[frontCenter,frontHalfWidth]=front_geometry(xp,frontSignal(positive),1,max(frontPeakIndex-1,1),xp(frontPeakIndex));
rawWall=max(omega(1,positive),0);
horizontal=ipm.diagnostics.peakResolution(rawWall,xp,[.1,.5,.9]);
assert(horizontal.peak>0&&xp(horizontal.peakIndex)>0, ...
    'ipm:AutonomousMeshFeature','A positive wall peak away from the origin is required.');
coreCenter=xp(horizontal.peakIndex);bounds=horizontal.bounds(3,:);
left=max(coreCenter-bounds(1),local_spacing(xp,horizontal.peakIndex));
right=max(bounds(2)-coreCenter,local_spacing(xp,horizontal.peakIndex));
indices=find(positive);peakColumn=indices(horizontal.peakIndex);
vertical=ipm.diagnostics.peakResolution(abs(omega(:,peakColumn))',y,[.1,.5,.9]);
yWidth=max(vertical.widths(3),local_spacing(y,vertical.peakIndex));
core=[max(0,coreCenter-left),min(x(end),coreCenter+right)];
front=[max(0,frontCenter-frontHalfWidth),min(x(end),frontCenter+frontHalfWidth)];
feature=struct('coreInterval',core,'frontInterval',front,'yCoreWidth',yWidth,'coreCenter',coreCenter, ...
    'actualCoreCells',[horizontal.gridPoints(3),vertical.gridPoints(3)], ...
    'leftFrontCells',interval_count(xp,front),'peakColumn',peakColumn, ...
    'rawCoreBounds',bounds,'sourcePairedExactly',true,'positivePeakSelection','largest positive wall nodal peak', ...
    'nodalVerticalColumnDefinition',true,'continuousWidthClaim',false,'gaugeRatesModified',false);
end

function v=smooth_signal(v)
span=min(5,numel(v));if mod(span,2)==0,span=span-1;end
if span>1,v=movmean(v,span,'Endpoints','shrink');end
end

function [center,halfWidth]=front_geometry(s,indicator,first,last,fallback)
if last<first||isempty(indicator(first:last))||max(indicator(first:last))<=0
    center=fallback;halfWidth=local_spacing(s,min(max(first,1),numel(s)));return
end
[peak,relative]=max(indicator(first:last));index=first+relative-1;center=s(index);
left=crossing_left(s,indicator,index,.5*peak);right=crossing_right(s,indicator,index,.5*peak);
halfWidth=max(.5*(right-left),local_spacing(s,index));
end

function crossing=crossing_left(s,signal,index,level)
crossing=s(1);
for probe=index-1:-1:1
    if signal(probe)<=level,crossing=linear_crossing(s(probe),s(probe+1),signal(probe),signal(probe+1),level);return;end
end
end

function crossing=crossing_right(s,signal,index,level)
crossing=s(end);
for probe=index+1:numel(s)
    if signal(probe)<=level,crossing=linear_crossing(s(probe-1),s(probe),signal(probe-1),signal(probe),level);return;end
end
end

function x=linear_crossing(x1,x2,y1,y2,target)
if y2==y1,x=.5*(x1+x2);else,fraction=(target-y1)/(y2-y1);fraction=min(max(fraction,0),1);x=x1+fraction*(x2-x1);end
end

function spacing=local_spacing(s,index)
if index<=1,spacing=s(2)-s(1);elseif index>=numel(s),spacing=s(end)-s(end-1);else,spacing=.5*(s(index+1)-s(index-1));end
end

function count=interval_count(axis,interval)
v=axis(:)';count=sum(max(0,min(v(2:end),interval(2))-max(v(1:end-1),interval(1)))./diff(v));
end
