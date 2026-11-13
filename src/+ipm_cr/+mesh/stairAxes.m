function [out,info]=stairAxes(action,varargin)
%STAIRAXES Smooth point allocation on the existing physical coordinate maps.
switch action
    case 'core',[out,info]=core(varargin{:});
    case 'radial',[out,info]=radial(varargin{:});
    otherwise,error('ipm_cr:StairAxes','Unknown staircase axis operation.');
end
end

function [a,info]=core(a,scale,fine)
% The reference dense rectangle ends at 35% of each positive logical axis.
% A scalar scales the 2D share isotropically. A pair scales the X1 and X2
% interval shares independently. Each axis keeps its total interval budget.
cut=.35;
if isscalar(scale)
    validateattributes(scale,{'numeric'},{'real','finite','>=',.1,'<=',1});
    directionScale=sqrt(scale)*[1,1];
else
    validateattributes(scale,{'numeric'},{'vector','numel',2,'real','finite','>=',sqrt(.1),'<=',1.8});
    directionScale=reshape(scale,1,2);
end
extra=1.35*double(fine);baseTotal=1+extra*(1-cut);
m=(numel(a.x)+1)/2;
[x,dx,xcut]=redistribute(a.x(m:end),a.dxMetric(m:end),directionScale(1));
[y,dy,ycut]=redistribute(a.y,a.dyMetric,directionScale(2));
a.x=[-fliplr(x(2:end)),x];a.dxMetric=[fliplr(dx(2:end)),dx];
a.y=y(:);a.dyMetric=dy(:);
info=struct('targetRelativeCoreShare',prod(directionScale),'coreBounds',[xcut,ycut], ...
    'targetAxisScale',directionScale, ...
    'referenceLogicalCut',cut,'referenceAxisShare',cut/baseTotal, ...
    'targetAxisShare',directionScale*cut/baseTotal);

    function [v,d,physicalCut]=redistribute(v,d,axisScale)
        n=numel(v)-1;newN=round(n*baseTotal);
        if axisScale<=1
            outer=baseTotal*(1-axisScale)/(axisScale-baseTotal*cut^2);
            total=baseTotal+outer;
        else
            % Add a positive monitor over the core and taper it to zero.
            % This avoids driving the exterior node density toward zero.
            share=axisScale*cut/baseTotal;
            coreMass=1.5*cut;
            gain=(axisScale-1)*cut/(cut-share*coreMass);
            total=baseTotal+gain*coreMass;
        end
        target=linspace(0,total,newN+1);
        lo=zeros(size(target));hi=ones(size(target));
        for k=1:46
            s=(lo+hi)/2;
            mass=s+extra*max(0,s-cut).^3/(1-cut)^2;
            if axisScale<=1
                mass=mass+outer*s.^3;
            else
                mass=mass+gain*core_weight(s,cut);
            end
            left=mass<target;lo(left)=s(left);hi(~left)=s(~left);
        end
        s=(lo+hi)/2;s([1,end])=[0,1];
        [physicalCut,~]=hermite(v,d,n*cut);
        [v,d]=hermite(v,d,n*s);
        density=1+3*extra*max(0,s-cut).^2/(1-cut)^2;
        if axisScale<=1
            density=density+3*outer*s.^2;
        else
            [~,weight]=core_weight(s,cut);
            density=density+gain*weight;
        end
        d=d.*(n*total/newN)./density;
    end
end

function [mass,weight]=core_weight(s,cut)
t=max(0,min(1,(s-cut)/cut));
mass=s-cut*t.^4.*(2.5-3*t+t.^2)-max(0,s-2*cut);
weight=1-t.^3.*(10-15*t+6*t.^2);
end

function [r,dr]=radial(r,dr,interfaceRadius)
% Double resolution near the interface with a smooth transition to the
% original far axis. In reversed indices the enriched region is at the end.
n=numel(r)-1;
first=find(r>=1.18*interfaceRadius,1);
stop=max(1,n+1-first);start=max(0,stop-24);
span=stop-start;
added=ceil(n-(start+stop)/2);gain=added/(n-(start+stop)/2);
target=0:n+added;
lo=zeros(size(target));hi=n+lo;
for k=1:52
    z=(lo+hi)/2;
    t=max(0,min(1,(z-start)/span));
    mass=z+gain*(span*t.^4.*(2.5-3*t+t.^2)+max(0,z-stop));
    left=mass<target;lo(left)=z(left);hi(~left)=z(~left);
end
z=(lo+hi)/2;z(target<=start)=target(target<=start);z(end)=n;
t=max(0,min(1,(z-start)/span));
[rev,d]=hermite(-flip(r(:)),flip(dr(:)),z);
d=d./(1+gain*t.^3.*(10-15*t+6*t.^2));
r=-flip(rev(:));dr=flip(d(:));
end

function [v,d]=hermite(v,d,z)
% Interpolate physical coordinates using the supplied metric as index slope.
v=v(:).';d=d(:).';n=numel(v)-1;
j=min(floor(z),n-1)+1;t=z-(j-1);
gap=v(j+1)-v(j);m0=d(j);m1=d(j+1);
c2=3*gap-2*m0-m1;c3=m0+m1-2*gap;
first=v(1);last=v(end);
v=v(j)+t.*(m0+t.*(c2+t.*c3));v(z==0)=first;v(z==n)=last;
d=m0+t.*(2*c2+3*t.*c3);
end
