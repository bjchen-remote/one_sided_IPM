function [values,diagnostic]=sampleState(state,g,method,q)
%SAMPLESTATE Read complete physical fields at a fixed common query cloud.
% STATE must contain an accepted residual and its corresponding FLOW; this
% function does not advance, align, project, or modify either transport grid.
% Background terms are evaluated analytically at the queries, while only
% residuals and native velocities are interpolated. In particular, F and Q
% include the physical-time background terms timeH and timeX.
%
% Q is the full comparison-cloud structure: finite column vectors x/y and
% matching logical masks near, interfaceBand, wall, remoteWallStrip, hotspot.
% Each region must contain a query. Reported infinity norms apply to these
% finite samples and are not continuum bounds over the corresponding region.
%
% VALUES preserves the IVP study contract: W, q, Q, F, u1Centered, u2, rates,
% scalars. Uref remains a separate rate; it is subtracted after interpolating
% physical u1. Scalar clockHigh/clockLow retain the compensated clock parts.
% DIAGNOSTIC compares reconstruction with interpolation of stored full arrays.
required={'x','y','near','interfaceBand','wall','remoteWallStrip','hotspot'};
assert(isstruct(q)&&isscalar(q)&&all(isfield(q,required)), ...
    'ipm_cr:SampleStateQueries','Use a complete comparison query cloud.');
assert(isnumeric(q.x)&&isreal(q.x)&&iscolumn(q.x)&&~isempty(q.x)&& ...
    isnumeric(q.y)&&isreal(q.y)&&isequal(size(q.x),size(q.y))&& ...
    all(isfinite(q.x))&&all(isfinite(q.y)), ...
    'ipm_cr:SampleStateQueries','Use matching finite real query columns.');
for name=required(3:end)
    mask=q.(name{1});
    assert(islogical(mask)&&isequal(size(mask),size(q.x))&&any(mask), ...
        'ipm_cr:SampleStateQueries','Every comparison region requires a nonempty matching logical mask.');
end

b=ipm_cr.explicit.background(method.background,state.scale,q.x,q.y);
eta=ipm_cr.explicit.interpolate(state.residual,g,q.x,q.y);
d=ipm_cr.explicit.differentiate(state.residual,g);
qr=ipm_cr.explicit.synchronize(struct('core',d.core.x,'polar',d.polar.x),g);
d=ipm_cr.explicit.differentiate(state.flow.residualRhs,g);
Qr=ipm_cr.explicit.synchronize(struct('core',d.core.x,'polar',d.polar.x),g);
f=state.flow;Ht=f.c_omega*b.H-(f.c_l*q.x+f.c_r).*b.x-f.c_l*q.y.*b.y+f.c_r*b.originX;
Htx=(f.c_omega-f.c_l)*b.x-(f.c_l*q.x+f.c_r).*b.xx-f.c_l*q.y.*b.xy;
% These explicit-background time terms already include dt_physical/dtau.
Ht=Ht+b.timeH;Htx=Htx+b.timeX;
u1=struct('core',f.physical.core.u1,'polar',f.physical.polar.u1);
u2=struct('core',f.physical.core.u2,'polar',f.physical.polar.u2);
values=struct('W',b.H+eta,'q',b.x+ipm_cr.explicit.interpolate(qr,g,q.x,q.y), ...
    'Q',Htx+ipm_cr.explicit.interpolate(Qr,g,q.x,q.y), ...
    'F',Ht+ipm_cr.explicit.interpolate(f.residualRhs,g,q.x,q.y), ...
    'u1Centered',ipm_cr.explicit.interpolate(u1,g,q.x,q.y)-f.Uref, ...
    'u2',ipm_cr.explicit.interpolate(u2,g,q.x,q.y), ...
    'rates',rate_values(state),'scalars',scalar_values(state));
diagnostic=struct('rawFullArrayReinterpolationDifference',struct());
for name={'W','q','Q'}
    key=name{1};field=state.rho;if ~strcmp(key,'W'),field=f.(key);end
    raw=ipm_cr.explicit.interpolate(field,g,q.x,q.y);
    diagnostic.rawFullArrayReinterpolationDifference.(key)=region_metrics(raw-values.(key),q);
end
end

function info=region_metrics(d,q)
assert(all(isfinite(d)),'ipm_cr:IVPQueries');
info=struct('global',metric(d,true(size(d)),q),'near',metric(d,q.near,q), ...
    'interfaceBand',metric(d,q.interfaceBand,q),'wall',metric(d,q.wall,q), ...
    'remoteWallStrip',metric(d,q.remoteWallStrip,q),'bareNullspaceHotspots',metric(d,q.hotspot,q));
end

function info=metric(d,mask,q)
indices=find(mask);[value,k]=max(abs(d(indices)));j=indices(k);
info=struct('sampleCount',numel(indices),'absoluteInfinity',value,'rms',sqrt(mean(d(indices).^2)), ...
    'maximumLocation',struct('X',q.x(j),'Y',q.y(j),'r',hypot(q.x(j),q.y(j))));
end

function v=rate_values(s)
v=struct('names',{{'c_l','c_omega','c_r','delta','Uref','transportOrigin'}}, ...
    'values',[s.flow.c_l,s.flow.c_omega,s.flow.c_r,s.flow.delta,s.flow.Uref,s.flow.transportOrigin]);
end

function v=scalar_values(s)
v=struct('names',{{'logC_l','logC_omega','physicalTimeCompensated','X_shift','densityOrigin'}}, ...
    'values',[s.scale.logC_l,s.scale.logC_omega,s.scale.physicalTime+s.scale.physicalTimeLow,s.scale.X_shift,s.scale.densityOrigin], ...
    'clockHigh',s.scale.physicalTime,'clockLow',s.scale.physicalTimeLow);
end
