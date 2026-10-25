function [values,diagnostic]=compare(state,g,method,q)
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

values=ipm_cr.diag.fields(state,g,method,q.x,q.y);f=state.rhs;
diagnostic=struct('rawFullArrayReinterpolationDifference',struct());
for name={'W','q','Q'}
    key=name{1};field=state.W;if ~strcmp(key,'W'),field=f.(key);end
    raw=ipm_cr.mesh.sample(field,g,q.x,q.y);
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
