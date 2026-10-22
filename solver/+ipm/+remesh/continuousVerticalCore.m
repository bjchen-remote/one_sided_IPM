function observation=continuousVerticalCore(flow,x,y,levels)
%IPM.REMESH.CONTINUOUSVERTICALCORE Observe the Y core at the gauge's
% continuous wall peak. PCHIP needs only the four neighboring X columns.
% This is a diagnostic; it does not change the wall gauge or the flow.
% The quadratic amplitude gauge exposes its own continuous vertex. Other
% amplitude gauges still have the independently tracked wall-gradient peak.
peakX=flow.omegaGaugeQuadraticPeakX;
if ~isfinite(peakX)
    peakX=flow.trackedPeakX;
end
x=x(:)';y=y(:)';
assert(isfinite(peakX) && peakX>=x(1) && peakX<=x(end) && ...
    numel(x)>=4 && size(flow.source,2)==numel(x) && ...
    size(flow.source,1)==numel(y), ...
    'ipm:ContinuousVerticalCoreInput','The continuous peak must lie on the current grid.');
interval=find(x<=peakX,1,'last');
interval=min(interval,numel(x)-1);
first=max(1,interval-1);last=min(numel(x),interval+2);
column=interp1(x(first:last),flow.source(:,first:last)',peakX,'pchip');
resolution=ipm.diagnostics.peakResolution(abs(column),y,levels);
assert(all(isfinite(resolution.gridPoints)) && ...
    all(isfinite(resolution.widths)) && all(resolution.gridPoints>0) && ...
    all(resolution.widths>0), ...
    'ipm:ContinuousVerticalCoreInput','Continuous Y feature is unresolved.');
observation=struct('peakX',peakX,'coreCells',resolution.gridPoints(end), ...
    'coreWidth',resolution.widths(end),'peakCells',resolution.areaPoints(1), ...
    'levelCells',resolution.gridPoints,'bracket',interval);
end
