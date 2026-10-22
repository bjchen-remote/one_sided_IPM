function wallCore = recordWallCore(wallCore,state,historyIndex)
%IPM.OUTPUT.RECORDWALLCORE Save a compact moving wall window.
%   The record contains at most 129 wall nodes centered on the tracked
%   singular feature.  It is small enough to keep in every checkpoint and
%   avoids storing the full two-dimensional state at every output time.

if nargin < 1 || isempty(wallCore) || ~isstruct(wallCore) || ...
        ~isscalar(wallCore) || isempty(fieldnames(wallCore))
    wallCore = ipm.output.emptyWallCoreHistory();
end

rho = state.rho;
flow = state.flow;
ops = state.ops;
scale = state.scale;
Cx = exp(scale.logC_l);
Comega = exp(scale.logC_omega);

trackedPeakX = flow.trackedPeakX;
if ~isfinite(trackedPeakX)
    [~,peakIndex] = max(abs(flow.source(1,:)));
    trackedPeakX = ops.x(peakIndex);
else
    [~,peakIndex] = min(abs(ops.x-trackedPeakX));
end

halfNodeCount = 64;
firstIndex = max(1,peakIndex-halfNodeCount);
lastIndex = min(ops.nx,peakIndex+halfNodeCount);
index = firstIndex:lastIndex;

physicalX = (ops.x(index)-scale.X_shift)/Cx;
physicalTrackedPeakX = (trackedPeakX-scale.X_shift)/Cx;
physicalRho = rho(1,index)/Comega;
physicalRhoX = (Cx/Comega)*flow.source(1,index);
physicalCoreHalfWidth = flow.trackedWallCoreWidth/Cx;
if ~isfinite(physicalCoreHalfWidth) || physicalCoreHalfWidth <= 0
    physicalCoreHalfWidth = max(abs(physicalX-physicalTrackedPeakX));
end

wallCore.historyIndex(end+1,1) = historyIndex;
wallCore.acceptedStep(end+1,1) = state.step;
wallCore.canonicalTime(end+1,1) = scale.canonicalTime;
wallCore.physicalTime(end+1,1) = scale.physicalTime;
wallCore.trackedPeakX(end+1,1) = trackedPeakX;
wallCore.physicalTrackedPeakX(end+1,1) = physicalTrackedPeakX;
wallCore.physicalCoreHalfWidth(end+1,1) = physicalCoreHalfWidth;
wallCore.physicalRhoXMaximum(end+1,1) = ...
    ipm.diagnostics.physicalRhoXInf(flow,scale);
wallCore.gridIndex{end+1,1} = index;
wallCore.rescaledX{end+1,1} = ops.x(index);
wallCore.physicalX{end+1,1} = physicalX;
wallCore.comovingPhysicalX{end+1,1} = ...
    physicalX-physicalTrackedPeakX;
wallCore.physicalRho{end+1,1} = physicalRho;
wallCore.physicalRhoX{end+1,1} = physicalRhoX;
end
