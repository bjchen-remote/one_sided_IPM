function measurement = measureFeatures(values,axis,levels,mask)
%IPM.DIAGNOSTICS.MEASUREFEATURES Measure one nonnegative peak profile without policy.

feature = max(values(:)',0);
if nargin >= 4 && ~isempty(mask)
    feature(~mask) = 0;
end
measurement = ipm.diagnostics.peakResolution(feature,axis,levels);
measurement.feature = feature;
end
