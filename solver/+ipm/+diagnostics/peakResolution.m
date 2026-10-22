function metrics = peakResolution(values,x,levels)
%IPM.DIAGNOSTICS.PEAKRESOLUTION Resolution of the positive peak on a nonuniform line.
%   Only the above-level connected component containing the global positive
%   peak contributes. This prevents tails or separated scales from hiding an
%   under-resolved locked peak.

values = max(values(:)',0);
x = x(:)';
levels = levels(:)';
if numel(values) ~= numel(x)
    error('ipm:PeakResolutionSize', ...
        'The feature values and grid coordinates must have equal length.');
end
[peak,peakIndex] = max(values);
metrics = struct('peak',peak,'peakIndex',peakIndex, ...
    'gridPoints',zeros(size(levels)),'areaPoints',zeros(size(levels)), ...
    'widths',NaN(size(levels)),'bounds',NaN(numel(levels),2));
if peak <= 0
    return;
end

for levelIndex = 1:numel(levels)
    threshold = levels(levelIndex)*peak;
    leftNode = peakIndex;
    while leftNode > 1 && values(leftNode-1) >= threshold
        leftNode = leftNode-1;
    end
    rightNode = peakIndex;
    while rightNode < numel(values) && values(rightNode+1) >= threshold
        rightNode = rightNode+1;
    end

    fractions = zeros(1,numel(values)-1);
    normalizedArea = zeros(size(fractions));
    if rightNode > leftNode
        faces = leftNode:rightNode-1;
        fractions(faces) = 1;
        normalizedArea(faces) = ...
            0.5*(values(faces)+values(faces+1))/peak;
    end
    if leftNode > 1
        face = leftNode-1;
        fraction = (values(leftNode)-threshold) / ...
            (values(leftNode)-values(face));
        fractions(face) = fraction;
        normalizedArea(face) = fraction*0.5* ...
            (values(leftNode)+threshold)/peak;
    end
    if rightNode < numel(values)
        face = rightNode;
        fraction = (values(rightNode)-threshold) / ...
            (values(rightNode)-values(rightNode+1));
        fractions(face) = fraction;
        normalizedArea(face) = fraction*0.5* ...
            (values(rightNode)+threshold)/peak;
    end
    metrics.gridPoints(levelIndex) = sum(fractions);
    metrics.areaPoints(levelIndex) = sum(normalizedArea);
    metrics.widths(levelIndex) = sum(diff(x).*fractions);
    leftCrossing = x(leftNode);
    rightCrossing = x(rightNode);
    if leftNode > 1
        leftCrossing = leftCrossing- ...
            fractions(leftNode-1)*(x(leftNode)-x(leftNode-1));
    end
    if rightNode < numel(values)
        rightCrossing = rightCrossing+ ...
            fractions(rightNode)*(x(rightNode+1)-x(rightNode));
    end
    metrics.bounds(levelIndex,:) = [leftCrossing,rightCrossing];
end
end
