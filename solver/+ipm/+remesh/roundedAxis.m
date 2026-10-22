function [axis,info]=roundedAxis( ...
    anchor,halfWidth,nodeCount,fineSpacing,fineCells, ...
    anchorFraction,roundingCells,positiveOnly)
%IPM.REMESH.ROUNDEDAXIS Candidate-only rounded-log axis at arbitrary units.
% nodeCount counts stored nodes: [0,H] for positiveOnly, [-H,H] otherwise.
% The rounded bridge/geometric laws and split objective are from the maintained
% fourth_order_smooth_peak_grid coordinated constructor. This independent
% direct family removes its q256/power-of-two restriction. It is not a claim
% to reproduce a historical cropped grid or the factory's nested q512 axis.
if nargin<8
    positiveOnly=false;
end
validateattributes(anchor,{'numeric'},{'scalar','finite','positive'});
validateattributes(halfWidth,{'numeric'},{'scalar','finite','>',anchor});
if positiveOnly
    validateattributes(nodeCount,{'numeric'},{'scalar','integer','>=',9});
    cells=nodeCount-1;
else
    validateattributes(nodeCount,{'numeric'},{'scalar','integer','>=',17});
    assert(mod(nodeCount,2)==1);
    cells=(nodeCount-1)/2;
end
validateattributes(fineSpacing,{'numeric'},{'scalar','finite','positive'});
validateattributes(fineCells,{'numeric'},{'scalar','integer','>=',8});
validateattributes(roundingCells,{'numeric'},{'scalar','integer','>=',2});
validateattributes(anchorFraction,{'numeric'},{'scalar','>',0,'<',1});
H=halfWidth/anchor;
leftFineCells=round(fineCells*anchorFraction);
assert(leftFineCells>=1 && leftFineCells<fineCells);
rightFineCells=fineCells-leftFineCells;
requestedSpacing=fineSpacing/anchor;
designSpacing=requestedSpacing;
spacingCap=Inf;
spacingLimited=false;
if positiveOnly
    minimumBranch=roundingCells+2;
    outerCells=cells-fineCells;
    if outerCells<2*minimumBranch
        error('ipm:RoundedAxisInfeasible', ...
            'Insufficient rounded-branch cell budget.');
    end
    % Largest constant fine width that leaves a feasible integer split.
    % The continuous balance is exact; only its two adjacent integers can
    % maximize the minimum left/right capacity.
    idealLeft=cells/H-leftFineCells;
    candidates=unique(min(outerCells-minimumBranch, ...
        max(minimumBranch,[floor(idealLeft),ceil(idealLeft)])));
    candidateCaps=min(1./(candidates+leftFineCells), ...
        (H-1)./(outerCells-candidates+rightFineCells));
    spacingCap=max(candidateCaps);
    safeCap=spacingCap-128*eps(spacingCap);
    designSpacing=min(designSpacing,safeCap);
    spacingLimited=designSpacing<requestedSpacing;
else
    assert(fineCells+2*(roundingCells+2)<=cells, ...
        'Insufficient rounded-branch cell budget.');
end
interval=1+[-leftFineCells,rightFineCells]*designSpacing;
assert(interval(1)>0 && interval(2)<H);
% Preserve the established interval-to-width floating-point operation.
spacing=diff(interval)/fineCells;
leftLength=interval(1);
rightLength=H-interval(2);
best=struct();
bestObjective=Inf;
bestMismatch=Inf;
if positiveOnly
    % The level-set path has one geometry proposal. Allocate its outer-cell
    % budget from the two measured spans, then solve only two scalar slopes.
    % Projection enforces the minimum rounded branch and spacing feasibility.
    outerCells=cells-fineCells;
    lower=max(minimumBranch,outerCells-floor(rightLength/spacing));
    upper=min(outerCells-minimumBranch,floor(leftLength/spacing));
    if lower>upper
        error('ipm:RoundedAxisInfeasible', ...
            'No feasible direct rounded-log cell split.');
    end
    spanWeight=log1p([leftLength,rightLength]/spacing);
    leftCount=min(upper,max(lower,round(outerCells*spanWeight(1)/sum(spanWeight))));
    rightCount=outerCells-leftCount;
    leftSlope=growth(leftLength,spacing,leftCount,roundingCells,true);
    rightSlope=growth(rightLength,spacing,rightCount,roundingCells,false);
    best=struct('leftCells',leftCount,'rightCells',rightCount, ...
        'leftSlope',leftSlope,'rightSlope',rightSlope);
    splitStrategy='direct_log_span_projection';
    splitEvaluations=1;
else
    for leftCount=roundingCells+2:cells-fineCells-roundingCells-2
        rightCount=cells-fineCells-leftCount;
        if leftCount*spacing>leftLength || rightCount*spacing>rightLength
            continue
        end
        leftSlope=growth(leftLength,spacing,leftCount,roundingCells,true);
        rightSlope=growth(rightLength,spacing,rightCount,roundingCells,false);
        objective=max(leftSlope,rightSlope);
        mismatch=abs(leftSlope-rightSlope);
        tolerance=100*eps(max(1,objective));
        if isempty(fieldnames(best)) || objective<bestObjective-tolerance || ...
                (abs(objective-bestObjective)<=tolerance && mismatch<bestMismatch)
            bestObjective=objective;
            bestMismatch=mismatch;
            best=struct('leftCells',leftCount,'rightCells',rightCount, ...
                'leftSlope',leftSlope,'rightSlope',rightSlope);
        end
    end
    splitStrategy='minimum_worst_slope_scan';
    splitEvaluations= ...
        max(0,cells-fineCells-2*(roundingCells+2)+1);
end
assert(~isempty(fieldnames(best)),'No feasible rounded-log cell split.');
left=spacing*exp(best.leftSlope*fliplr(bridge(best.leftCells,roundingCells)));
right=spacing*exp(best.rightSlope*geometric(best.rightCells,roundingCells));
left(1)=left(1)+(leftLength-sum(left));
right(end)=right(end)+(rightLength-sum(right));
positive=[0,cumsum([left,spacing*ones(1,fineCells),right])];
positive(best.leftCells+1)=interval(1);positive(best.leftCells+fineCells+1)=interval(2);
positive(end)=H;
anchorPositiveIndex=best.leftCells+leftFineCells+1;
before=positive(anchorPositiveIndex);
assert(abs(before-1)<=512*eps(1));
positive(anchorPositiveIndex)=1;
if positiveOnly
    axis=anchor*positive;
    axis(1)=0;
    axis(end)=halfWidth;
    anchorIndex=anchorPositiveIndex;
    axis(anchorIndex)=anchor;
    assert(numel(axis)==nodeCount && all(diff(axis)>0) && ...
        axis(1)==0 && axis(end)==halfWidth);
else
    axis=anchor*[-fliplr(positive(2:end)),positive];
    axis(1)=-halfWidth;
    axis(end)=halfWidth;
    anchorIndex=cells+anchorPositiveIndex;
    axis(anchorIndex)=anchor;
    axis(nodeCount+1-anchorIndex)=-anchor;
    assert(isequal(axis,-fliplr(axis)) && all(diff(axis)>0) && axis(cells+1)==0);
end
info=struct('kind','direct_generic_rounded_log_candidate','anchor',anchor, ...
    'normalization','u=x/anchor','normalizedHalfWidth',H,'nodeCount',nodeCount, ...
    'fineCells',fineCells,'leftFineCells',leftFineCells,'roundingCells',roundingCells, ...
    'normalizedFineInterval',interval,'normalizedFineSpacing',spacing, ...
    'requestedNormalizedFineSpacing',requestedSpacing, ...
    'normalizedFineSpacingCap',spacingCap,'fineSpacingLimited',spacingLimited, ...
    'split',best,'anchorIndex',anchorIndex,'normalizedAnchorBeforeSnap',before, ...
    'splitStrategy',splitStrategy,'splitEvaluations',splitEvaluations, ...
    'historicalAxisReproductionClaim',false,'nestedLegacyFactory',false, ...
    'qualityNotYetAdmitted',true);
if positiveOnly
    info.conceptualFullNodeCount=2*nodeCount-1;
    info.positiveOnly=true;
end
end

function slope=growth(length,spacing,count,rounding,bothEnds)
if bothEnds
    shape=bridge(count,rounding);
else
    shape=geometric(count,rounding);
end
minimum=spacing*count;
tolerance=100*eps(max(1,length));
assert(length>=minimum-tolerance);
if length<=minimum+tolerance
    slope=0;
    return
end
lower=0;
upper=1/max(1,count);
while spacing*sum(exp(upper*shape))<length
    upper=2*upper;
    assert(upper*max(shape)<=log(realmax)-4);
end
for k=1:80
    mid=.5*(lower+upper);
    % Further midpoint steps are identical after floating-point stagnation.
    if mid==lower || mid==upper
        slope=mid;
        return
    end
    if spacing*sum(exp(mid*shape))<length
        lower=mid;
    else
        upper=mid;
    end
end
slope=.5*(lower+upper);
end

function shape=geometric(count,rounding)
index=0:count-1;
shape=index-rounding/2;
mask=index<=rounding;
shape(mask)=index(mask).^2/(2*rounding);
end

function shape=bridge(count,rounding)
index=1:count-1;
increments=min([ones(size(index));(index-.5)/rounding;(count-index-.5)/rounding],[],1);
shape=[0,cumsum(increments)];
end
