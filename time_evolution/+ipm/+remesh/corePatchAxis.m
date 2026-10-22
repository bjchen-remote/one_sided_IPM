function [x,info]=corePatchAxis(anchor,coreCenter,halfWidth,nodeCount,spacing,fineCells,fineFraction,roundingCells,warpRadius,positiveOnly)
%IPM.REMESH.COREPATCHAXIS Core-centered candidate; exact anchor via a C2 warp.
% The anchor may lie beyond the fine patch. A compact polynomial displacement
% fixes one existing transition node; full real-unit quality gates remain
% mandatory. nodeCount counts stored nodes on the selected domain.
% This factory constructs no historical grid or solver state.
if nargin<10,positiveOnly=false;end
validateattributes(warpRadius,{'numeric'},{'scalar','finite','positive','<',.9});
[base,baseInfo]=ipm.remesh.roundedAxis(coreCenter/anchor,halfWidth/anchor, ...
    nodeCount,spacing/anchor,fineCells,fineFraction,roundingCells,positiveOnly);
if positiveOnly,positive=base;else,positive=base(base>=0);end
[~,j]=min(abs(positive-1));
center=positive(j);delta=1-center;
assert(center-warpRadius>0 && center+warpRadius<positive(end));
u=(positive-center)/warpRadius;bump=zeros(size(u));inside=abs(u)<1;
bump(inside)=(1-u(inside).^2).^3;
moved=positive+delta*bump;moved(j)=1;
assert(all(diff(moved)>0));
if positiveOnly
    x=anchor*moved;x(1)=0;x(end)=halfWidth;
    anchorIndex=j;x(anchorIndex)=anchor;
    assert(all(diff(x)>0) && any(x==anchor));
else
    x=anchor*[-fliplr(moved(2:end)),moved];x(1)=-halfWidth;x(end)=halfWidth;
    anchorIndex=(nodeCount-1)/2+j;x(anchorIndex)=anchor;x(nodeCount+1-anchorIndex)=-anchor;
    assert(all(diff(x)>0) && isequal(x,-fliplr(x)) && any(x==anchor) && any(x==-anchor));
end
fineInterval=baseInfo.normalizedFineInterval*(coreCenter/anchor);
info=struct('kind','core_centered_rounded_patch_with_compact_C2_anchor_warp', ...
    'normalization','all construction coordinates divided by actual transport anchor', ...
    'actualTransportAnchor',anchor,'coreCenter',coreCenter,'baseConstruction',baseInfo, ...
    'warpRadiusInAnchorUnits',warpRadius,'normalizedAnchorNodeBeforeWarp',center, ...
    'normalizedAnchorDisplacement',delta,'anchorIndex',anchorIndex, ...
    'normalizedFineIntervalBeforeWarp',fineInterval, ...
    'anchorOutsideFineIntervalBeforeWarp',1<fineInterval(1) || 1>fineInterval(2), ...
    'maximumPhysicalDisplacement',anchor*max(abs(delta*bump)), ...
    'warpPolynomial','(1-u^2)^3 on abs(u)<1, otherwise0; value and first two derivatives vanish at endpoints', ...
    'historicalAxisReproductionClaim',false,'fullQualityNotYetAdmitted',true);
end
