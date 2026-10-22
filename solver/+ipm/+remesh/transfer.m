function [rhoNew,opsNew,metrics] = ...
    transfer(rho,ops,config,proposal)
%IPM.REMESH.TRANSFER Transfer a field and rebuild its numerical operators.

xNew = proposal.x;
yNew = proposal.y;
registeredTarget = [];
if isfield(config.remesh,'autonomousMesh') && ...
        config.remesh.autonomousMesh.enabled && ...
        any(config.remesh.autonomousMesh.version == [2,3,4,5])
    % The frozen config still describes the original node counts.  Only a
    % reproducible, explicitly registered member may override runtime N.
    registeredTarget = registered_target(rho,ops,config,proposal);
end
sixthOrderTransfer = ...
    strcmp(config.remesh.remeshTransferScheme,'sixth_order');
highOrderTransfer = strcmp(config.remesh.remeshTransferScheme,'high_order');
if sixthOrderTransfer
    % Build and validate the candidate first, then conserve against the
    % exact one-dimensional factors of the solver's old and new tensor
    % norms.  Analytic sinh grids use an exact metric, whereas remeshed
    % custom grids use the paired discrete metric; recomputing the old norm
    % from coordinates would therefore target a subtly different mass.
    gridOverride = config.grid;
    gridOverride.customX = xNew;
    gridOverride.customY = yNew;
    if ~isempty(registeredTarget)
        gridOverride.nx = registeredTarget.nodeCount(1);
        gridOverride.ny = registeredTarget.nodeCount(2);
    end
    opsNew = ipm.mesh.build(config,gridOverride);
    rhoX = ipm.remesh.interpolate(ops.x,rho,xNew, ...
        struct('sampleDimension',2,'conservation','none', ...
        'stencilWidth',8));
    rhoX = conserve_rows_smooth( ...
        rho,rhoX,ops.integrationHx,opsNew.integrationHx);
elseif highOrderTransfer
    if ipm.mesh.isQuadrant(ops)
        gridOverride = config.grid;
        gridOverride.customX = xNew;
        gridOverride.customY = yNew;
        opsNew = ipm.mesh.build(config,gridOverride);
        rhoX = ipm.remesh.interpolate(ops.x,rho,xNew, ...
            struct('sampleDimension',2,'conservation','none'));
        rhoX = conserve_rows_smooth( ...
            rho,rhoX,ops.integrationHx,opsNew.integrationHx);
    else
        rhoX = ipm.remesh.interpolate(ops.x,rho,xNew, ...
            struct('sampleDimension',2,'conservation','constant'));
    end
else
    rhoX = interp1(ops.x,rho',xNew,'pchip')';
    if ~strcmp(ops.rescalingMode,'physical')
        rhoX = conserve_rows(rho,rhoX,ops.hx,xNew);
    end
end
if strcmp(ops.symmetryMode,'double_odd_omega') && ~ipm.mesh.isQuadrant(ops)
    rhoX = 0.5*(rhoX+fliplr(rhoX));
end
if sixthOrderTransfer
    rhoNew = ipm.remesh.interpolate(ops.y,rhoX,yNew, ...
        struct('sampleDimension',1,'conservation','none', ...
        'stencilWidth',8));
    rhoNew = conserve_columns_with_boundary_trace( ...
        rhoX,rhoNew,ops.y,yNew,8, ...
        ops.integrationHy,opsNew.integrationHy);
elseif highOrderTransfer
    rhoNew = ipm.remesh.interpolate(ops.y,rhoX,yNew, ...
        struct('sampleDimension',1,'conservation','none'));
    rhoNew = conserve_columns_with_boundary_trace( ...
        rhoX,rhoNew,ops.y,yNew);
else
    rhoNew = interp1(ops.y,rhoX,yNew,'pchip');
end
% The PCHIP path is range preserving.  The high-order path preserves both
% y-boundary traces during its smooth wall-normal conservation correction.
if strcmp(ops.symmetryMode,'double_odd_omega') && ~ipm.mesh.isQuadrant(ops)
    rhoNew = 0.5*(rhoNew+fliplr(rhoNew));
end
if ~sixthOrderTransfer && ~(highOrderTransfer && ipm.mesh.isQuadrant(ops))
    gridOverride = config.grid;
    gridOverride.customX = xNew;
    gridOverride.customY = yNew;
    if ~isempty(registeredTarget)
        gridOverride.nx = registeredTarget.nodeCount(1);
        gridOverride.ny = registeredTarget.nodeCount(2);
    end
    opsNew = ipm.mesh.build(config,gridOverride);
end
rescaling = ops.rescaling;
[~,rescaling.originIndex] = min(abs(xNew));
[~,rescaling.pinIndex] = min(abs(xNew-rescaling.pinX));
opsNew.rescaling = rescaling;
opsNew.baseX = ops.baseX;
opsNew.baseY = ops.baseY;
if ~isempty(registeredTarget)
    opsNew.baseX = registeredTarget.baseX;
    opsNew.baseY = registeredTarget.baseY;
end
opsNew.remeshCount = ops.remeshCount;

newOmega = rhoNew*opsNew.Dx';
xMetrics = ipm.diagnostics.peakResolution(max(newOmega(1,:),0),opsNew.x, ...
    opsNew.rescaling.adaptiveLevels);
xCore = xMetrics.gridPoints(end);
peak = xMetrics.peak;
yMetrics = ipm.diagnostics.peakResolution(abs(newOmega(:,xMetrics.peakIndex))', ...
    opsNew.y',opsNew.rescaling.adaptiveLevels);
yCore = yMetrics.gridPoints(end);
oldMinimum = min(rho,[],'all');
oldMaximum = max(rho,[],'all');
rangeScale = max([oldMaximum-oldMinimum,abs(oldMinimum), ...
    abs(oldMaximum),eps]);
rangeViolation = max([oldMinimum-min(rhoNew,[],'all'), ...
    max(rhoNew,[],'all')-oldMaximum,0]);
metrics = struct('xCore',xCore,'yCore',yCore,'peak',peak, ...
    'rangeViolation',rangeViolation, ...
    'relativeRangeViolation',rangeViolation/rangeScale);
end

function target = registered_target(rho,ops,config,proposal)
% Pure validation precedes interpolation or operator construction.  The
% caller owns the accepted ledger and budget; the candidate audit validates
% those before commit.  No duplicate reference-axis arguments are trusted.
identifier = 'ipm:AutonomousMeshFamilyTransfer';
required = {'referenceFamily','sourceLevelId','targetLevelId'};
assert(all(isfield(proposal,required)),identifier, ...
    'A version-2 transfer must identify its registered source and target levels.');
policy = config.remesh.autonomousMesh;
family = proposal.referenceFamily;
assert(isstruct(family) && isscalar(family) && ...
    all(isfield(family,{'rootX','rootY','anchor','members'})),identifier, ...
    'The complete immutable reference family is required.');
assert(strcmp(config.transport.spatialDiscretization,'high_order') && ...
    strcmp(config.remesh.remeshTransferScheme,'high_order'),identifier, ...
    'Version-2 node-family transfers use the maintained high-order path.');
assert(numel(family.rootX)==config.grid.nx && ...
    numel(family.rootY)==config.grid.ny && ...
    family.anchor==config.scaling.transportAnchorX && ...
    isequal(family.rootX([1,end]),config.grid.xlim) && ...
    isequal(family.rootY([1,end]),[0;config.grid.ymax]),identifier, ...
    'Family roots and anchor must match the frozen original configuration.');
expected = ipm.remesh.referenceAxisFamily( ...
    family.rootX,family.rootY,config.scaling.transportAnchorX,policy);
assert(isequaln(family,expected),identifier, ...
    'Every registered family member must exactly reproduce from its root axes.');
assert(isnumeric(proposal.sourceLevelId) && isscalar(proposal.sourceLevelId) && ...
    isnumeric(proposal.targetLevelId) && isscalar(proposal.targetLevelId),identifier, ...
    'Source and target level identifiers must be numeric scalars.');
ids = [proposal.sourceLevelId,proposal.targetLevelId];
assert(isreal(ids) && ...
    all(isfinite(ids)) && all(ids==fix(ids)) && ...
    all(ids>=1) && all(ids<=numel(family.members)),identifier, ...
    'Source and target must name registered levels.');
source = family.members(ids(1));target = family.members(ids(2));
assert(isequal([numel(ops.x),numel(ops.y)],source.nodeCount) && ...
    isequal(ops.baseX,source.baseX) && isequal(ops.baseY,source.baseY) && ...
    isequal(ops.x([1,end]),source.baseX([1,end])) && ...
    isequal(ops.y([1,end]),source.baseY([1,end])),identifier, ...
    'The actual source axes and base references must match the source member.');
assert(isequal(size(rho),[source.nodeCount(2),source.nodeCount(1)]) && ...
    isnumeric(rho) && isreal(rho) && all(isfinite(rho),'all'),identifier, ...
    'The source field must be finite and paired with its registered axes.');
assert(source.resourceAdmitted && source.qualityPassed && ...
    target.resourceAdmitted && target.qualityPassed && ...
    all(target.cellFactors>=source.cellFactors),identifier, ...
    'A transfer may only stay in or grow monotonically to an admitted member.');
assert(isrow(proposal.x) && iscolumn(proposal.y) && ...
    isreal(proposal.x) && isreal(proposal.y) && ...
    all(isfinite(proposal.x)) && all(isfinite(proposal.y)) && ...
    all(diff(proposal.x)>0) && all(diff(proposal.y)>0) && ...
    isequal([numel(proposal.x),numel(proposal.y)],target.nodeCount) && ...
    isequal(proposal.x([1,end]),ops.x([1,end])) && ...
    isequal(proposal.y([1,end]),ops.y([1,end])) && ...
    isequal(proposal.x,-fliplr(proposal.x)) && proposal.y(1)==0 && ...
    any(proposal.x==family.anchor) && any(proposal.x==-family.anchor),identifier, ...
    'Target coordinates must match registered N, preserve the box and exact anchors.');
end

function corrected = conserve_columns_with_boundary_trace( ...
        oldField,newField,oldAxis,newAxis,stencilWidth, ...
        oldWeights,newWeights)
% A smooth bubble correction preserves each wall-normal column integral
% without changing either boundary trace.  Its amplitude has the same
% high-order size as the interpolation/quadrature defect.
if nargin < 5
    oldWeights = ipm.mesh.quadrature(oldAxis);
    newWeights = ipm.mesh.quadrature(newAxis);
elseif stencilWidth == 8
    if nargin < 7
        oldWeights = ipm.mesh.quadrature6(oldAxis);
        newWeights = ipm.mesh.quadrature6(newAxis);
    end
else
    error('ipm:HighOrderRemeshStencil', ...
        'The remesh correction stencil width must be 8 when specified.');
end
oldMass = oldWeights(:)'*oldField;
newMass = newWeights(:)'*newField;
bubble = (newAxis(:)-newAxis(1)).*(newAxis(end)-newAxis(:));
bubbleMass = newWeights(:)'*bubble;
if bubbleMass <= 0
    error('ipm:HighOrderRemeshBubble', ...
        'The wall-normal conservative correction has zero capacity.');
end
correction = (oldMass-newMass)/bubbleMass;
corrected = newField+bubble*correction;
% Remove the final dot-product roundoff while retaining zero boundary data.
residual = oldMass-newWeights(:)'*corrected;
corrected = corrected+bubble*(residual/bubbleMass);
end

function corrected = conserve_rows_smooth( ...
        oldField,newField,oldWeights,newWeights)
% A constant-in-x correction preserves every row integral in the exact
% tensor factors supplied by the accepted and candidate operators.
oldWeights = oldWeights(:);
newWeights = newWeights(:);
oldMass = oldField*oldWeights;
newMass = newField*newWeights;
weightSum = sum(newWeights);
if weightSum <= 0
    error('ipm:SixthOrderQuadratureWeights', ...
        'The sixth-order row-conservation norm must remain positive.');
end
correction = (oldMass-newMass)/weightSum;
corrected = newField+correction;
residual = oldMass-corrected*newWeights;
corrected = corrected+residual/weightSum;
end

function corrected = conserve_rows(oldField,newField,oldWeights,newAxis)
% Match each slice integral while keeping every value in its old range.
newWeights = zeros(size(newAxis));
newWeights([1,end]) = [newAxis(2)-newAxis(1), ...
    newAxis(end)-newAxis(end-1)]/2;
newWeights(2:end-1) = (newAxis(3:end)-newAxis(1:end-2))/2;
corrected = newField;
for row = 1:size(oldField,1)
    oldMass = sum(oldField(row,:).*oldWeights);
    newMass = sum(corrected(row,:).*newWeights);
    defect = oldMass-newMass;
    rowMinimum = min(oldField(row,:));
    rowMaximum = max(oldField(row,:));
    rowScale = max(rowMaximum-rowMinimum,1);
    tolerance = 100*eps(rowScale);
    if abs(defect) <= tolerance*sum(newWeights)
        continue;
    end
    if defect > 0
        capacity = rowMaximum-corrected(row,:);
        transition = corrected(row,:) > rowMinimum+tolerance;
    else
        capacity = corrected(row,:)-rowMinimum;
        transition = corrected(row,:) < rowMaximum-tolerance;
    end
    capacity(~transition) = 0;
    capacityMass = sum(capacity.*newWeights);
    if capacityMass < abs(defect)
        if defect > 0
            capacity = rowMaximum-corrected(row,:);
        else
            capacity = corrected(row,:)-rowMinimum;
        end
        capacityMass = sum(capacity.*newWeights);
    end
    if capacityMass < abs(defect)-tolerance*sum(newWeights)
        error('ipm:ConservativeRemesh', ...
            'Insufficient bounded capacity for conservative remapping.');
    end
    fraction = min(abs(defect)/max(capacityMass,eps),1);
    if defect > 0
        corrected(row,:) = corrected(row,:)+fraction*capacity;
    else
        corrected(row,:) = corrected(row,:)-fraction*capacity;
    end
end
end
