function family=referenceAxisFamily(rootX,rootY,anchor,policy)
%IPM.REMESH.REFERENCEAXISFAMILY Immutable, registered index-PCHIP base axes.
% Only geometry is constructed. Old knots, box endpoints, 0 and +/-anchor
% are exact. Member indices are registration indices, not execution ranks.
% A resource cap is a node budget, not an estimate of available LU memory.
identifier='ipm:AutonomousMeshReferenceFamily';
assert(isstruct(policy)&&isscalar(policy)&&isfield(policy,'version')&& ...
    any(policy.version == [2,3,4,5]),identifier,'A resolved variable-node policy is required.');
normalized=ipm.config.autonomousMeshPolicy(policy);
assert(isequaln(normalized,policy),identifier,'Pass the normalized registered policy.');
validateattributes(rootX,{'double'},{'row','real','finite','increasing','>=',-realmax});
validateattributes(rootY,{'double'},{'column','real','finite','increasing','nonnegative'});
validateattributes(anchor,{'double'},{'scalar','real','finite','positive'});
assert(numel(rootX)>=7&&mod(numel(rootX),2)==1&&numel(rootY)>=7&& ...
    isequal(rootX,-fliplr(rootX))&&rootY(1)==0&&rootY(end)>0&& ...
    anchor<rootX(end)&&any(rootX==anchor)&&any(rootX==-anchor), ...
    identifier,'Root axes must have exact symmetry, origin, interior anchors and at least seven nodes.');
registration=policy.nodeFamily;
assert(numel(rootX)*numel(rootY)<=registration.maximumTotalNodes, ...
    'ipm:AutonomousMeshResourceCap','The root member must fit maximumTotalNodes.');

positive=rootX(rootX>=0);
refinedPositive=pchip(linspace(0,1,numel(positive)),positive, ...
    linspace(0,1,2*numel(positive)-1));
refinedPositive(1:2:end)=positive;
refinedPositive(1)=0;
xAxes={rootX,[-fliplr(refinedPositive(2:end)),refinedPositive]};
refinedY=pchip(linspace(0,1,numel(rootY)),rootY, ...
    linspace(0,1,2*numel(rootY)-1))';
refinedY(1:2:end)=rootY;
refinedY(1)=0;
yAxes={rootY,refinedY};
xQuality=cell(1,2);yQuality=cell(1,2);
for factor=1:2
    x=xAxes{factor};y=yAxes{factor};
    assert(isequal(x(1:factor:end),rootX)&&isequal(y(1:factor:end),rootY)&& ...
        isequal(x,-fliplr(x))&&any(x==anchor)&&any(x==-anchor)&& ...
        isequal(x([1,end]),rootX([1,end]))&&isequal(y([1,end]),rootY([1,end])), ...
        identifier,'Index refinement must preserve every root knot and endpoint exactly.');
    xQuality{factor}=checked_quality(x,policy.qualityLimits,'x',factor);
    yQuality{factor}=checked_quality(y,policy.qualityLimits,'y',factor);
end
% Version three adds directional factor three, independently from the immutable
% root. Root knots are exact; factor-two intermediate knots need not nest.
if any(policy.version == [3,4,5])
    refinedThird=pchip(linspace(0,1,numel(positive)),positive, ...
        linspace(0,1,3*(numel(positive)-1)+1));
    refinedThird(1:3:end)=positive;
    refinedThird(1)=0;
    xAxes{3}=[-fliplr(refinedThird(2:end)),refinedThird];
    x=xAxes{3};
    assert(isequal(x(1:3:end),rootX)&&isequal(x,-fliplr(x))&& ...
        any(x==anchor)&&any(x==-anchor)&&isequal(x([1,end]),rootX([1,end])), ...
        identifier,'Factor three must preserve original root knots and anchors exactly.');
    xQuality{3}=checked_quality(x,policy.qualityLimits,'x',3);
    yAxes{3}=pchip(linspace(0,1,numel(rootY)),rootY, ...
        linspace(0,1,3*(numel(rootY)-1)+1))';
    yAxes{3}(1:3:end)=rootY;yAxes{3}(1)=0;
    assert(isequal(yAxes{3}(1:3:end),rootY)&& ...
        isequal(yAxes{3}([1,end]),rootY([1,end])), ...
        identifier,'Factor three must preserve original y root knots and endpoints exactly.');
    yQuality{3}=checked_quality(yAxes{3},policy.qualityLimits,'y',3);
end
if policy.version==5
    xFactors=unique(registration.cellFactors(:,1))';
    yFactors=unique(registration.cellFactors(:,2))';
    for factor=xFactors(xFactors>3)
        refined=pchip(linspace(0,1,numel(positive)),positive, ...
            linspace(0,1,factor*(numel(positive)-1)+1));
        refined(1:factor:end)=positive;refined(1)=0;
        xAxes{factor}=[-fliplr(refined(2:end)),refined];
        x=xAxes{factor};
        assert(isequal(x(1:factor:end),rootX) && isequal(x,-fliplr(x)) && ...
            any(x==anchor)&&any(x==-anchor) && ...
            isequal(x([1,end]),rootX([1,end])),identifier, ...
            'Automatic X factor must preserve root knots, symmetry, and anchor.');
        xQuality{factor}=checked_quality(x,policy.qualityLimits,'x',factor);
    end
    for factor=yFactors(yFactors>3)
        refined=pchip(linspace(0,1,numel(rootY)),rootY, ...
            linspace(0,1,factor*(numel(rootY)-1)+1))';
        refined(1:factor:end)=rootY;refined(1)=0;
        yAxes{factor}=refined;
        assert(isequal(refined(1:factor:end),rootY)&& ...
            isequal(refined([1,end]),rootY([1,end])),identifier, ...
            'Automatic Y factor must preserve root knots and endpoints.');
        yQuality{factor}=checked_quality(refined,policy.qualityLimits,'y',factor);
    end
end
members=struct('index',{},'cellFactors',{},'nodeCount',{},'baseX',{},'baseY',{}, ...
    'resourceAdmitted',{},'qualityPassed',{},'xQuality',{},'yQuality',{});
for index=1:size(registration.cellFactors,1)
    factors=registration.cellFactors(index,:);
    x=xAxes{factors(1)};y=yAxes{factors(2)};
    members(index)=struct('index',index,'cellFactors',factors, ...
        'nodeCount',[numel(x),numel(y)],'baseX',x,'baseY',y, ...
        'resourceAdmitted',numel(x)*numel(y)<=registration.maximumTotalNodes, ...
        'qualityPassed',true,'xQuality',xQuality{factors(1)},'yQuality',yQuality{factors(2)});
end
family=struct('version',1,'generator',registration.generator,'rootX',rootX, ...
    'rootY',rootY,'anchor',anchor,'members',members);
if any(policy.version == [3,4]),family.version=2;end
if policy.version==5,family.version=3;end
end

function q=checked_quality(axis,limits,label,factor)
q=ipm.mesh.quality(axis,7,ipm.mesh.quadrature(axis));
values=[q.maximumAdjacentCellRatio,q.maximumLogSpacingCurvature,q.minimumStencilRcond, ...
    q.minimumQuadratureWeightRatio,q.minimumQuadratureWeightToControlWidthRatio, ...
    q.maximumQuadratureWeightToControlWidthRatio];
passed=all(isfinite(values))&&q.quadratureWeightsStrictlyPositive&& ...
    values(1)<=limits.maxAdjacentCellRatio&&values(2)<=limits.maxLogSpacingCurvature&& ...
    values(3)>=limits.minStencilRcond&&values(4)>=limits.minQuadratureWeightRatio&& ...
    values(5)>=limits.minWeightToControlWidth&&values(6)<=limits.maxWeightToControlWidth;
assert(passed,'ipm:AutonomousMeshReferenceQuality', ...
    'Registered %s reference factor %d fails original axis quality; no family may be admitted.',label,factor);
end
