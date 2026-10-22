function [dt,stopReason,details] = selectTimestep(flow,scale,ops,time)
%IPM.EVOLVE.SELECTTIMESTEP Apply CFL and canonical/physical clock limits.
%   DETAILS records the transport rate, every candidate limit, the selected
%   limiter, and the CFL actually realized by DT.

spatial=isfield(ops,'spatialOperators');
if spatial
    assert(isfield(flow,'generatorCflRate')&&isscalar(flow.generatorCflRate)&& ...
        isfinite(flow.generatorCflRate)&&flow.generatorCflRate>=0, ...
        'ipm:SpatialTransportRate','The spatial backend must supply its complete transport rate.');
    rate=flow.generatorCflRate;
elseif strcmp(ops.spatialDiscretization,'sixth_order')
    % Mapped WENO evolves on the uniform computational grid. The minimum of
    % J*dXi and the adjacent physical faces avoids optimistic CFL estimates
    % near a rapidly changing but still admitted metric.
    spacingX=ops.transportSpacingX;spacingY=ops.transportSpacingY;
elseif any(strcmp(ops.transportScheme,{'weno5_fd','weno7_fd'}))
    % The mapped WENO operator uses (1/J)*d/dXi. After remeshing, J*dXi
    % can be smaller than the physical centered gaps ops.hx/ops.hy.
    spacingX=ops.metricX*ops.computationalSpacingX;
    spacingY=ops.metricY*ops.computationalSpacingY;
else
    spacingX=ops.hx;spacingY=ops.hy;
end
if ~spatial
    rate = max(abs(flow.transportU1)./spacingX + ...
        abs(flow.transportU2)./spacingY,[],'all');
end
% Optional experimental CFL mask; it is off by default.
% WENO5-FD uses up to a three-node halo in each direction. The source
% threshold is only a timestep experiment: every stage still evaluates the
% complete, unmodified semi-discrete RHS on the original grid.
generatorSplit=ops.rescaling.enabled && ...
    any(strcmp(ops.rescaling.scalingContract, ...
    {'wall_integral_generators_v2','wall_integral_centered_generators_v3', ...
    'wall_integral_centroid_practical_v1'}));
maskSetting = getenv('IPM_SOURCE_CFL_REL');
if spatial&&~isempty(maskSetting)
    error('ipm:GeneratorCflMask','The spatial backend requires its complete physical transport rate.');
end
if generatorSplit
    if ~isempty(maskSetting)
        error('ipm:GeneratorCflMask', ...
            'The fixed-generator contract requires the full absolute generator CFL.');
    end
    if ~isfield(flow,'generatorCflRate') || ...
            ~isnumeric(flow.generatorCflRate) || ~isscalar(flow.generatorCflRate) || ...
            ~isreal(flow.generatorCflRate) || ~isfinite(flow.generatorCflRate) || ...
            flow.generatorCflRate<0
        error('ipm:GeneratorCflRate','The complete generator CFL rate is missing or invalid.');
    end
    rate=flow.generatorCflRate;
end
if ~isempty(maskSetting)
    relativeThreshold = str2double(maskSetting);
    if ~(isfinite(relativeThreshold) && relativeThreshold > 0 && ...
            relativeThreshold < 1)
        error('ipm:ExperimentalCflMask', ...
            'IPM_SOURCE_CFL_REL must lie strictly between zero and one.');
    end
    if ~isfield(flow,'source') || ...
            ~isequal(size(flow.source),[ops.ny,ops.nx])
        error('ipm:ExperimentalCflMask', ...
            'The experimental CFL mask requires the complete rho_x field.');
    end
    rateField = abs(flow.transportU1)./spacingX + ...
        abs(flow.transportU2)./spacingY;
    sourceMagnitude = abs(flow.source);
    sourceMaximum = max(sourceMagnitude,[],'all');
    if sourceMaximum > 0
        active = sourceMagnitude >= relativeThreshold*sourceMaximum;
        active = conv2(double(active),ones(7,7),'same') > 0;
        rate = max(rateField(active),[],'all');
    end
end
if rate > 0
    dtCfl = time.cfl/rate;
else
    dtCfl = time.maxDt;
end

physicalClockSpeed = exp(scale.logC_omega-scale.logC_l);
implicit=any(strcmp(time.timeIntegrator,{'sdirk3_alexander_practical_v1', ...
    'sdirk3_alexander_practical_v2','sdirk3_alexander_practical_v3','sdirk2_alexander_practical_v1','sdirk2_alexander_practical_v2'}));
stabilityLimit=dtCfl;
if implicit,stabilityLimit=Inf;end
limits = [time.maxDt,stabilityLimit, ...
    time.finalTime-scale.canonicalTime, ...
    (time.physicalFinalTime-scale.physicalTime)/physicalClockSpeed];
[dt,activeIndex] = min(limits);
activeLimiters = {'max_dt','cfl','canonical_final_time', ...
    'physical_final_time'};
details = struct('dt',dt,'rate',rate,'dtCfl',dtCfl, ...
    'maxDtLimit',limits(1),'canonicalLimit',limits(3), ...
    'physicalLimit',limits(4),'physicalClockSpeed',physicalClockSpeed, ...
    'realizedCfl',dt*rate, ...
    'activeLimiter',activeLimiters{activeIndex});
if implicit
    details.stabilityMode=['implicit_',time.timeIntegrator];
    details.explicitCflLimitApplied=false;
    details.fullTransportRateUnmodified=true;
    if any(strcmp(time.timeIntegrator,{'sdirk3_alexander_practical_v2','sdirk3_alexander_practical_v3','sdirk2_alexander_practical_v1','sdirk2_alexander_practical_v2'}))
        details.initialCfl=details.realizedCfl;
    end
end

stopReason = '';
if dt >= time.minDt
    return;
end
if time.physicalFinalTime-scale.physicalTime <= ...
        time.minDt*physicalClockSpeed
    stopReason = 'physical_final_time';
elseif time.finalTime-scale.canonicalTime <= time.minDt
    stopReason = 'final_time';
else
    stopReason = 'time_step_below_minimum';
end
end
