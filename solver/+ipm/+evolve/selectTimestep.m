function [dt,stopReason,details] = selectTimestep(flow,scale,ops,time)
%IPM.EVOLVE.SELECTTIMESTEP Apply CFL and canonical/physical clock limits.
%   DETAILS records the transport rate, every candidate limit, the selected
%   limiter, and the CFL actually realized by DT.

if strcmp(ops.spatialDiscretization,'sixth_order')
    % Mapped WENO evolves on the uniform computational grid. The minimum of
    % J*dXi and the adjacent physical faces avoids optimistic CFL estimates
    % near a rapidly changing but still admitted metric.
    rate = max(abs(flow.transportU1)./ops.transportSpacingX + ...
        abs(flow.transportU2)./ops.transportSpacingY,[],'all');
else
    rate = max(abs(flow.transportU1)./ops.hx + ...
        abs(flow.transportU2)./ops.hy,[],'all');
end
% Isolated CFL-mask audit only. The production package has no such branch.
% WENO5-FD uses up to a three-node halo in each direction. The source
% threshold is only a timestep experiment: every stage still evaluates the
% complete, unmodified semi-discrete RHS on the original grid.
maskSetting = getenv('IPM_SOURCE_CFL_REL');
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
    if strcmp(ops.spatialDiscretization,'sixth_order')
        rateField = abs(flow.transportU1)./ops.transportSpacingX + ...
            abs(flow.transportU2)./ops.transportSpacingY;
    else
        rateField = abs(flow.transportU1)./ops.hx + ...
            abs(flow.transportU2)./ops.hy;
    end
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
limits = [time.maxDt,dtCfl, ...
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
