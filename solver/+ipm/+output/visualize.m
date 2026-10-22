function viz = visualize( ...
        viz,~,flow,ops,scale,history,output,runMetadata)
%IPM.OUTPUT.VISUALIZE Update live diagnostics and optionally append an MP4 frame.

if output.livePlot && ~output.writeVideo && ...
        isfield(output.caseMetadata,'study') && ...
        strcmp(output.caseMetadata.study,'one_sided_main')
    viz = boundary_center_live(viz,flow,ops,scale);
    return;
end

initializing = isempty(viz);
common = history.common;
mesh = history.mesh;
gauge = history.gauge;
safetyHistory = optional_history_series( ...
    mesh,'safetyFactor',common.physicalTime);
safetyFactor = safetyHistory(end);
verticalPeakWidth = optional_history_value( ...
    mesh,'trackedVerticalPeakWidth',NaN);
horizontalCorePoints = optional_history_value( ...
    mesh,'coreGridPoints',NaN);
verticalCorePoints = optional_history_value( ...
    mesh,'verticalCoreGridPoints',NaN);
peakGridPoints = optional_history_value( ...
    mesh,'peakGridPoints',NaN);
supportGridPoints = optional_history_value( ...
    mesh,'supportGridPoints',NaN);
trackedWallCoreWidth = optional_history_value( ...
    mesh,'trackedWallCoreWidth',NaN);
if isfield(scale,'logC_y')
    Cy = exp(scale.logC_y);
else
    Cy = exp(scale.logC_l);
end
physicalOmega = exp(scale.logC_l-scale.logC_omega)*flow.source;
if ops.rescaling.enabled
    displayPeakX = flow.trackedPeakX;
else
    displayPeakX = common.wallPeakX(end);
end
physicalWallX = (ops.x-displayPeakX)/exp(scale.logC_l);
physicalY = ops.y/Cy;
physicalGridX = repmat(physicalWallX,ops.ny,1);
physicalGridY = repmat(physicalY,1,ops.nx);
if strcmp(ops.symmetryMode,'double_odd_omega')
    observationMask = ops.x >= 0;
else
    observationMask = true(size(ops.x));
end
observedX = ops.x(observationMask);
observedPhysicalWallX = physicalWallX(observationMask);

try
if initializing
    visibility = 'off';
    if output.livePlot
        visibility = 'on';
    end
    viz.figure = figure('Name','Live symmetric IPM evolution', ...
        'Color','w','Visible',visibility);
    layout = tiledlayout(viz.figure,2,2,'TileSpacing','compact');
    viz.rhoAxes = nexttile(layout);
    viz.rhoSurface = surf(viz.rhoAxes,ops.X(:,observationMask), ...
        ops.Y(:,observationMask),flow.source(:,observationMask), ...
        'EdgeColor','none');
    view(viz.rhoAxes,2); colorbar(viz.rhoAxes);
    if ops.rescaling.enabled
        xlabel(viz.rhoAxes,'rescaled X_1');
        ylabel(viz.rhoAxes,'rescaled X_2');
    else
        xlabel(viz.rhoAxes,'physical x_1');
        ylabel(viz.rhoAxes,'physical x_2');
    end

    viz.omegaAxes = nexttile(layout);
    viz.omegaSurface = surf(viz.omegaAxes, ...
        physicalGridX(:,observationMask),physicalGridY(:,observationMask), ...
        physicalOmega(:,observationMask), ...
        'EdgeColor','none');
    view(viz.omegaAxes,2); colorbar(viz.omegaAxes);
    xlabel(viz.omegaAxes,'physical x_1 relative to peak');
    ylabel(viz.omegaAxes,'physical x_2');

    viz.wallAxes = nexttile(layout);
    viz.wallLine = plot(viz.wallAxes,observedPhysicalWallX, ...
        physicalOmega(1,observationMask), ...
        'LineWidth',1.5);
    grid(viz.wallAxes,'on');
    xlabel(viz.wallAxes,'physical x_1 relative to peak');
    ylabel(viz.wallAxes,'physical d_{x_1}rho');

    viz.rateAxes = nexttile(layout);
    yyaxis(viz.rateAxes,'left');
    viz.cLLine = plot(viz.rateAxes,NaN,NaN,'LineWidth',1.4);
    hold(viz.rateAxes,'on');
    viz.cOmegaLine = plot(viz.rateAxes,NaN,NaN,'LineWidth',1.4);
    viz.cRLine = plot(viz.rateAxes,NaN,NaN,'LineWidth',1.4);
    if ops.rescaling.enabled
        ylabel(viz.rateAxes,'canonical scaling rates');
    else
        ylabel(viz.rateAxes,'physical gradient growth');
    end
    yyaxis(viz.rateAxes,'right');
    viz.innerEnergyLine = plot(viz.rateAxes,NaN,NaN,'--','LineWidth',1.0);
    hold(viz.rateAxes,'on');
    viz.nearEnergyLine = plot(viz.rateAxes,NaN,NaN,'--','LineWidth',1.0);
    viz.outerEnergyLine = plot(viz.rateAxes,NaN,NaN,'--','LineWidth',1.0);
    ylabel(viz.rateAxes,'gradient-energy fraction');
    ylim(viz.rateAxes,[0,1]);
    yyaxis(viz.rateAxes,'left');
    grid(viz.rateAxes,'on');
    if ops.rescaling.enabled
        xlabel(viz.rateAxes,'canonical tau');
        legend(viz.rateAxes,'c_l','c_omega','c_r', ...
            'E inner','E near','E outer','Location','best');
    else
        xlabel(viz.rateAxes,'physical t');
        legend(viz.rateAxes,'||grad rho||_inf','wall max rho_{x_1}', ...
            'safety','E inner','E near','E outer','Location','best');
    end

    viz.writer = [];
    if output.writeVideo
        requestedVideoFile = output.videoFile;
        if strlength(runMetadata.videoFile) > 0
            requestedVideoFile = char(runMetadata.videoFile);
        end
        videoFile = ipm.output.uniquePath( ...
            requestedVideoFile,runMetadata.caseId);
        viz.videoFile = videoFile;
        videoDirectory = fileparts(videoFile);
        if ~isempty(videoDirectory) && ~isfolder(videoDirectory)
            mkdir(videoDirectory);
        end
        viz.writer = VideoWriter(videoFile,'MPEG-4');
        viz.writer.FrameRate = output.videoFrameRate;
        viz.writer.Quality = output.videoQuality;
        open(viz.writer);
    end
end

set(viz.rhoSurface,'XData',ops.X(:,observationMask), ...
    'YData',ops.Y(:,observationMask), ...
    'ZData',flow.source(:,observationMask), ...
    'CData',flow.source(:,observationMask));
set(viz.omegaSurface,'XData',physicalGridX(:,observationMask), ...
    'YData',physicalGridY(:,observationMask), ...
    'ZData',physicalOmega(:,observationMask), ...
    'CData',physicalOmega(:,observationMask));
set(viz.wallLine,'XData',observedPhysicalWallX, ...
    'YData',physicalOmega(1,observationMask));
if ops.rescaling.enabled
    set(viz.cLLine,'XData',common.canonicalTau,'YData',common.canonicalCL);
    set(viz.cOmegaLine,'XData',common.canonicalTau, ...
        'YData',common.canonicalCOmega);
    set(viz.cRLine,'XData',common.canonicalTau,'YData',common.canonicalCR);
else
    set(viz.cLLine,'XData',common.physicalTime, ...
        'YData',common.physicalGradInf);
    set(viz.cOmegaLine,'XData',common.physicalTime, ...
        'YData',common.physicalWallOmegaPeak);
    set(viz.cRLine,'XData',common.physicalTime,'YData',safetyHistory);
end
set(viz.innerEnergyLine,'XData',common.canonicalTau, ...
    'YData',common.gradientInnerFraction);
set(viz.nearEnergyLine,'XData',common.canonicalTau, ...
    'YData',common.gradientNearFraction);
set(viz.outerEnergyLine,'XData',common.canonicalTau, ...
    'YData',common.gradientOuterFraction);

if ops.rescaling.enabled
    peakWidth = mesh.trackedWallPeakWidth(end);
    displayedPhysicalPeak = common.physicalTrackedWallOmegaPeak(end);
else
    peakWidth = common.wallPeakWidth(end);
    displayedPhysicalPeak = common.physicalWallOmegaPeak(end);
end
if ~isfinite(peakWidth) || peakWidth <= 0
    peakWidth = 8*ops.dx;
end
if ops.rescaling.enabled
    center = ops.rescaling.pinX;
else
    center = common.wallPeakX(end);
end
xRadius = max(2.5*peakWidth,12*ops.dx);
yTop = max(2.5*peakWidth,12*ops.dy);
if isfinite(verticalPeakWidth) && verticalPeakWidth > 0
    yTop = max(2.5*verticalPeakWidth,12*min(ops.dyFaces));
end
xWindow = [max(center-xRadius,observedX(1)), ...
    min(center+xRadius,observedX(end))];
yWindow = [0,yTop];
physicalPeakWidth = peakWidth/exp(scale.logC_l);
physicalXRadius = max(2.5*physicalPeakWidth, ...
    12*min(diff(physicalWallX)));
physicalYTop = max(2.5*physicalPeakWidth,12*min(diff(physicalY)));
if isfinite(verticalPeakWidth) && verticalPeakWidth > 0
    physicalYTop = max(2.5*verticalPeakWidth/Cy, ...
        12*min(diff(physicalY)));
end
xlim(viz.rhoAxes,xWindow); ylim(viz.rhoAxes,yWindow);
physicalXWindow = [max(-physicalXRadius,observedPhysicalWallX(1)), ...
    min(physicalXRadius,observedPhysicalWallX(end))];
xlim(viz.omegaAxes,physicalXWindow);
ylim(viz.omegaAxes,[0,physicalYTop]);
xlim(viz.wallAxes,physicalXWindow);

symmetryDefect = NaN;
if isfield(common,'rhoEvenDefect') && isfield(common,'omegaOddDefect')
    symmetryDefect = max(common.rhoEvenDefect(end), ...
        common.omegaOddDefect(end));
end
title(viz.rhoAxes,sprintf(['rho_{x_1} (x_1>=0): t=%.4f, x_*=%.3g, ' ...
    'sym=%.1e'],scale.physicalTime,common.physicalTrackedPeakX(end), ...
    symmetryDefect), ...
    'FontSize',9);
title(viz.omegaAxes,sprintf('x/y core=%.1f/%.1f, safe=%.2f', ...
    horizontalCorePoints,verticalCorePoints,safetyFactor), ...
    'FontSize',9);
% Keep the full three-level x counts in UserData for interactive inspection
% without forcing a long title into every video frame.
viz.omegaAxes.UserData = [ ...
    peakGridPoints,supportGridPoints,horizontalCorePoints];
title(viz.wallAxes,sprintf('peak=%.3e, FWHM=%.2e, x core width=%.2e', ...
    displayedPhysicalPeak,physicalPeakWidth, ...
    trackedWallCoreWidth/exp(scale.logC_l)), ...
    'FontSize',9);
if ops.rescaling.enabled
    widthCorrection = 0;
    widthControlMode = 0;
    if isfield(gauge,'canonicalWidthCorrection')
        widthCorrection = gauge.canonicalWidthCorrection(end);
    end
    if isfield(gauge,'widthControlMode')
        widthControlMode = gauge.widthControlMode(end);
    end
    title(viz.rateAxes,sprintf(['corr=%.2g, mode=%+d, remesh=%d, ' ...
        'min dX/dY=%.1e/%.1e'], ...
        widthCorrection,widthControlMode, ...
        mesh.remeshCount(end),mesh.minimumDx(end), ...
        mesh.minimumDy(end)),'FontSize',9);
else
    title(viz.rateAxes,sprintf(['no-c: remesh=%d, dx/dy=%.1e/%.1e, ' ...
        'adj=%.1f/%.1f, global=%.1e/%.1e, range=%.1e'], ...
        mesh.remeshCount(end),mesh.minimumDx(end), ...
        mesh.minimumDy(end),mesh.maximumCellRatioX(end), ...
        mesh.maximumCellRatioY(end), ...
        mesh.globalGridRatioX(end),mesh.globalGridRatioY(end), ...
        common.physicalRangeViolation(end)), ...
        'FontSize',9);
end

drawnow;
if ~isempty(viz.writer)
    writeVideo(viz.writer,getframe(viz.figure));
end
catch exception
    if initializing
        ipm.output.closeLive(viz,true);
    end
    rethrow(exception);
end
end

function value = optional_history_value(history,name,fallback)
if isfield(history,name) && ~isempty(history.(name))
    value = history.(name)(end);
else
    value = fallback;
end
end

function values = optional_history_series(history,name,reference)
if isfield(history,name) && numel(history.(name)) == numel(reference)
    values = history.(name);
else
    values = NaN(size(reference));
end
end

function viz = boundary_center_live(viz,flow,ops,scale)
% One figure follows the wall and nearby core on the current, possibly
% remeshed, computational grid. Show R_X, the field used by both gauges.
if ~isempty(viz) && isfield(viz,'closed') && viz.closed
    return;
end
if ~isempty(viz) && ~isgraphics(viz.figure,'figure')
    viz.closed = true;
    return;
end

xHalfWidth = max([1.5*flow.trackedWallCoreWidth, ...
    12*min(ops.dxFaces),1]);
yHeight = max([1.5*flow.trackedVerticalCoreWidth, ...
    12*min(ops.dyFaces),1]);
if ~isfinite(xHalfWidth), xHalfWidth = 2; end
if ~isfinite(yHeight), yHeight = 2; end
center = 0;
if ops.rescaling.enabled
    center = ops.rescaling.pinX;
end
xWindow = [max(ops.x(1),center-xHalfWidth), ...
    min(ops.x(end),center+xHalfWidth)];
yWindow = [0,min(ops.y(end),yHeight)];
xSelection = ops.x >= xWindow(1) & ops.x <= xWindow(2);
ySelection = ops.y >= yWindow(1) & ops.y <= yWindow(2);
if nnz(xSelection) < 3 || nnz(ySelection) < 3
    error('ipm:LiveBoundaryCenterWindow', ...
        'The live boundary-center window has fewer than three grid nodes.');
end
x = ops.x(xSelection);
y = ops.y(ySelection);
wallGradient = flow.source(ySelection,xSelection);
[X,Y] = meshgrid(x,y);

if isempty(viz)
    viz = struct();
    viz.figure = figure('Name','IPM boundary center', ...
        'NumberTitle','off','Color','w','Visible','on');
    layout = tiledlayout(viz.figure,2,1,'TileSpacing','compact', ...
        'Padding','compact');
    viz.fieldAxes = nexttile(layout);
    viz.fieldSurface = surf(viz.fieldAxes,X,Y,wallGradient, ...
        'EdgeColor','none');
    view(viz.fieldAxes,2);
    colorbar(viz.fieldAxes);
    ylabel(viz.fieldAxes,'Y');
    viz.wallAxes = nexttile(layout);
    viz.wallLine = plot(viz.wallAxes,x,flow.source(1,xSelection), ...
        'LineWidth',1.5);
    grid(viz.wallAxes,'on');
    xlabel(viz.wallAxes,'X');
    ylabel(viz.wallAxes,'R_X(X,0)');
    viz.writer = [];
else
    set(viz.fieldSurface,'XData',X,'YData',Y, ...
        'ZData',wallGradient,'CData',wallGradient);
    set(viz.wallLine,'XData',x,'YData',flow.source(1,xSelection));
end
xlim(viz.fieldAxes,xWindow);
ylim(viz.fieldAxes,yWindow);
xlim(viz.wallAxes,xWindow);
title(viz.fieldAxes,sprintf( ...
    'Boundary center R_X   tau=%.4f   t=%.4f   remesh=%d', ...
    scale.canonicalTime,scale.physicalTime,ops.remeshCount));
drawnow limitrate;
end
