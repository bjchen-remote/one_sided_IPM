function outputFile = plot_streamfunction_full_20s(root, outputFile)
%PLOT_STREAMFUNCTION_FULL_20S Reconstruct and plot all stored tau=0..20 frames.
% The result file stores density snapshots; the streamfunction is rebuilt
% from each snapshot with the same Poisson operator as the forward run.

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(root, 'runs', 'streamfunction_geometry_20260928', ...
        'streamfunction_full_20s.png');
end

addpath(fullfile(root, 'solver'));
runFile = fullfile(root, 'runs', 'run_20260927_133548_651', ...
    'result_or_target.mat');
loaded = load(runFile, 'result');
result = loaded.result;
ops = ipm.mesh.build(result.config);

time = result.snapshots.canonicalTime(:);
numberOfFrames = numel(time);
assert(numberOfFrames == 201 && abs(time(1)) < 1e-12 && ...
    abs(time(end)-20) < 1e-10, 'Expected the complete tau=0..20 trajectory.');
assert(max(abs(ops.x-result.snapshots.x{1})) < 1e-8 && ...
    max(abs(ops.y-result.snapshots.y{1})) < 1e-8, ...
    'Rebuilt mesh does not match the stored snapshots.');

xi = linspace(-3,4,220);
eta = linspace(0,2,130);
[XI,ETA] = meshgrid(xi,eta);
midEta = 0.75;
history = NaN(numberOfFrames,numel(xi));
peakPosition = NaN(numberOfFrames,1);
curvatureWidth = NaN(numberOfFrames,1);

frameIndex = round(linspace(1,numberOfFrames,16));
frameImages = cell(numel(frameIndex),1);

for k = 1:numberOfFrames
    rho = result.snapshots.rho{k};
    [u1,~,psi] = ipm.field.velocity(rho,ops);
    wallGradient = rho(1,:)*ops.Dx';
    [~,index] = max(wallGradient);
    index = max(2,min(numel(ops.x)-1,index));
    localX = ops.x(index-1:index+1);
    localQ = wallGradient(index-1:index+1);
    polynomial = polyfit(localX,localQ,2);
    a = -polynomial(2)/(2*polynomial(1));
    P = polyval(polynomial,a);
    C = 2*polynomial(1);
    if ~(isfinite(a) && isfinite(P) && P > 0 && isfinite(C) && C < 0)
        a = ops.x(index);
        P = wallGradient(index);
        C = (wallGradient(index+1)-2*wallGradient(index)+ ...
            wallGradient(index-1))/((ops.x(index+1)-ops.x(index))^2);
    end
    ell = sqrt(max(-2*P/C,eps));
    speedAtPeak = interp1(ops.x,u1(1,:),a,'pchip','extrap');
    centredPsi = (psi + speedAtPeak*ops.y)/(P*ell^2);

    queryX = a + ell*xi;
    queryY = ell*midEta;
    history(k,:) = interp2(ops.x,ops.y,centredPsi,queryX, ...
        queryY*ones(size(queryX)),'linear',NaN);
    if any(frameIndex == k)
        imageIndex = find(frameIndex == k,1);
        frameImages{imageIndex} = interp2(ops.x,ops.y,centredPsi, ...
            a + ell*XI,ell*ETA,'linear',NaN);
    end
    peakPosition(k) = a;
    curvatureWidth(k) = ell;
end

finiteValues = [];
for k = 1:numel(frameImages)
    finiteValues = [finiteValues; frameImages{k}(isfinite(frameImages{k}))]; %#ok<AGROW>
end
finiteValues = sort(finiteValues);
lowIndex = max(1,round(0.01*numel(finiteValues)));
highIndex = min(numel(finiteValues),round(0.99*numel(finiteValues)));
colourLimits = [finiteValues(lowIndex),finiteValues(highIndex)];
colourScale = max(abs(colourLimits));
colourLimits = [-colourScale,colourScale];
levels = linspace(colourLimits(1),colourLimits(2),17);

figureHandle = figure('Visible','off','Color','w', ...
    'Position',[40,40,1600,2300]);
layout = tiledlayout(6,4,'TileSpacing','compact','Padding','compact');
colormap(turbo(256));
for j = 1:numel(frameIndex)
    axisHandle = nexttile(layout);
    contourf(axisHandle,XI,ETA,frameImages{j},levels,'LineColor','none');
    hold(axisHandle,'on');
    contour(axisHandle,XI,ETA,frameImages{j},levels(2:2:end-1), ...
        'LineColor',[0.2,0.2,0.2],'LineWidth',0.25);
    axis(axisHandle,'image');
    xlim(axisHandle,[xi(1),xi(end)]); ylim(axisHandle,[eta(1),eta(end)]);
    caxis(axisHandle,colourLimits);
    title(axisHandle,sprintf('\\tau = %.1f',time(frameIndex(j))));
    xlabel(axisHandle,'\\xi'); ylabel(axisHandle,'\\eta');
    set(axisHandle,'FontSize',8);
end

axisHandle = nexttile(layout,[2,4]);
imagesc(axisHandle,xi,time,history);
set(axisHandle,'YDir','normal');
axis(axisHandle,'tight');
caxis(axisHandle,colourLimits);
hold(axisHandle,'on');
plot(axisHandle,[xi(1),xi(end)], ...
    [time(frameIndex(1)),time(frameIndex(1))],'w-','LineWidth',0.5);
xlabel(axisHandle,'\\xi = (x-a)/\\ell');
ylabel(axisHandle,'\\tau (seconds)');
title(axisHandle,'Complete 0--20 s history: centred \\psi/(P\\ell^2) at \\eta=0.75');
cb = colorbar(axisHandle);
cb.Label.String = '\\psi/(P\\ell^2)';
set(axisHandle,'FontSize',10);

sgtitle(layout,'Full streamfunction evolution from the stored 20-second run', ...
    'FontSize',16,'FontWeight','bold');
exportgraphics(figureHandle,outputFile,'Resolution',180);
close(figureHandle);

metadataFile = replace(outputFile,'.png','_metadata.mat');
save(metadataFile,'time','xi','eta','history','frameIndex', ...
    'peakPosition','curvatureWidth','colourLimits','-v7');
end
