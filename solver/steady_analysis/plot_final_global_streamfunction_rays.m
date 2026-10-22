function [outputFile,fitTable] = plot_final_global_streamfunction_rays(root,outputFile)
%PLOT_FINAL_GLOBAL_STREAMFUNCTION_RAYS Global tau=20 field and ray fits.
% The rays start at the wall core (a,0). Fits are descriptive power laws
% on the largest sign-consistent segment in the measured outer window.

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(root,'runs','streamfunction_geometry_20260928', ...
        'streamfunction_final_global_rays.png');
end

addpath(fullfile(root,'solver'));
sourceFile = fullfile(root,'runs','steady_20260928_tau20','run_analysis.mat');
m = load(sourceFile);
x = m.xFull(:)';
y = m.yFull(:);
a = double(m.metadata.a);
P = double(m.metadata.P);
ell = double(m.metadata.ell);
U0 = double(m.metadata.uAtPeak);

% This is the centred streamfunction used by the independent polar audit.
% Inner units make the ray fits numerically readable while preserving the
% whole native computational box in the first two panels.
psiCentered = m.psi + U0*y;
psiInner = psiCentered/(P*ell^2);
eta = y/ell;
H = NaN(size(psiInner));
H(2:end,:) = psiInner(2:end,:)./eta(2:end);
u1Inner = (m.u1-U0)/(P*ell);
H(1,:) = -u1Inner(1,:);

% Rays are measured in inner radius r=R/ell and extend to the edge of the
% saved box. Avoid the exact coordinate axes, where the boundary trace is
% a limiting value rather than an interior interpolation point.
anglesDeg = [10,25,35,45,60,90,120,145,165,175];
angles = anglesDeg*pi/180;
numberOfRays = numel(angles);
ray = cell(numberOfRays,1);
fitRows = repmat(struct('angle_deg',NaN,'exponent_p',NaN, ...
    'coefficient_C',NaN,'r2_log_abs',NaN,'fit_rmin_inner',NaN, ...
    'fit_rmax_inner',NaN,'sign',NaN,'points',0,'status','no_fit'), ...
    numberOfRays,1);

fitLower = 4;
fitUpper = 4096;
for j = 1:numberOfRays
    theta = angles(j);
    c = cos(theta);
    s = sin(theta);
    horizontalLimits = [(x(end)-a)/c,(x(1)-a)/c];
    horizontalLimits = horizontalLimits(horizontalLimits > 0 & isfinite(horizontalLimits));
    verticalLimit = y(end)/s;
    limits = [horizontalLimits,verticalLimit];
    maxNativeRadius = 0.995*min(limits);
    maxInnerRadius = maxNativeRadius/ell;
    innerRadius = logspace(log10(0.05),log10(maxInnerRadius),600);
    queryX = a + ell*innerRadius*c;
    queryY = ell*innerRadius*s;
    value = interp2(x,y,psiInner,queryX,queryY,'linear',NaN);
    ray{j} = struct('theta',theta,'angleDeg',anglesDeg(j), ...
        'r',innerRadius,'x',queryX,'y',queryY,'psi',value, ...
        'rMax',maxInnerRadius);

    valid = innerRadius >= fitLower & innerRadius <= min(fitUpper,maxInnerRadius) & ...
        isfinite(value) & abs(value) > 0;
    validValues = value(valid);
    if any(valid)
        threshold = max(1e-12,1e-5*max(abs(validValues)));
        valid = valid & abs(value) >= threshold;
    end
    fitIndex = longest_sign_run(valid,value);
    if numel(fitIndex) >= 12
        logRadius = log(innerRadius(fitIndex));
        logValue = log(abs(value(fitIndex)));
        coefficient = polyfit(logRadius,logValue,1);
        predicted = polyval(coefficient,logRadius);
        residual = logValue-predicted;
        total = logValue-mean(logValue);
        r2 = 1-sum(residual.^2)/max(sum(total.^2),eps);
        fitRows(j).angle_deg = anglesDeg(j);
        fitRows(j).exponent_p = coefficient(1);
        fitRows(j).coefficient_C = sign(median(value(fitIndex)))*exp(coefficient(2));
        fitRows(j).r2_log_abs = r2;
        fitRows(j).fit_rmin_inner = innerRadius(fitIndex(1));
        fitRows(j).fit_rmax_inner = innerRadius(fitIndex(end));
        fitRows(j).sign = sign(median(value(fitIndex)));
        fitRows(j).points = numel(fitIndex);
        fitRows(j).status = 'power_fit';
        ray{j}.fitIndex = fitIndex;
        ray{j}.fit = struct('coefficient',coefficient,'r2',r2, ...
            'sign',fitRows(j).sign);
    else
        ray{j}.fitIndex = [];
        ray{j}.fit = [];
    end
end

fitTable = struct2table(fitRows);
fitCsv = replace(outputFile,'.png','_ray_fits.csv');
writetable(fitTable,fitCsv);
dataFile = replace(outputFile,'.png','_data.mat');
save(dataFile,'x','y','a','P','ell','U0','psiInner','H','anglesDeg','ray','fitTable', ...
    'fitLower','fitUpper','-v7');

figureHandle = figure('Visible','off','Color','w','Position',[30,30,1900,1300]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% Full native computational box.
axisHandle = nexttile(layout);
globalValues = psiInner(isfinite(psiInner));
globalLimits = percentile_limits(globalValues,[1,99]);
contourf(axisHandle,x/1e3,y/1e3,psiInner,35,'LineColor','none');
hold(axisHandle,'on');
caxis(axisHandle,globalLimits);
plot_rays(axisHandle,ray,1e3);
plot(axisHandle,a/1e3,0,'kp','MarkerFaceColor','w','MarkerSize',10);
axis(axisHandle,'image'); axis(axisHandle,'tight');
xlabel(axisHandle,'X/10^3'); ylabel(axisHandle,'Y/10^3');
title(axisHandle,'tau=20 全局中心化流函数 Psi/(P ell^2)');
cb = colorbar(axisHandle); cb.Label.String = 'Psi/(P ell^2)';

% H=Psi/Y removes the wall-forced zero and exposes the global transverse
% structure without changing the ray origin or the underlying field.
axisHandle = nexttile(layout);
hValues = H(isfinite(H));
hLimits = percentile_limits(hValues,[1,99]);
contourf(axisHandle,x/1e3,y/1e3,H,35,'LineColor','none');
hold(axisHandle,'on');
caxis(axisHandle,hLimits);
plot_rays(axisHandle,ray,1e3);
plot(axisHandle,a/1e3,0,'kp','MarkerFaceColor','w','MarkerSize',10);
axis(axisHandle,'image'); axis(axisHandle,'tight');
xlabel(axisHandle,'X/10^3'); ylabel(axisHandle,'Y/10^3');
title(axisHandle,'全局 H=Psi/(P ell Y)（壁面连续延拓）');
cb = colorbar(axisHandle); cb.Label.String = 'H';

% Signed radial profiles retain zero crossings and therefore show which
% rays can or cannot be represented by one signed power law.
axisHandle = nexttile(layout);
hold(axisHandle,'on');
colors = lines(numberOfRays);
for j = 1:numberOfRays
    plot(axisHandle,ray{j}.r,ray{j}.psi,'Color',colors(j,:), ...
        'LineWidth',1.05,'DisplayName',sprintf('theta=%g deg',anglesDeg(j)));
end
yline(axisHandle,0,'k:');
set(axisHandle,'XScale','log');
xlabel(axisHandle,'r=R/ell'); ylabel(axisHandle,'Psi/(P ell^2)');
title(axisHandle,'核心点出发的射线剖面（保留零线穿越）');
grid(axisHandle,'on'); legend(axisHandle,'Location','eastoutside');

% Absolute values and the fitted power law on the selected sign-consistent
% outer segment.
axisHandle = nexttile(layout);
hold(axisHandle,'on');
for j = 1:numberOfRays
    finite = isfinite(ray{j}.psi) & ray{j}.psi ~= 0;
    plot(axisHandle,ray{j}.r(finite),abs(ray{j}.psi(finite)), ...
        'Color',colors(j,:),'LineWidth',1.0, ...
        'DisplayName',sprintf('theta=%g deg, p=%s',anglesDeg(j), ...
            format_fit_exponent(fitRows(j))));
    if ~isempty(ray{j}.fitIndex)
        idx = ray{j}.fitIndex;
        fitValue = abs(fitRows(j).coefficient_C)*ray{j}.r(idx).^fitRows(j).exponent_p;
        plot(axisHandle,ray{j}.r(idx),fitValue,'--','Color',colors(j,:), ...
            'LineWidth',1.3,'HandleVisibility','off');
    end
end
set(axisHandle,'XScale','log','YScale','log');
xlabel(axisHandle,'r=R/ell'); ylabel(axisHandle,'abs(Psi)/(P ell^2)');
title(axisHandle,'中远区幂律尝试: abs(Psi) ~= abs(C) r^p');
grid(axisHandle,'on'); legend(axisHandle,'Location','eastoutside');

exportgraphics(figureHandle,outputFile,'Resolution',180);
close(figureHandle);
end

function fitIndex = longest_sign_run(valid,value)
% Return the longest contiguous run with valid, non-changing sign.
indices = find(valid);
if isempty(indices)
    fitIndex = [];
    return;
end
signs = sign(value(indices));
breaks = [1; find(diff(signs) ~= 0)+1; numel(indices)+1];
lengths = diff(breaks);
[~,which] = max(lengths);
fitIndex = indices(breaks(which):breaks(which+1)-1);
end

function limits = percentile_limits(values,percentiles)
values = sort(values(:));
if isempty(values)
    limits = [-1,1];
    return;
end
lo = values(max(1,round(percentiles(1)/100*numel(values))));
hi = values(min(numel(values),round(percentiles(2)/100*numel(values))));
if ~(isfinite(lo) && isfinite(hi)) || lo == hi
    span = max(abs([lo,hi,1]));
    limits = [-span,span];
else
    limits = [lo,hi];
end
end

function plot_rays(axisHandle,ray,scale)
colors = lines(numel(ray));
for j = 1:numel(ray)
    plot(axisHandle,ray{j}.x/scale,ray{j}.y/scale,'-', ...
        'Color',colors(j,:),'LineWidth',0.8,'HandleVisibility','off');
end
end

function textValue = format_fit_exponent(row)
if strcmp(row.status,'power_fit')
    textValue = sprintf('%.2f',row.exponent_p);
else
    textValue = 'n/a';
end
end
