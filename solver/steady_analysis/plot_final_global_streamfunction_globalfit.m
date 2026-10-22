function [outputFile,fitTable] = plot_final_global_streamfunction_globalfit(root,outputFile)
%PLOT_FINAL_GLOBAL_STREAMFUNCTION_GLOBALFIT Fit the whole ray, not only its tail.
% The model is bounded at the core by construction:
%
%   Psi(r) = C r^2 [1+(r/rc)^q]^((p-2)/q)
%            [2 exp(-(r/r0)^s)-1],
%
% with the last factor included only for rays that cross the curved zero
% line. Thus Psi~C r^2 at the core and Psi~const*r^p in the far field.
% The fit objective uses every logarithmically spaced sample through the
% edge of the saved box, with an asinh residual to retain both core and
% far-field information.

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(root,'runs','streamfunction_geometry_20260928', ...
        'streamfunction_final_global_rays_globalfit.png');
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

psiCentered = m.psi + U0*y;
psiInner = psiCentered/(P*ell^2);
eta = y/ell;
H = NaN(size(psiInner));
H(2:end,:) = psiInner(2:end,:)./eta(2:end);
u1Inner = (m.u1-U0)/(P*ell);
H(1,:) = -u1Inner(1,:);

anglesDeg = [10,25,35,45,60,90,120,145,165,175];
angles = anglesDeg*pi/180;
numberOfRays = numel(angles);
fitRows = repmat(struct('angle_deg',NaN,'near_exponent',2, ...
    'far_exponent',NaN,'crossover_rc_inner',NaN, ...
    'transition_order_q',NaN,'zero_radius_inner',NaN, ...
    'zero_transition_order_s',NaN,'zero_crossings',0, ...
    'global_asinh_rms',NaN,'global_relative_rms',NaN, ...
    'global_asinh_r2',NaN,'points',0,'status',""),numberOfRays,1);
ray = cell(numberOfRays,1);

for j = 1:numberOfRays
    theta = angles(j);
    c = cos(theta);
    s = sin(theta);
    horizontalLimits = [(x(end)-a)/c,(x(1)-a)/c];
    horizontalLimits = horizontalLimits(horizontalLimits > 0 & isfinite(horizontalLimits));
    verticalLimit = y(end)/s;
    maxNativeRadius = 0.995*min([horizontalLimits,verticalLimit]);
    maxInnerRadius = maxNativeRadius/ell;

    % The smallest radius is above the first reliably resolved interior
    % shell; the analytic model itself supplies the exact r=0 limit.
    innerRadius = logspace(log10(0.05),log10(maxInnerRadius),700);
    queryX = a + ell*innerRadius*c;
    queryY = ell*innerRadius*s;
    value = interp2(x,y,psiInner,queryX,queryY,'linear',NaN);

    nearMask = innerRadius <= 0.3 & isfinite(value);
    if nnz(nearMask) >= 5
        nearCoeff = median(value(nearMask)./innerRadius(nearMask).^2);
    else
        nearCoeff = value(1)/innerRadius(1)^2;
    end
    dataScale = max(abs(value(isfinite(value))));
    if ~(isfinite(nearCoeff) && abs(nearCoeff) > 0)
        nearCoeff = sign(value(find(isfinite(value),1)))*1e-6;
    end
    if ~(isfinite(dataScale) && dataScale > 0)
        dataScale = 1;
    end
    zeroCrossings = count_sign_crossings(innerRadius,value);
    hasZero = zeroCrossings >= 1;
    start = initial_parameters(innerRadius,value,nearCoeff,dataScale,hasZero);
    [parameters,model,objective] = fit_global_ray( ...
        innerRadius,value,nearCoeff,dataScale,hasZero,start);

    finite = isfinite(value) & isfinite(model);
    asinhScale = 1e-3;
    observedAsinh = asinh(value(finite)/dataScale/asinhScale);
    fittedAsinh = asinh(model(finite)/dataScale/asinhScale);
    residual = fittedAsinh-observedAsinh;
    total = observedAsinh-mean(observedAsinh);
    globalAsinhRms = sqrt(mean(residual.^2));
    globalAsinhR2 = 1-sum(residual.^2)/max(sum(total.^2),eps);
    globalRelativeRms = sqrt(mean((model(finite)-value(finite)).^2))/ ...
        max(sqrt(mean(value(finite).^2)),eps);

    decoded = decode_parameters(parameters,nearCoeff/dataScale,hasZero);
    fitRows(j).angle_deg = anglesDeg(j);
    fitRows(j).far_exponent = decoded.p;
    fitRows(j).crossover_rc_inner = decoded.rc;
    fitRows(j).transition_order_q = decoded.q;
    fitRows(j).zero_radius_inner = decoded.zeroRadius;
    fitRows(j).zero_transition_order_s = decoded.s;
    fitRows(j).zero_crossings = zeroCrossings;
    fitRows(j).global_asinh_rms = globalAsinhRms;
    fitRows(j).global_relative_rms = globalRelativeRms;
    fitRows(j).global_asinh_r2 = globalAsinhR2;
    fitRows(j).points = nnz(finite);
    if hasZero && zeroCrossings > 1
        fitRows(j).status = "one-zero model; multiple crossings";
    elseif hasZero
        fitRows(j).status = "global crossover + one zero";
    else
        fitRows(j).status = "global crossover";
    end
    ray{j} = struct('theta',theta,'angleDeg',anglesDeg(j), ...
        'r',innerRadius,'x',queryX,'y',queryY,'psi',value, ...
        'rMax',maxInnerRadius,'model',model,'parameters',decoded, ...
        'objective',objective,'zeroCrossings',zeroCrossings);
end

fitTable = struct2table(fitRows);
fitCsv = replace(outputFile,'.png','_ray_global_fits.csv');
writetable(fitTable,fitCsv);
dataFile = replace(outputFile,'.png','_data.mat');
save(dataFile,'x','y','a','P','ell','U0','psiInner','H','anglesDeg','ray', ...
    'fitTable','-v7');

figureHandle = figure('Visible','off','Color','w','Position',[30,30,1900,1300]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

axisHandle = nexttile(layout);
globalLimits = percentile_limits(psiInner(isfinite(psiInner)),[1,99]);
contourf(axisHandle,x/1e3,y/1e3,psiInner,35,'LineColor','none');
hold(axisHandle,'on'); caxis(axisHandle,globalLimits);
plot_rays(axisHandle,ray,1e3);
plot(axisHandle,a/1e3,0,'kp','MarkerFaceColor','w','MarkerSize',10);
axis(axisHandle,'image'); axis(axisHandle,'tight');
xlabel(axisHandle,'X/10^3'); ylabel(axisHandle,'Y/10^3');
title(axisHandle,'tau=20 全局中心化流函数 Psi/(P ell^2)');
cb = colorbar(axisHandle); cb.Label.String = 'Psi/(P ell^2)';

axisHandle = nexttile(layout);
hLimits = percentile_limits(H(isfinite(H)),[1,99]);
contourf(axisHandle,x/1e3,y/1e3,H,35,'LineColor','none');
hold(axisHandle,'on'); caxis(axisHandle,hLimits);
plot_rays(axisHandle,ray,1e3);
plot(axisHandle,a/1e3,0,'kp','MarkerFaceColor','w','MarkerSize',10);
axis(axisHandle,'image'); axis(axisHandle,'tight');
xlabel(axisHandle,'X/10^3'); ylabel(axisHandle,'Y/10^3');
title(axisHandle,'全局 H=Psi/(P ell Y)（壁面连续延拓）');
cb = colorbar(axisHandle); cb.Label.String = 'H';

colors = lines(numberOfRays);
axisHandle = nexttile(layout);
hold(axisHandle,'on');
for j = 1:numberOfRays
    plot(axisHandle,ray{j}.r,ray{j}.psi,'Color',colors(j,:), ...
        'LineWidth',1.0,'DisplayName',sprintf('theta=%g deg',anglesDeg(j)));
    plot(axisHandle,ray{j}.r,ray{j}.model,'--','Color',colors(j,:), ...
        'LineWidth',1.2,'HandleVisibility','off');
end
yline(axisHandle,0,'k:');
set(axisHandle,'XScale','log');
xlabel(axisHandle,'r=R/ell'); ylabel(axisHandle,'Psi/(P ell^2)');
title(axisHandle,'全半径射线：实线数据，虚线全局有界拟合');
grid(axisHandle,'on'); legend(axisHandle,'Location','eastoutside');

axisHandle = nexttile(layout);
hold(axisHandle,'on');
for j = 1:numberOfRays
    finite = isfinite(ray{j}.psi) & ray{j}.psi ~= 0;
    finiteModel = isfinite(ray{j}.model) & ray{j}.model ~= 0;
    label = sprintf('theta=%g, p=%.2f, err=%.1f%%',anglesDeg(j), ...
        fitRows(j).far_exponent,100*fitRows(j).global_relative_rms);
    plot(axisHandle,ray{j}.r(finite),abs(ray{j}.psi(finite)), ...
        'Color',colors(j,:),'LineWidth',1.0,'DisplayName',label);
    plot(axisHandle,ray{j}.r(finiteModel),abs(ray{j}.model(finiteModel)), ...
        '--','Color',colors(j,:),'LineWidth',1.2,'HandleVisibility','off');
end
set(axisHandle,'XScale','log','YScale','log');
xlabel(axisHandle,'r=R/ell'); ylabel(axisHandle,'abs(Psi)/(P ell^2)');
title(axisHandle,'全局拟合的幅度对照（原点强制为 O(r^2)）');
grid(axisHandle,'on'); legend(axisHandle,'Location','eastoutside');

exportgraphics(figureHandle,outputFile,'Resolution',180);
close(figureHandle);
end

function start = initial_parameters(r,value,nearCoeff,dataScale,hasZero)
nearCoeff = nearCoeff/dataScale;
outer = isfinite(value) & r >= 4;
if nnz(outer) >= 10
    outerIndex = find(outer);
    outerIndex = outerIndex(abs(value(outerIndex)) > 1e-8*max(abs(value(outer))));
else
    outerIndex = find(isfinite(value));
end
if numel(outerIndex) >= 10
    p0 = polyfit(log(r(outerIndex)),log(abs(value(outerIndex))/dataScale),1);
    p0 = max(-0.5,min(3.2,p0(1)));
else
    p0 = 1.5;
end
rc0 = 1;
q0 = 2;
if hasZero
    transition = first_transition(r,value);
    if isfinite(transition)
        r00 = transition/(log(2)^(1/2));
    else
        r00 = 100;
    end
    s0 = 2;
    start = [log(max(abs(nearCoeff),1e-12)),log(rc0),log(q0),p0, ...
        log(max(r00,0.05)),log(s0)];
else
    start = [log(max(abs(nearCoeff),1e-12)),log(rc0),log(q0),p0];
end
end

function [best,model,bestObjective] = fit_global_ray(r,value,nearCoeff,dataScale,hasZero,start)
starts = start;
if hasZero
    starts = [starts; start+[0,log(4),0,0,log(4),0]; ...
        start+[0,log(.25),0,0,log(.25),0]];
else
    starts = [starts; start+[0,log(4),0,0];start+[0,log(.25),0,0]];
end
bestObjective = inf;
best = start;
for k = 1:size(starts,1)
    objective = @(z) global_objective(z,r,value/dataScale,nearCoeff,hasZero);
    options = optimset('Display','off','MaxIter',12000,'MaxFunEvals',40000, ...
        'TolX',1e-7,'TolFun',1e-8);
    [candidate,candidateObjective] = fminsearch(objective,starts(k,:),options);
    if candidateObjective < bestObjective
        best = candidate;
        bestObjective = candidateObjective;
    end
end
decoded = decode_parameters(best,nearCoeff,hasZero);
model = dataScale*global_model(r,decoded,hasZero);
end

function loss = global_objective(z,r,target,nearCoeff,hasZero)
decoded = decode_parameters(z,nearCoeff,hasZero);
if decoded.invalid
    loss = 1e8 + decoded.penalty;
    return;
end
model = global_model(r,decoded,hasZero);
sigma = 1e-3;
residual = asinh(model/sigma)-asinh(target/sigma);
loss = mean(residual(isfinite(residual)).^2);
if ~isfinite(loss)
    loss = 1e8;
end
end

function decoded = decode_parameters(z,nearCoeff,hasZero)
decoded.invalid = false;
decoded.penalty = 0;
decoded.C = sign(nearCoeff)*exp(z(1));
decoded.rc = exp(z(2));
decoded.q = exp(z(3));
decoded.p = z(4);
decoded.zeroRadius = Inf;
decoded.s = Inf;
if hasZero
    decoded.zeroRadius = exp(z(5));
    decoded.s = exp(z(6));
end
bounds = [decoded.rc,decoded.q,decoded.p];
if hasZero
    bounds = [bounds,decoded.zeroRadius,decoded.s];
end
lower = [0.02,0.2,-1.0];
upper = [1e7,12,3.5];
if hasZero
    lower = [lower,0.02,0.2];
    upper = [upper,1e8,12];
end
violation = max(lower-bounds,0).^2+max(bounds-upper,0).^2;
decoded.penalty = sum(violation);
if decoded.penalty > 0 || any(~isfinite(bounds))
    decoded.invalid = true;
end
end

function model = global_model(r,decoded,hasZero)
transition = r.^2.*(1+(r/decoded.rc).^decoded.q).^((decoded.p-2)/decoded.q);
if hasZero
    signFactor = 2*exp(-(r/decoded.zeroRadius).^decoded.s)-1;
else
    signFactor = 1;
end
model = decoded.C*transition.*signFactor;
end

function count = count_sign_crossings(r,value)
valid = isfinite(value) & abs(value) > 1e-5*max(abs(value(isfinite(value))));
indices = find(valid);
if numel(indices) < 2
    count = 0;
    return;
end
count = sum(diff(sign(value(indices))) ~= 0);
end

function transition = first_transition(r,value)
valid = isfinite(value) & abs(value) > 1e-5*max(abs(value(isfinite(value))));
indices = find(valid);
if numel(indices) < 2
    transition = NaN;
    return;
end
which = find(diff(sign(value(indices))) ~= 0,1);
if isempty(which)
    transition = NaN;
else
    transition = sqrt(r(indices(which))*r(indices(which+1)));
end
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
