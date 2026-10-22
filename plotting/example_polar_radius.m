% Plot an angular profile on the saved radial ring nearest the requested r.
addpath(fileparts(mfilename('fullpath')));

tau = 50;                         % Change to [] for the latest saved frame.
radius = 1;                       % Requested radius.
data = load_polar_data([], tau);
p = data.polar;
field = p.rho;                    % Change to another saved polar field.
[~, i] = min(abs(p.r - radius));
% Use the actual nearest radial node: no smoothing or interpolation.

figure('Color','w','Name','Polar ring profile');
plot(p.theta, field(i,:), 'LineWidth',1);
xlabel('\theta'); ylabel('R');
title(sprintf('r = %.4g, \\tau = %.4f',p.r(i),data.selected.tau),'Color','k');
set(gca,'Color','w','FontName','Times New Roman','FontSize',11, ...
    'Box','off','TickDir','out','LineWidth',1,'XColor','k','YColor','k');
