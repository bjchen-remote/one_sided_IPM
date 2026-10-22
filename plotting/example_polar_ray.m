% Plot a radial profile on the saved angular ray nearest the requested theta.
addpath(fileparts(mfilename('fullpath')));

tau = 50;                         % Change to [] for the latest saved frame.
theta = pi/4;                     % Requested angle in radians.
rWindow = [0 4];                  % Displayed radial interval.
data = load_polar_data([], tau);
p = data.polar;
field = p.rho;                    % Change to another saved polar field.
[~, j] = min(abs(p.theta - theta));
mask = p.r >= rWindow(1) & p.r <= rWindow(2);
% Use the actual nearest angular node: no smoothing or interpolation.

figure('Color','w','Name','Polar ray profile');
plot(p.r(mask), field(mask,j), 'LineWidth',1);
xlim(rWindow);
xlabel('r'); ylabel('R');ch
title(sprintf('\\theta = %.4g, \\tau = %.4f',p.theta(j),data.selected.tau),'Color','k');
set(gca,'Color','w','FontName','Times New Roman','FontSize',11, ...
    'Box','off','TickDir','out','LineWidth',1,'XColor','k','YColor','k');
