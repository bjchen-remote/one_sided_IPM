% Plot the complete rescaled density R on a local Cartesian box.
addpath(fileparts(mfilename('fullpath')));

tau = 50;                         % Change to [] for the latest saved frame.
box = [-2 2 0 4];                 % [Xmin Xmax Ymin Ymax].
data = load_polar_data([], tau, box, 201);

figure('Color','w','Name','Local density');
imagesc(data.x, data.y, data.rho);
axis xy equal;
xlim(box(1:2)); ylim(box(3:4));
xlabel('X'); ylabel('Y');
title(sprintf('R, \\tau = %.4f', data.selected.tau),'Color','k');
colorbar;
set(gca,'Color','w','FontName','Times New Roman','FontSize',11, ...
    'Box','off','TickDir','out','LineWidth',1,'XColor','k','YColor','k');
% data.rho is R = W + data.rho_origin, not the centered field data.W.
