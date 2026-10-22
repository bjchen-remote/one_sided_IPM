% Plot Darcy velocity after subtracting its reference-point translation.
addpath(fileparts(mfilename('fullpath')));

tau = [];                         % Latest saved frame; set a number if needed.
box = [-2 2 0 4];
data = load_polar_data([], tau, box, 201);
rows = 1:10:size(data.X,1);        % Arrow spacing only; the data stay intact.
cols = 1:10:size(data.X,2);

figure('Color','w','Name','Relative Darcy velocity');
quiver(data.X(rows,cols), data.Y(rows,cols), ...
    data.u1_relative(rows,cols), data.u2_relative(rows,cols), ...
    'Color',[0 .447 .741],'LineWidth',1);
axis xy equal;
xlim(box(1:2)); ylim(box(3:4));
xlabel('X'); ylabel('Y');
title(sprintf('Relative Darcy velocity, \\tau = %.4f', data.selected.tau),'Color','k');
set(gca,'Color','w','FontName','Times New Roman','FontSize',11, ...
    'Box','off','TickDir','out','LineWidth',1,'XColor','k','YColor','k');
