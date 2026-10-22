% Keep the full c_r trajectory when selecting one spatial frame.
addpath(fileparts(mfilename('fullpath')));
data = load_polar_data([], [], [-2 2 0 4], 201);

figure('Color','w','Name','Translation rate versus time');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
curve = plot(data.t,data.c_l);
set(curve,'LineWidth',1,'Color',[0 .447 .741]);
xlabel('physical time t'); ylabel('c_r');
set(gca,'Color','w','FontName','Times New Roman','FontSize',11, ...
    'Box','off','TickDir','out','LineWidth',1,'XColor','k','YColor','k');

nexttile;
curve = plot(data.tau,data.c_r);
set(curve,'LineWidth',1,'Color',[0 .447 .741]);
xlabel('rescaled time \tau'); ylabel('c_r');
set(gca,'Color','w','FontName','Times New Roman','FontSize',11, ...
    'Box','off','TickDir','out','LineWidth',1,'XColor','k','YColor','k');
% Physical t approaches its limiting value late in the run; tau keeps growing.
% Both panels connect every recorded point directly, without smoothing.
