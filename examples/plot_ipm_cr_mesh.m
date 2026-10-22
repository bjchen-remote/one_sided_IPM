function report=plot_ipm_cr_mesh(outputDir,checkpointFolder)
%PLOT_IPM_CR_MESH Native-node surfaces of the accepted hybrid trajectory.
% No time advance, spatial reconstruction or uniform display resampling.
% Global density is shown with its SAVED scalar origin removed for display;
% the exact scalar needed to reconstruct rho is recorded alongside the plot.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
stylePath=fullfile(fileparts(mfilename('fullpath')),'private','mesh_style.m');
styleIdentity=struct('path',stylePath,'sha256',ipm_cr.io.sha(stylePath));
if nargin<1||isempty(outputDir)
    outputDir=fullfile(root,'test','results','mesh_views_tau200');
end
if nargin<2||isempty(checkpointFolder)
    checkpointFolder=fullfile(root,'test','results','explicit_tau200_20261021_v1','native');
end
assert(~isfile(fullfile(outputDir,'plots.json')),'ipm_cr:PlotPreserve','Use a fresh output directory.');
if ~isfolder(outputDir),mkdir(outputDir);end
% Historical v1 is imported explicitly and only for read-only plotting.
% Current v2 uses a portable slot basename validated by the checkpoint API.
header=jsondecode(fileread(fullfile(checkpointFolder,'latest.json')));
if isfield(header,'kind')&&strcmp(header.kind,'ipm_cr_explicit_accepted_pointer_v1')
    mode='legacy';
else
    mode='load';
end
[cp,pointer]=ipm_cr.io.checkpoint(mode,checkpointFolder);
if strcmp(mode,'legacy'),checkpointPath=pointer.path;
else,checkpointPath=fullfile(checkpointFolder,pointer.path);end
% This lossless local adapter changes field names only; no context is
% constructed and no RHS, background, interpolation or time step is run.
state=ipm_cr.internal.state(cp.state,'new');g=cp.geometry;c=g.core;p=g.polar;
assert(state.scale.canonicalTime==pointer.tau&&state.step==pointer.sequence, ...
    'ipm_cr:PlotState','The displayed state must equal the accepted pointer.');
blue=[.12,.31,.48];orange=[.64,.29,.16];gray=[.62,.65,.68];
rho0=state.scale.densityOrigin;W=state.W;q=state.rhs.q;
rhoRange=[min([W.core(:);W.polar(:)]),max([W.core(:);W.polar(:)])]+rho0;
rows=p.activeRows;cols=p.activeColumns;
wc=W.core;wc(~c.active)=NaN;
wp=W.polar(rows,cols);
allW=[W.core(c.active);W.polar(p.active)];zlimW=[min(allW),max(allW)];
zlimW=zlimW+[-1,1]*.035*diff(zlimW);

% Two honest global-coordinate views: a true physical disk and a clearly
% labelled log-radius display. Every annular row and angular column is kept.
fig=canvas([11.2,5.6]);layout=panels(fig);
for k=1:2
    ax=nexttile(layout);hold(ax,'on');style(fig,ax);
    if k==1
        XC=c.X;YC=c.Y;XP=p.X(rows,cols);YP=p.Y(rows,cols);
    else
        [XC,YC]=compress_radius(c.X,c.Y);
        [XP,YP]=compress_radius(p.X(rows,cols),p.Y(rows,cols));
    end
    hc=mesh(ax,XC,YC,wc,'FaceColor','none','EdgeColor',blue,'LineWidth',.45);
    hp=mesh(ax,XP,YP,wp,'FaceColor','none','EdgeColor',orange,'LineWidth',.22);
    xlabel(ax,'X');ylabel(ax,'Y');zlabel(ax,'\rho-\rho(0,0)');
    if k==1
        title(ax,'(a) Physical coordinates','FontWeight','normal');
        xlim(ax,[-p.r(end),p.r(end)]);ylim(ax,[0,p.r(end)]);
        ax.XAxis.Exponent=10;ax.YAxis.Exponent=10;
    else
        title(ax,'(b) Log-radius coordinates','FontWeight','normal');
        xlabel(ax,'s cos\theta');ylabel(ax,'s sin\theta');
        S=log10(1+p.r(end));xlim(ax,[-S,S]);ylim(ax,[0,S]);
        ax.XTick=-10:5:10;ax.YTick=0:2:10;
    end
    zlim(ax,zlimW);pbaspect(ax,[2,1,1]);view(ax,-39,29);grid(ax,'on');
    ax.GridAlpha=.12;ax.ZAxis.Exponent=3;
    if k==2,mesh_legend(ax,blue,orange);end
    style(fig,ax);
end
annotation(fig,'textbox',[.06,.014,.88,.07],'String', ...
    ['\rho(0,0)=',sprintf('%.11g',rho0),';   s=log_{10}(1+r).   All active native nodes are retained.'], ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Times New Roman','FontSize',10);
globalFiles=save_plot(fig,outputDir,'rho_global_mesh');close(fig);

% q=rho_X is already stored on both actual patches. Here mesh() plots
% every selected native node; overlap is deliberately shown TWICE, exactly
% as both active interiors are independently evolved.
fig=canvas([11.2,5.6]);layout=panels(fig);
ax=nexttile(layout);hold(ax,'on');style(fig,ax);
qcore=q.core;qcore(~c.active)=NaN;
nearRows=rows(p.r(rows)<=3);qValues=[q.core(c.active);reshape(q.polar(nearRows,:),[],1)];
qRange=[min(qValues),max(qValues)];
qZlim=[min(0,qRange(1))-.025*diff(qRange),qRange(2)+.06*diff(qRange)];
hc=mesh(ax,c.X,c.Y,qcore,'FaceColor','none','EdgeColor',blue,'LineWidth',.43);
title(ax,'(a) Cartesian core','FontWeight','normal');
xlabel(ax,'X');ylabel(ax,'Y');zlabel(ax,'\rho_X');
xlim(ax,[-1.2,1.2]);ylim(ax,[0,1.2]);zlim(ax,qZlim);
view(ax,-38,30);pbaspect(ax,[2,1,1]);grid(ax,'on');ax.GridAlpha=.12;style(fig,ax);
ax=nexttile(layout);hold(ax,'on');style(fig,ax);
hp=mesh(ax,p.X(nearRows,:),p.Y(nearRows,:),q.polar(nearRows,:), ...
    'FaceColor','none','EdgeColor',orange,'LineWidth',.30);
hc=mesh(ax,c.X,c.Y,qcore,'FaceColor','none','EdgeColor',blue,'LineWidth',.48);
title(ax,'(b) Core and annulus','FontWeight','normal');
xlabel(ax,'X');ylabel(ax,'Y');zlabel(ax,'\rho_X');
xlim(ax,[-2.7,2.7]);ylim(ax,[0,2.7]);zlim(ax,qZlim);
view(ax,-38,31);pbaspect(ax,[2,1,1]);grid(ax,'on');ax.GridAlpha=.12;
mesh_legend(ax,blue,orange);style(fig,ax);
annotation(fig,'textbox',[.06,.014,.88,.07],'String', ...
    'Physical coordinates; no display resampling. Both independently evolved overlap meshes are shown.', ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Times New Roman','FontSize',10);
qFiles=save_plot(fig,outputDir,'rhox_core_transition_mesh');close(fig);

% An actual-grid plan view makes the collars and the fact that there is
% no single matching seam visible, without using a schematic replacement.
fig=canvas([9.0,5.6]);layout=panels(fig);
for k=1:2
    ax=nexttile(layout);hold(ax,'on');style(fig,ax);
    z=zeros(size(c.X));z(c.active)=NaN;
    mesh(ax,c.X,c.Y,z,'FaceColor','none','EdgeColor',gray,'LineWidth',.35,'LineStyle',':');
    mesh(ax,p.X(1:g.ghostWidth,:),p.Y(1:g.ghostWidth,:),zeros(g.ghostWidth,p.ntheta), ...
        'FaceColor','none','EdgeColor',gray,'LineWidth',.40,'LineStyle',':');
    z=zeros(size(c.X));z(~c.active)=NaN;
    mesh(ax,c.X,c.Y,z,'FaceColor','none','EdgeColor',blue,'LineWidth',.40);
    rr=rows(p.r(rows)<=2.3);
    mesh(ax,p.X(rr,:),p.Y(rr,:),zeros(numel(rr),p.ntheta), ...
        'FaceColor','none','EdgeColor',orange,'LineWidth',.35);
    xb=c.activeXBounds;yt=c.activeYBounds(end);
    bc=plot(ax,[xb(1),xb(1),xb(2),xb(2)],[0,yt,yt,0], ...
        'Color',blue,'LineWidth',1.8);
    t=linspace(0,pi,401);ra=p.innerActiveRadius;
    bp=plot(ax,ra*cos(t),ra*sin(t),'Color',orange,'LineWidth',1.8);
    plot(ax,[-1.2,-1.2,1.2,1.2],[0,1.2,1.2,0],':','Color',gray,'LineWidth',1.2);
    plot(ax,p.r(1)*cos(t),p.r(1)*sin(t),':','Color',gray,'LineWidth',1.2);
    plot(ax,[-2,2],[0,0],'Color',[.18,.18,.18],'LineWidth',1.1);
    view(ax,2);axis(ax,'equal');xlabel(ax,'X');ylabel(ax,'Y');
    if k==1
        title(ax,'(a) Active grids and ghost collars','FontWeight','normal');
        xlim(ax,[-1.5,1.5]);ylim(ax,[0,1.5]);
        legend(ax,[bc,bp],{'Active core boundary','First active polar ring'}, ...
            'Location','northoutside','Orientation','vertical','Box','off');
    else
        title(ax,'(b) Right overlap detail','FontWeight','normal');
        xlim(ax,[.60,1.26]);ylim(ax,[0,.66]);
    end
    style(fig,ax);
end
annotation(fig,'textbox',[.05,.014,.9,.07],'String', ...
    sprintf('h_X=%.6g; h_Y=%.6g; active core edge=%.6g; polar active radius=%.6f; gray dotted: %d ghost layers.', ...
        c.dx,c.dy,c.activeXBounds(end),p.innerActiveRadius,g.ghostWidth), ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Times New Roman','FontSize',9.5);
gridFiles=save_plot(fig,outputDir,'transport_grid_detail');close(fig);

report=struct('kind','accepted_native_hybrid_grid_mesh_figures_v1', ...
    'checkpoint',pointer,'tau',state.scale.canonicalTime,'rhoOrigin',rho0, ...
    'totalDensityRange',rhoRange,'centeredDensityRange',[min(allW),max(allW)], ...
    'source','Only saved W and rhs.q (or their historical rho/flow.q aliases); no RHS or time evolution is run.', ...
    'checkpointReadMode',mode, ...
    'style','journal-minimal; canvas enlarged for two native-grid panels', ...
    'styleHelper',styleIdentity, ...
    'globalDisplay','Height W=rho-rho(0,0); left actual X,Y, right s*cos(theta),s*sin(theta), s=log10(1+r). All active polar rows/angles and active core nodes are retained.', ...
    'nearDisplay','Saved q=rho_X at native nodes. Only physical view limits and r<=3 annular row selection; no interpolation, decimation or fabricated rectangular display grid.', ...
    'overlapDisplay','Both active patches shown simultaneously; their values are not averaged or shifted apart. Gray dotted collars are algebraic ghost nodes, not evolved DOFs.', ...
    'poissonGridShown',false,'plotChangesNumericalState',false, ...
    'coreActiveBounds',[c.activeXBounds(:)',c.activeYBounds(:)'], ...
    'displayedRhoXRange',qRange, ...
    'polarInnerHaloRadius',p.r(1),'polarInnerActiveRadius',p.innerActiveRadius, ...
    'ghostWidth',g.ghostWidth,'global',globalFiles,'rhox',qFiles,'grid',gridFiles, ...
    'captions',{{ ...
        sprintf('Global density at tau=%.17g. The vertical offset is rho(0,0)=%.17g, so rho equals the plotted height plus this constant. Left: physical coordinates over the full radius %.17g. Right: explicitly compressed log-radius coordinates. No node decimation.',state.scale.canonicalTime,rho0,p.r(end)), ...
        sprintf('rho_X at tau=%.17g on the actual uniform Cartesian core and mapped polar annulus. Blue: active Cartesian mesh. Orange: active polar mesh. The right panel includes the complete overlap twice because both patches evolve there.',state.scale.canonicalTime), ...
        'Top view of actual transport nodes: solid blue/orange boundaries delimit the active Cartesian rectangle and first active annular ring. Gray dotted lines show three algebraic ghost layers. The overlap is a two-dimensional region, not a single fitted seam. Chebyshev Poisson nodes are a separate discretization and are not drawn.'}});
report.script=struct('path',mfilename('fullpath'),'sha256',ipm_cr.io.sha([mfilename('fullpath'),'.m']));
report.acceptedCheckpointUnchanged=strcmp(ipm_cr.io.sha(checkpointPath),pointer.sha256);
assert(report.acceptedCheckpointUnchanged,'ipm_cr:PlotMutation');
assert(strcmp(ipm_cr.io.sha(stylePath),styleIdentity.sha256), ...
    'ipm_cr:PlotStyleChanged','Local style helper changed during plotting.');
fid=fopen(fullfile(outputDir,'plots.json'),'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(report,'PrettyPrint',true),'char');
fprintf('MESH FIGURES tau=%.17g rho0=%.17g rho range=[%.17g,%.17g] %s\n', ...
    report.tau,rho0,rhoRange,outputDir);
end

function fig=canvas(sizeInches)
fig=figure('Visible','off','Color','white','Units','inches', ...
    'Position',[1,1,sizeInches],'InvertHardcopy','off');
end
function layout=panels(fig)
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
layout.OuterPosition=[.025,.14,.95,.81];
end
function mesh_legend(ax,blue,orange)
hc=plot3(ax,NaN,NaN,NaN,'Color',blue,'LineWidth',1.3);
hp=plot3(ax,NaN,NaN,NaN,'Color',orange,'LineWidth',1.3);
legend(ax,[hc,hp],{'Cartesian active','Polar active'}, ...
    'Location','northoutside','Orientation','horizontal','Box','off');
end
function style(fig,ax)
mesh_style(fig,ax);
ax.FontSize=10;ax.XLabel.FontSize=11;ax.YLabel.FontSize=11;ax.ZLabel.FontSize=11;
ax.Title.FontSize=11;ax.ZColor=[.2,.2,.2];ax.Box='on';ax.LineWidth=.7;
end
function [x,y]=compress_radius(X,Y)
r=hypot(X,Y);scale=ones(size(r));positive=r>0;
scale(positive)=log10(1+r(positive))./r(positive);
x=X.*scale;y=Y.*scale;
end
function files=save_plot(fig,folder,name)
drawnow;
files=struct('png',fullfile(folder,[name,'.png']), ...
    'pdf',fullfile(folder,[name,'.pdf']),'fig',fullfile(folder,[name,'.fig']));
% Export the complete fixed-size canvas, including the reserved caption
% margin. Tight export can clip 3-D axis/legend text at the figure border.
fig.PaperUnits='inches';fig.PaperSize=fig.Position(3:4);
fig.PaperPosition=[0,0,fig.Position(3:4)];fig.PaperPositionMode='manual';
print(fig,files.png,'-dpng','-r240');
print(fig,files.pdf,'-dpdf','-r240');
savefig(fig,files.fig);
end
