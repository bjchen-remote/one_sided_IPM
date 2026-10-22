function mesh_style(fig,ax)
%MESH_STYLE Local serif style for native-grid figures, without global defaults.
% Only presentation properties are changed: never data, coordinates, limits,
% ticks, scales, color limits, labels, titles or existing curve attributes.
assert(isscalar(fig)&&isgraphics(fig,'figure')&&isscalar(ax)&&isgraphics(ax,'axes'), ...
    'ipm_cr:MeshStyleTarget','Use one figure and its Cartesian axes.');
assert(ancestor(ax,'figure')==fig,'ipm_cr:MeshStyleTarget','Axes belong to another figure.');
candidates={'STIXGeneral','DejaVu Serif','Times New Roman','Times'};
available=listfonts;font=candidates{1};
for j=1:numel(candidates)
    if any(strcmpi(candidates{j},available)),font=candidates{j};break;end
end
fig.Color='white';
set(ax,'FontName',font,'FontSize',9,'LineWidth',.7,'TickDir','out', ...
    'Color','white','XColor',[.2,.2,.2],'YColor',[.2,.2,.2], ...
    'Box','off','XGrid','off','YGrid','off');
labels=[ax.XLabel,ax.YLabel,ax.ZLabel];
for j=1:numel(labels)
    set(labels(j),'FontName',font,'FontSize',10,'Color',[34,34,34]/255);
end
set(ax.Title,'FontName',font,'FontSize',10);
end
