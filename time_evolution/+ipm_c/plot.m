function handles=plot(cart,field,options)
%PLOT Display one Cartesian slice produced by ipm_c.cartesian.
% h=ipm_c.plot(c,'rho'); h=ipm_c.plot(c,'velocity') overlays Darcy arrows.
% options: axes, limits=[cmin cmax], colormap, velocity={'u1','u2'},
% quiverStride=16, label, fontSize=10. Returns ordinary graphics handles;
% callers retain full control of contours, streamlines, layout and export.
% This routine displays data only: it never evolves, extrapolates or changes
% the selected time, coordinate convention or interpolation method.
if nargin<3,options=struct();end
assert(isstruct(cart)&&isfield(cart,'kind')&&strcmp(cart.kind,'ipm_c_cartesian_slice_v1'), ...
    'ipm_c:PlotSlice','Pass a Cartesian slice returned by ipm_c.cartesian.');
field=char(field);velocity=option(options,'velocity',{});
if strcmp(field,'velocity')
    if isempty(velocity),velocity={'u1','u2'};end
    assert(all(isfield(cart,velocity)),'ipm_c:PlotVelocity','The requested velocity components are absent.');
    values=hypot(cart.(velocity{1}),cart.(velocity{2}));label='|u|';
else
    assert(isfield(cart,field),'ipm_c:PlotField','The Cartesian slice does not contain %s.',field);
    values=cart.(field);label=strrep(field,'_','\_');
    if isfield(cart.metadata.fieldUnits,field)
        unit=cart.metadata.fieldUnits.(field);
        if startsWith(unit,'canonical'),label=[label,' (canonical)'];end
    end
end
assert(any(cart.mask,'all'),'ipm_c:PlotEmptyDomain','The requested box contains no mesh interior.');
if isfield(options,'axes'),ax=options.axes;
else,fig=figure('Color','white');ax=axes(fig);end
assert(isgraphics(ax,'axes'),'ipm_c:PlotAxes','options.axes must be an axes handle.');
im=imagesc(ax,cart.x,cart.y,values);im.AlphaData=cart.mask&isfinite(values);
set(ax,'YDir','normal','Color','white','FontName','Times New Roman', ...
    'FontSize',option(options,'fontSize',10),'TickDir','in','LineWidth',1,'Box','on');
axis(ax,'equal');xlim(ax,cart.box(1:2));ylim(ax,cart.box(3:4));
if strcmp(cart.coordinates,'physical'),xlabel(ax,'x');ylabel(ax,'y');
else,xlabel(ax,'X');ylabel(ax,'Y');end
finite=values(cart.mask&isfinite(values));
if isfield(options,'limits'),clim(ax,options.limits);
elseif min(finite)<0&&max(finite)>0,range=max(abs(finite));clim(ax,[-range,range]);end
if isfield(options,'colormap'),colormap(ax,options.colormap);
elseif min(finite)<0&&max(finite)>0
    anchors=[.12,.30,.62;.96,.97,.98;.65,.12,.17];
    colormap(ax,interp1([0,.5,1],anchors,linspace(0,1,257),'linear'));
else,colormap(ax,parula(256));end
cb=colorbar(ax);cb.Label.String=option(options,'label',label);cb.Label.Interpreter='tex';
arrows=gobjects(0);
if ~isempty(velocity)
    assert(iscell(velocity)&&numel(velocity)==2&&all(isfield(cart,velocity)), ...
        'ipm_c:PlotVelocity','velocity must name two interpolated components.');
    stride=option(options,'quiverStride',16);
    validateattributes(stride,{'numeric'},{'scalar','integer','positive'});
    rows=1:stride:numel(cart.y);columns=1:stride:numel(cart.x);
    wasHeld=ishold(ax);hold(ax,'on');
    arrows=quiver(ax,cart.X(rows,columns),cart.Y(rows,columns), ...
        cart.(velocity{1})(rows,columns),cart.(velocity{2})(rows,columns), ...
        'Color',[.1,.1,.1],'LineWidth',.75);
    if ~wasHeld,hold(ax,'off');end
end
handles=struct('figure',ancestor(ax,'figure'),'axes',ax,'image',im,'colorbar',cb,'quiver',arrows);
end

function value=option(options,name,fallback)
if isfield(options,name),value=options.(name);else,value=fallback;end
end
