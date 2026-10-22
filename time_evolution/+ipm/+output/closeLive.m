function closeLive(viz,closeFigure)
%IPM.OUTPUT.CLOSELIVE Release optional video and hidden graphics resources.

if nargin < 2
    closeFigure = false;
end
if ~isempty(viz) && isfield(viz,'writer') && ~isempty(viz.writer)
    try
        close(viz.writer);
    catch exception
        fprintf(2,'IPM warning: could not close the live video: %s\n', ...
            exception.message);
    end
end
if ~isempty(viz) && isfield(viz,'figure') && ...
        isgraphics(viz.figure,'figure')
    hidden = strcmpi(viz.figure.Visible,'off');
    if closeFigure || hidden
        try
            delete(viz.figure);
        catch exception
            fprintf(2,'IPM warning: could not close the live figure: %s\n', ...
                exception.message);
        end
    end
end
end
