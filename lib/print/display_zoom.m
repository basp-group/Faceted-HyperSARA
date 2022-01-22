function [f, h] = display_zoom(x, fig_size, pixel_shift_colorbar, clim, cmap, fontsize, tickLabel)
% Display a the input image :math:`x` (zoomed-in regions from the MNRAS
% paper).
%
% Parameters
% ----------
% x : array
%     Input image.
% fig_size : array (int)
%     Size of the figure (in pixels).
% pixel_shift_colorbar : array
%     Shift applied to the position of the colorbar (in pixels).
% clim : array
%     Values associated to the limits of the colorbar.
% cmap : array
%     Colormap.
% fontsize : int
%     Size of the font.
% tickLabel : array
%     Values for the ticks reported on the colorbar.
%
% Returns
% -------
% f : handle
%     Figure handle.
% h : handle
%     Colorbar handle.
%

%%
    f = figure('visible', 'on');
    % sets axes background
    set(gca, 'Color', 'none');
    set(f, 'PaperUnits', 'centimeters');
    set(f, 'PaperType', 'A4');
    set(f, 'PaperOrientation', orient);
    set(f, 'units', 'pixel', 'outerposition', [0 0 fig_size]);
    im = imagesc(x, clim);
    colormap(gca, cmap);
    axis image;
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter = 'latex';
    h = colorbar;
    set(h, 'Fontsize', fontsize);
    set(h, 'color', 'white');
    set(gca, 'ColorScale', 'log');
    set(h, 'XTick', tickLabel);
    pause(0.5);
    % Get the new aspect ratio data
    aspect = get(ax, 'PlotBoxAspectRatio');
    % Change axes Units property (this only works with non-normalized units)
    set(ax, 'Units', 'pixels');
    % Get and update the axes width to match the plot aspect ratio.
    pos = get(ax, 'Position');
    pos(3) = aspect(1) / aspect(2) * pos(4);
    set(ax, 'Position', pos);
    pause(0.1);
    set(h, 'Units', 'pixels');
    pos_h = get(h, 'Position');
    pos_h = [pos(1) + pos(3) + pixel_shift_colorbar, pos_h(2), 2 * pos_h(3), pos(4)];
    set(h, 'Position', pos_h);
    set(ax, 'Units', 'normalized');
    pause(0.1);
    pos = get(ax, 'Position');
    annotation(f, 'rectangle', ...
               [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)], ...
               'Color', [1 1 1], ...
               'LineWidth', 4);

end
