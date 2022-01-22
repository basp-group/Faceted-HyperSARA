function [f, h] = display_image(x, fig_size, shift_colorbar, clim, map_img, fontsize, bool_log, location)
% Display the input image :math:`x`.
%
% Parameters
% ----------
% x : array
%     Input image.
% fig_size : array (int)
%     Size of the figure (in pixels).
% shift_colorbar : array
%     Shift applied to the position of the colorbar.
% clim : array
%     Values associated to the limits of the colorbar.
% map_img : array
%     Colormap.
% fontsize : int
%     Size of the font.
% bool_log : bool
%     Flag to activate log scale for the displayed image.
% location : string
%     String specifynig the location of the colorbar within the figure.
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
    set(f, 'PaperOrientation', orient);
    set(f, 'units', 'pixel', 'outerposition', [0 0 fig_size]);

    % t = tiledlayout(1,1,'Padding','none'); % 'none' before R2021, 'tight'after
    % t.Units = 'inches';
    % t.OuterPosition = [0.25 0.25 4 4];
    % nexttile;

    imagesc(x, clim);
    colormap(gca, map_img);

    axis image;
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter = 'latex';
    if bool_log
        set(ax, 'ColorScale', 'log');
    end
    h = colorbar;
    set(h, 'Location', location); % 'southoutside'
    set(h, 'Fontsize', fontsize);
    pause(0.5);
    h.Position = h.Position + shift_colorbar;

end
