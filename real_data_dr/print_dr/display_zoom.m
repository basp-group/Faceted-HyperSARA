function [f, h] = display_zoom(x, fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel)

    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im = imagesc(x, clim_log);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    set(gca,'ColorScale','log');
    set(h,'XTick',tickLabel);
    pause(0.5)
    % Get the new aspect ratio data
    aspect = get(ax,'PlotBoxAspectRatio');
    % Change axes Units property (this only works with non-normalized units)
    set(ax,'Units','pixels');
    % Get and update the axes width to match the plot aspect ratio.
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    pause(0.1)
    set(h,'Units', 'pixels')
    pos_h = get(h,'Position');
    pos_h = [pos(1)+pos(3)+pixel_shift_colorbar, pos_h(2), 2*pos_h(3), pos(4)];
    set(h,'Position',pos_h);
    pause(0.1)
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    pause(0.5)
end