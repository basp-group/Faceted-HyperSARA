function [f, h] = display_real_images(x, fig_size, shift_colorbar, clim_log, map_img, fontsize)

    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    imagesc(x, clim_log);
    colormap(gca, map_img);
    
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    set(gca,'ColorScale','log')
    h = colorbar;
    set(h,'Location','southoutside');
    set(h,'Fontsize',fontsize);
    pause(0.5)
    h.Position = h.Position + shift_colorbar;
%     aspect = get(ax,'PlotBoxAspectRatio');
%     set(ax,'Units','pixels');
%     pos = get(ax,'Position');
%     pos(3) = aspect(1)/aspect(2)*pos(4);
%     set(ax,'Position',pos);
    % to be corrected!!
%     rectangle('Position',[cen-floor(len/2), len, len],'EdgeColor',[1 1 1],...
%     'LineWidth',2)
%     rectangle('Position',[cen0-floor(len0/2), len0, len0],'EdgeColor',[1 1 1],...
%     'LineWidth',2)
%     rectangle('Position',[0,0, len0, len0],'EdgeColor',[1 1 1],...
%     'LineWidth',2)
end