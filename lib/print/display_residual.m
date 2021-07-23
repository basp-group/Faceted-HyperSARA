function [f, h] = display_residual(x, fig_size, clim_log, map_img, fontsize, bool_real)

    f=figure('visible','on');
    % sets axes background
    set(gca, 'Color', 'none');
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    imagesc(x, clim_log);
    % imagesc(x);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    if bool_real
        set(h,'color','white')
    end
end
