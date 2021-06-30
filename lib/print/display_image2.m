function [f, h] = display_image2(x, clim, map_img, fontsize, log)

    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 600 600])
    imagesc(x, clim);
    colormap(gca, map_img);

    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    if log
        set(ax,'ColorScale','log'); 
    end
    h = colorbar;
    set(h,'Location','southoutside');
    set(h,'Fontsize',fontsize);
    pause(0.5)
    h.Position = h.Position + shift_colorbar;
