function [f, h] = display_image(x, fig_size, shift_colorbar, clim, map_img, fontsize, bool_log, location)

    f=figure('visible','on');
    % sets axes background
    set(gca, 'Color', 'none');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    
    % t = tiledlayout(1,1,'Padding','none'); % 'none' before R2021, 'tight'after
    % t.Units = 'inches';
    % t.OuterPosition = [0.25 0.25 4 4];
    % nexttile;
    
    imagesc(x, clim);
    colormap(gca, map_img);

    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    if bool_log
        set(ax,'ColorScale','log'); 
    end
    h = colorbar;
    set(h,'Location', location); % 'southoutside'
    set(h,'Fontsize',fontsize);
    pause(0.5)
    h.Position = h.Position + shift_colorbar;
