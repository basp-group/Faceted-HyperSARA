function display_image(x, clim, map_img, fontsize, log)

f=figure('visible','on');
% sets axes background
set(gca, 'Color', 'none');
set(f,'PaperUnits','centimeters')
set(f,'PaperType','A4');
set(f,'PaperOrientation',orient);
set(f,'units','pixel','outerposition',[0 0 600 600])
imagesc(x, clim);
colormap(gca, map_img);
if log
    set(gca,'ColorScale','log'); 
end
axis image
ax = gca;
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';
h = colorbar;
set(h,'Fontsize',fontsize)
