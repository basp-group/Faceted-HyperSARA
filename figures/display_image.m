function display_image(x, clim, map_img, fontsize)

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
h = colorbar;
set(h,'Fontsize',fontsize)
