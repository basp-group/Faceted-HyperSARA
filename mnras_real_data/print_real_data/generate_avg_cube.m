clc; 
format compact

%% generate average cube

% % 16 sub-cubes, composed of 30 channels each
% % average 16 consecutive channels for compraison with DR
% ref_name = '/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data/xsol_FacetedHyperSARA_4-8GHz_NEW.fits';
% number_effective_channels = 30;
% step = 16;
% info        = fitsinfo(ref_name);
% rowend      = info.PrimaryData.Size(1);
% colend      = info.PrimaryData.Size(2);
% nchannels   = info.PrimaryData.Size(3);
% xl1 = zeros(rowend, colend, 3);
% x_avg = zeros(rowend, colend,number_effective_channels);
% for n = 1:number_effective_channels
%     x = fitsread(ref_name,'primary',...
%     'Info', info,...
%     'PixelRegion',{[1 1 rowend], [1 1 colend], [(n-1)*step+1 1 n*step]});
%     x_avg(:,:,n) = mean(x, 3);
% end
% fitswrite(x_avg,'xsol_FacetedHyperSARA_average.fits');
% 
% filename = 'xsol_FacetedHyperSARA_average.fits';
% x = fitsread(filename);
% xhs = zeros(size(x,1),size(x,2),3);
% xhs(:,:,1:2) = x(:,:,[1,end]);
% xhs(:,:,3) = mean(x, 3);
% fitswrite(xhs, 'x_fhs_avg_reduced.fits')

%%

% 16 sub-cubes, composed of 30 channels each
% average 16 consecutive channels for compraison with DR
ref_name = '/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data/res_FacetedHyperSARA_4-8GHz_NEW.fits';
number_effective_channels = 30;
step = 16;
info        = fitsinfo(ref_name);
rowend      = info.PrimaryData.Size(1);
colend      = info.PrimaryData.Size(2);
nchannels   = info.PrimaryData.Size(3);
xl1 = zeros(rowend, colend, 3);
x_avg = zeros(rowend, colend,number_effective_channels);
for n = 1:number_effective_channels
    x = fitsread(ref_name,'primary',...
    'Info', info,...
    'PixelRegion',{[1 1 rowend], [1 1 colend], [(n-1)*step+1 1 n*step]});
    x_avg(:,:,n) = mean(x, 3);
end
fitswrite(x_avg,'res_FacetedHyperSARA_average.fits');

filename = 'res_FacetedHyperSARA_average.fits';
x = fitsread(filename);
xhs = zeros(size(x,1),size(x,2),3);
xhs(:,:,1:2) = x(:,:,[1,end]);
xhs(:,:,3) = mean(x, 3);
fitswrite(xhs, 'r_fhs_avg_reduced.fits')