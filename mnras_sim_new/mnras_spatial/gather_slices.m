clc; clear all; close all

% fileroot = @(ind) strcat('x_cygASband_Cube_1024_2048_20_sara_srf=2_none_Qy=1_Qx=1_Qc=20_ind=', ...
%     num2str(ind), ...
%     '_gam=1_homotopy=0_rwtype=dirty_updatereg=_regtype=inv_snr=40.fits');

% fileroot = @(ind) fullfile('results/cygASband_Cube_1024_2048_20_spatial/sara', strcat('spatial_cygASband_Cube_1024_2048_20_sara_none_srf=2_Ny=1024_Nx=2048_L=20_Qy=1_Qx=1_Qc=20_ind=', num2str(ind),'_g=1_gb=1_overlap=0_0_hom=0_rwt=heuristic_updreg=0_regtype=heuristic_snr=40_rw=25.mat'));

% N = [1024, 2048];
% nChannels = 20;

% x = zeros([N, nChannels]);
% res = zeros([N, nChannels]);
% for l = [1:4, 6, 8:nChannels]
%     m = matfile(fileroot(l));
%     % x_ = fitsread(fileroot(l));
%     x(:,:,l) = m.xsol;
%     res(:,:,l) = m.res;
% end

% fitswrite(x, "x_sara_heuristic.fits")
% fitswrite(res, "res_sara_heuristic.fits")

%%
% x0 = fitsread('../data/cygASband_Cube_1024_2048_20.fits');
% asnr = mean(10*log10(sum(x0.^2, [1,2])./sum((x0 - xsara).^2, [1,2])));

%% HS

filename = fullfile('results/cygASband_Cube_1024_2048_20_spatial/hypersara/heuristic', 'spatial_cygASband_Cube_1024_2048_20_hypersara_none_srf=2_Ny=1024_Nx=2048_L=20_Qy=1_Qx=1_Qc=1_ind=1_g=1_gb=1_overlap=0_0_hom=0_rwt=heuristic_updreg=0_regtype=heuristic_snr=40_rw=30.mat');

load(filename, 'xsol', 'res')
fitswrite(xsol, "x_fhs_heuristic.fits")
fitswrite(res, "res_fhs_heuristic.fits")

%% FHS

filename = fullfile('results/cygASband_Cube_1024_2048_20_spatial/cw/heuristic', 'spatial_cygASband_Cube_1024_2048_20_cw_triangular_srf=2_Ny=1024_Nx=2048_L=20_Qy=4_Qx=4_Qc=1_ind=1_g=1_gb=1_overlap=0.5_0.5_hom=0_rwt=heuristic_updreg=0_regtype=heuristic_snr=40_rw=30.mat');

load(filename, 'xsol', 'res')
fitswrite(xsol, "x_hs_heuristic.fits")
fitswrite(res, "res_hs_heuristic.fits")