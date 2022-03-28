% Main script to generate the reference synthetic wideband image cubes
% used in :cite:p:`Thouvenin2021`. Requires the
% `S_DDE_MODEL.fits <https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits>`_
% image file associated with :cite:p:`Dabbech2021`, which can be retrieved
% as follows from a unix (MAC, Linux, Windows WSL or MinGw) terminal.
%
% .. code-block:: bash
%
%    # if on MAC:
%    # brew install wget
%    cd path/to/Faceted-HyperSARA
%    mkdir data && cd data
%    wget -P . https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits
%
% .. note::
%
%    Cubes and additional informations for the synthetic data experiments 
%    are saved by default in the folder ``Faceted-HyperSARA/data/``.
%

clear;
clc;
close all;

data_path = strcat('..', filesep, 'data', filesep);
mkdir(data_path);
addpath(data_path);
addpath(strcat('..', filesep, 'lib', filesep, 'generate_data'));

%% Freq info
nu_1 = 2.052e9;  % starting freq
dnu = 16e6;  % freq step
L = 100;  % number of channels
nu_vect = [nu_1 (dnu * (1:L - 1) + nu_1)];

%% for info
dfreqSamples = 5;
nu_2k = nu_vect(1:dfreqSamples:end);
nu_512 = nu_vect;

%% Cyg A images at different resolutions
% Note: dataset to be retrived at
% https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits,
% DOI: 10.17861/529cdcbc-7c18-47a6-970f-755a5da19071
im_ref = fitsread(fullfile(data_path, 'S_DDE_MODEL.fits'));
% adjust the FoV desired and the image size 2kx1k
Nh = 3000;
imr = imresize(im_ref, [Nh, Nh]);
im_2k = imr(Nh / 2 + 1 - 512:Nh / 2 + 1 + 511, Nh / 2 + 50 - 1024:Nh / 2 + 50 + 1023);
im_2k = im_2k .* (im_2k > 2e-5);
fitswrite(im_2k, [data_path, 'cygASband_', num2str(size(im_2k, 1)), '_', num2str(size(im_2k, 2)), '.fits']);
% image  down size 512x256
im_512 = imresize(im_2k, [256, 512]);
im_512 = im_512 .* (im_512 > 8e-5);
fitswrite(im_512, [data_path, 'cygASband_', num2str(size(im_512, 1)), '_', num2str(size(im_512, 2)), '.fits']);

%% generate spectral maps
im_512_init = im_512;
im_2k_init = im_2k;
fact_dim = size(im_2k, 1) ./ [size(im_2k, 1) size(im_512, 1)];

% spectral index map: 1st component:smooth image and then take the log10
N_ = [64, 64];
[X, Y] = meshgrid(-N_(1) / 2:N_(1) / 2 - 1, -N_(1) / 2:N_(2) / 2  - 1);
h = [X(:) Y(:)];
theta = 0;

% 2kx1k
sig_x0 = 10;
sig_y0 = 10;
rho = gauss2d(sig_x0, sig_y0, theta);
beam_2k = reshape(rho(h), N_(1), N_(2));
beam_2k_nz = beam_2k ./ sum(sum(reshape(rho(h), N_(1), N_(2))));
im_2kSLog = log10(conv2(im_2k, beam_2k_nz, 'same') .* (im_2k_init > 0));

% 512x256
sig_x = sig_x0 / fact_dim(2);
sig_y = sig_y0 / fact_dim(2);
rho = gauss2d(sig_x, sig_y, theta);
beam_512 = reshape(rho(h), N_(1), N_(2));
beam_512_nz = beam_512 ./ sum(sum(beam_512));
im_512SLog = log10(conv2(im_512, beam_512_nz, 'same') .* (im_512_init > 0));

% spectral index map: 2nd component : gaussian process
sig_gprocess = 10; %% sigma value for Gaussian process

% 2kx1k
N2k = size(im_2k);
theta = max(min(randn(1, 1), 1), -1) * pi;
sig_x = sig_gprocess / fact_dim(1);
sig_y = sig_gprocess / fact_dim(1);
rho = gauss2d(sig_x, sig_y, theta);

[f1Index_2k, ~] = stationary_Gaussian_process(N2k(1), N2k(2), rho);
f1Index_2k = f1Index_2k ./ max(max(abs(f1Index_2k)))  .* (im_2k_init > 0);

% 512x256
N512 = size(im_512);
sig_x = sig_gprocess / fact_dim(2);
sig_y = sig_gprocess / fact_dim(2);
rho  = gauss2d(sig_x, sig_y, theta);
[f1Index_512, ~] = stationary_Gaussian_process(N512(1), N512(2), rho);
f1Index_512  = f1Index_512 ./ max(max(abs(f1Index_512)))  .* (im_512_init > 0);

% spectral index map: adjusting the range
spectralIds2k = (f1Index_2k .* (im_2k_init > 0) + im_2kSLog);
spectralIds512 = (f1Index_512 .* (im_512_init > 0) + im_512SLog);

spectralIds2k = (spectralIds2k - max(max(spectralIds2k)) - 0.3);
spectralIds512 = (spectralIds512 - max(max(spectralIds512)) - 0.3);

spectralIds2k(spectralIds2k < -5) = NaN;
spectralIds512(spectralIds512 < -5) = NaN;

%% generate curvature map
scale = 0.5;
theta = max(min(randn(1, 1), 1), -1) * pi;
% 2kx1k
sig_x = sig_gprocess / fact_dim(1);
sig_y = sig_gprocess / fact_dim(1);
rho  = gauss2d(sig_x, sig_y, theta);
[f1Index_2k, ~] = stationary_Gaussian_process(N2k(1), N2k(2), rho);
Curv2k = scale * f1Index_2k ./ max(max(abs(f1Index_2k))) .* (im_2k_init > 0);
Curv2k(Curv2k < -scale) = NaN;

% 512x256
sig_x = sig_gprocess / fact_dim(2);
sig_y = sig_gprocess / fact_dim(2);
rho  = gauss2d(sig_x, sig_y, theta);
[f1Index_512, ~] = stationary_Gaussian_process(N512(1), N512(2), rho);
Curv512 = scale * f1Index_512 ./ max(max(abs(f1Index_512))) .* (im_512_init > 0);
Curv512(Curv512 < -scale) = NaN;

% figure,subplot 211, imagesc(Curv2k),axis image,colorbar,title('Curv')
% subplot 212, imagesc(Curv512),axis image,colorbar,title('Curv')
%% Build cube

cube2k = zeros(N2k(1), N2k(2), L);
cube512 = zeros(N512(1), N512(2), L);

spectralIds2k(isnan(spectralIds2k)) = 0;
spectralIds512(isnan(spectralIds512)) = 0;

Curv2k(isnan(Curv2k)) = 0;
Curv512(isnan(Curv512)) = 0;

for i = 1:numel(nu_vect)
    cube2k(:, :, i) = im_2k_init .* (nu_vect(i) ./ nu_1).^(spectralIds2k + Curv2k * log(nu_vect(i) ./ nu_1));
    cube512(:, :, i) = im_512_init .* (nu_vect(i) ./ nu_1).^(spectralIds512 + Curv512 * log(nu_vect(i) ./ nu_1));
end
cube2k = cube2k(:, :, 1:dfreqSamples:end);

fitswrite(cube2k, [data_path, 'cygASband_Cube_', num2str(size(cube2k, 1)), '_' num2str(size(cube2k, 2)), '_', num2str(size(cube2k, 3)), '.fits']);
fitswrite(cube512, [data_path, 'cygASband_Cube_', num2str(size(cube512, 1)), '_' num2str(size(cube512, 2)), '_', num2str(size(cube512, 3)), '.fits']);

fitswrite(spectralIds2k, [data_path, 'SpectralIdxMap', num2str(size(cube2k, 1)), '_', num2str(size(cube2k, 2)), '_', num2str(size(cube2k, 3)), '.fits']);
fitswrite(spectralIds512, [data_path, 'SpectralIdxMap', num2str(size(cube512, 1)), '_', num2str(size(cube512, 2)), '_', num2str(size(cube512, 3)), '.fits']);

fitswrite(Curv2k, [data_path, 'CurvMap', num2str(size(cube2k, 1)), '_', num2str(size(cube2k, 2)), '_', num2str(size(cube2k, 3)), '.fits']);
fitswrite(Curv512, [data_path, 'CurvMap', num2str(size(cube512, 1)), '_', num2str(size(cube512, 2)), '_', num2str(size(cube512, 3)), '.fits']);

function rho = gauss2d(sig_x, sig_y, theta)

a = cos(theta).^2 / (2 * sig_x^2) + sin(theta).^2 / (2 * sig_y^2);
b = sin(2 * theta).^2 / (4 * sig_x^2) + sin(2 * theta).^2 / (4 * sig_y^2);
c = sin(theta).^2 / (2 * sig_x^2) + cos(theta).^2 / (2 * sig_y^2);
rho = @(h) exp(-(a .* h(:, 1).^2 + 2 * b .* h(:, 1) .* h(:, 2) + c * h(:, 2).^2));

end
