function [x0_new, X0_new, f_new, c_new] = generate_undersampled_cube(x0, f, Ny, Nx, c, unds)
% Function returning a spectrally undersampled version of an input
% wideband image cube.
%
% Parameters
% ----------
% x0 : double[:, :, :]
%     Synthetic wideband image (3D) [Ny, Nx, c].
% f : double[:]
%     Frequencies over which the wideband image `x0` is defined.
% Ny : int
%     Spatial dimension of `x0` along axis y.
% Nx : int
%     Spatial dimension of `x0` along axis y.
% c : int
%     Spectrsl dimension of `x0`.
% unds : int
%     Spectral downsampling factor.
% 
% Returns
% -------
% x0_new : double[:, :, :]
%     Output wideband image.
% X0_new : double[:, :, :]
%     Matrix-format version of `x0_new`.
% f_new : double[:]
%     Frequencies over which `x0_new` is defined.
% c_new : int
%     Number of spectral channels after spectral downsmapling.
%

    % Subsample the maps and reduce to desired size
    x0_new = zeros(Ny, Nx, c / unds);
    X0_new = zeros(Ny * Nx, c / unds);
    f_new = zeros(c / unds, 1);

    counter = 1;
    for i = round(linspace(1, c, c / unds))
        x0_new(:, :, counter) = imresize(x0(:, :, i), [Ny, Nx], 'nearest');
        a = x0_new(:, :, counter);
        X0_new(:, counter) = a(:);
        f_new(counter) = f(i);
        counter = counter + 1;
    end
    c_new = counter - 1;
