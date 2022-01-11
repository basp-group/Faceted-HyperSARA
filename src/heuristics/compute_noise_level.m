function [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
    compute_noise_level(Ny, Nx, n_channels, std_noise, algo_version, ...
    Qx, Qy, overlap_size, squared_operator_norm)
% [summary]
%
% [extended_summary]
%
% Parameters
% ----------
% Ny : [type]
%     [description]
% Nx : [type]
%     [description]
% n_channels : [type]
%     [description]
% std_noise : [type]
%     [description]
% algo_version : [type]
%     [description]
% Qy : [type]
%     [description]
% overlap_size : [type]
%     [description]
% squared_operator_norm : [type]
%     [description]
%
% Returns
% -------
% [type]
%     [description]

% compute sig and sig_bar
sig_sara = sqrt(mean((std_noise.^2) ./ squared_operator_norm));
mu_chi = sqrt(2) * gamma((n_channels + 1) / 2) / gamma(n_channels / 2);
sig_chi = sqrt(n_channels - mu_chi^2);
sig = sig_sara * (mu_chi + sig_chi);

% compute sig_bar
if strcmp(algo_version, 'hs')
    sig_bar = sqrt(Nx * Ny * sum(std_noise.^2 ./ squared_operator_norm) / min(Nx * Ny, n_channels));
else
    Q = Qx * Qy;
    rg_y = split_range(Qy, Ny);
    rg_x = split_range(Qx, Nx);
    I = zeros(Q, 2);
    dims = zeros(Q, 2);
    for qx = 1:Qx
        for qy = 1:Qy
            q = (qx - 1) * Qy + qy;
            I(q, :) = [rg_y(qy, 1) - 1, rg_x(qx, 1) - 1];
            dims(q, :) = [rg_y(qy, 2) - rg_y(qy, 1) + 1, rg_x(qx, 2) - rg_x(qx, 1) + 1];
        end
    end
    clear rg_y rg_x;

    rg_yo = split_range(Qy, Ny, overlap_size(1));
    rg_xo = split_range(Qx, Nx, overlap_size(2));
    Io = zeros(Q, 2);
    dims_o = zeros(Q, 2);
    for qx = 1:Qx
        for qy = 1:Qy
            q = (qx - 1) * Qy + qy;
            Io(q, :) = [rg_yo(qy, 1) - 1, rg_xo(qx, 1) - 1];
            dims_o(q, :) = [rg_yo(qy, 2) - rg_yo(qy, 1) + 1, rg_xo(qx, 2) - rg_xo(qx, 1) + 1];
        end
    end
    clear rg_yo rg_xo;

    sig_bar = zeros(Q, 1);
    for q = 1:Q
        Noq = prod(dims_o(q, :));
        sig_bar(q) = sqrt(Noq * sum(std_noise.^2 ./ squared_operator_norm) / min(Noq, n_channels));
    end
end
