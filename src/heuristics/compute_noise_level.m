function [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
    compute_noise_level(Ny, Nx, n_channels, std_noise, algo_version, ...
    Qx, Qy, overlap_size, squared_operator_norm)
% Compute heuristic order of magnitude for the regularization parameters
% involved in SARA :cite:p:`Abdulaziz2019` or Faceted HyperSARA
% :cite:p:`Thouvenin2021`.
%
% Parameters
% ----------
% Ny : int
%     Spatial image size along axis y.
% Nx : int
%     Spatial image size along axis y.
% n_channels : int
%     Number of frequency channels.
% std_noise : double
%     Estimate of the standard deviation of the white Gaussian noise
%     affecting the data.
% algo_version : string
%     Imaging problem considered (HyperSARA, ``'hs'`` or Faceted
%     HyperSARA ``'fhs'``),
% Qx : int
%     Number of spatial facets along spatial axis x.
% Qy : int
%     Number of spatial facets along spatial axis y.
% overlap_size : int[2]
%     Overlap size between consecutive facets along each axis (y and x).
% squared_operator_norm : double
%     Square of the measurement operator norm, :math:`\|\Phi\|_2^2`.
%
% Returns
% -------
% sig : double
%     Heuristic value sparsity regularization parameter.
% sig_bar : double
%     Heursitic value low-rankness regularization parameter.
% mu_chi : double
%     Heuristic estimate of the mean of the noise using a
%     :math:`\chi_2` approximation. (?)
% sig_chi : double
%     Heuristic estimate of the standard deviation of the noise using a
%     :math:`\chi_2` approximation. (?)
% sig_sara : double
%     Heuristic estimate of the noise level (standard deviation) when
%     transferred to the SARA dictionary domain. (?)
%

% compute sig and sig_bar
sig_sara = full(sqrt(mean((std_noise.^2) ./ squared_operator_norm)));
mu_chi = full(sqrt(2) * gamma((n_channels + 1) / 2) / gamma(n_channels / 2));
sig_chi = full(sqrt(n_channels - mu_chi^2));
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
        sig_bar(q) = full(sqrt(Noq * sum(std_noise.^2 ./ squared_operator_norm) / min(Noq, n_channels)));
    end
end
