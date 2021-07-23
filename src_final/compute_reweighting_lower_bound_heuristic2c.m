function [sig, sig_bar, mu, mu_bar, mu_c, sig_c, sig_w] = ...
    compute_reweighting_lower_bound_heuristic2c(Ny, Nx, ...
    nChannels, filters_length, nlevel, sigma_noise, ...
    algo_version, Qx, Qy, overlap_size, window_type, ...
    squared_operator_norm, alph)
 
N = Ny*Nx;    % number of image pixels

% compute number of wavelet coefficients
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel);
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)

% compute sig and sig_bar
% sig = sqrt(N*sum((sigma_noise.^2)./squared_operator_norm)/s);
sig_w = sqrt(mean((sigma_noise.^2)./squared_operator_norm)); %/numel(filters_length));
mu_c = sqrt(2)*gamma((nChannels+1)/2)/gamma(nChannels/2);
sig_c = sqrt(nChannels - mu_c^2); 
sig = sig_w*(mu_c + alph*sig_c);

% compute sig_bar
sig_bar = sig_w;
mu = sig;
mu_bar = sig_bar;
