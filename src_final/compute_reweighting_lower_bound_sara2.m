function [sig, mu] = ...
    compute_reweighting_lower_bound_sara2(Ny, Nx, sigma_noise, squared_operator_norm)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
% No = N*oy*ox; % size of oversampled Fourier space

% set-up global SARA dictionary
% dwtmode('zpd')
% [~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
% [~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel); % suppose Dirac is the last entry
% s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
% Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);

% compute sig (estimate of the "noise level" in SARA space) involved in the
% reweighting scheme

% sig = sqrt(N*(sigma_noise^2)/(squared_operator_norm*s));
sig = sqrt((sigma_noise^2)/squared_operator_norm);
mu = sig;
