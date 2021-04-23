function [sig, max_psf, mu, l21_norm, l21_norm_x0, dirty_image] = ...
    compute_reweighting_lower_bound_sara(y, W, G, A, At, Ny, Nx, oy, ox, ...
    wavelet_basis, filters_length, nlevel, sigma_noise, rw_type, x0, Anorm)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% estimate mu: ratio between nuclear and l21 norm priors applied to
% the dirty image
temp = zeros(No, 1);
for b = 1:numel(G{1})
    temp(W{1}{b}) = temp(W{1}{b}) + G{1}{b}' * y{1}{b};
end
dirty_image = At(temp);
[B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No, sigma_noise, 1234);
% dirty_image = dirty_image/max_psf;
dirty_image = dirty_image/Anorm;

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
% B = B/max_psf;
B = B/Anorm;

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel); % suppose Dirac is the last entry
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
d0 = sqrt(sum(Psit_full(x0).^2, 2));
d1 = sqrt(sum(Psit_full(dirty_image).^2, 2));
l21_norm = sum(d1); 
l21_norm_x0 = sum(d0); 

%! wavelet prior
sig = std(abs(Psit_full(reshape(B, [Ny, Nx]))));
if strcmp(rw_type, "ground_truth")
    mu = 1/sum(sig*log(d0/sig + 1));
else
    mu = 1/sum(sig*log(d1/sig + 1));
end

end
