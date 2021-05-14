function [sig, max_psf, mu, l11_norm, l11_norm_x0, dirty_image_precond] = ...
    compute_reweighting_lower_bound_sara(y, W, G, aW, A, At, Ny, Nx, oy, ox, ...
    wavelet_basis, filters_length, nlevel, sigma_noise, rw_type, x0, Anorm, squared_operator_norm)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% estimate mu: ratio between nuclear and l21 norm priors applied to
% the dirty image
temp = zeros(No, 1);
for b = 1:numel(G{1})
    temp(W{1}{b}) = temp(W{1}{b}) + G{1}{b}' * (aW{1}{b}.*y{1}{b});
end
dirty_image_precond = At(temp);
[B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No, sigma_noise, 1234);
% dirty_image_precond = dirty_image_precond/max_psf;
dirty_image_precond = dirty_image_precond/Anorm;

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
% B = B/max_psf;
B = B/squared_operator_norm;

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel); % suppose Dirac is the last entry
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
l11_norm = sum(abs(Psit_full(x0)));
l11_norm_x0 = sum(abs(Psit_full(dirty_image_precond)));

%! wavelet prior
sig = std(abs(Psit_full(reshape(B, [Ny, Nx]))));
if strcmp(rw_type, "ground_truth")
    mu = 1/l11_norm_x0;
else
    mu = 1/l11_norm;
end

% if strcmp(rw_type, "ground_truth")
%     mu = 1/sum(sig*log(d0/sig + 1));
% else
%     mu = 1/sum(sig*log(d1/sig + 1));
% end

end
