function [sig, max_psf, dirty_image] = ...
    compute_reweighting_lower_bound_sara(y, W, G, A, At, Ny, Nx, oy, ox, wavelet_basis, filters_length, nlevel)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% estimate mu: ratio between nuclear and l21 norm priors applied to
% the dirty image
dirty_image = zeros([Ny, Nx]);
temp = zeros(No, 1);
for b = 1:numel(G{1})
    temp(W{1}{b}) = temp(W{1}{b}) + G{1}{b}' * y{1}{b};
end
dirty_image(:,:) = At(temp);

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
[B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No);
B = B/max_psf;

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel); % suppose Dirac is the last entry
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
sig = std(abs(Psit_full(reshape(B, [Ny, Nx]))));

end
