function [sig, sig_bar, max_psf, l21_norm, nuclear_norm, dirty_image] = ...
    compute_reweighting_lower_bound(y, W, G, A, At, Ny, Nx, oy, ox, ...
    nChannels, wavelet_basis, filters_length, nlevel)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% estimate mu: ratio between nuclear and l21 norm priors applied to
% the dirty image
dirty_image = zeros([Ny, Nx, nChannels]);
for l = 1:nChannels
    temp = zeros(No, 1);
    for b = 1:numel(G{l})
        temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * y{l}{b};
    end
    dirty_image(:,:,l) = At(temp);
end
[~,S0,~] = svd(reshape(dirty_image, [N, nChannels]),'econ');
nuclear_norm = sum(abs(diag(S0)));

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel);
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
l21_norm = sum(sqrt(sum(Psit_full(dirty_image).^2, 2))); 

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
[B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No);
B = B./reshape(max_psf, [1, nChannels]);
[~,S0,~] = svd(B,'econ');
sig_bar = std(diag(S0));
sig = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));

end
