function [sig, mu0, mu, mu_bar, max_psf, l21_norm, nuclear_norm, l21_norm_x0, nuclear_norm_x0, dirty_image_precond] = ...
    compute_reweighting_lower_bound_test(y, W, G, aW, A, At, Ny, Nx, oy, ox, ...
    nChannels, wavelet_basis, filters_length, nlevel, sigma_noise, x0, Anorm, squared_operator_norm)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% form dirty image (no normalization)
dirty_image_precond = zeros([Ny, Nx, nChannels]);
for l = 1:nChannels
    temp = zeros(No, 1);
    for b = 1:numel(G{l})
        temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * (aW{l}{b}.*y{l}{b}); % sqrt(aW{i}{j}) .* sqrt(aW{i}{j}) .* y{l}{b}
    end
    dirty_image_precond(:,:,l) = At(temp);
end

% generate noise transferred from data to image domain, compute maximum of 
% the psf in each channel
[B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No, sigma_noise, 1234);

% evaluate nuclear norm of the dirty image -> regularization parameter
% dirty_image_precond = dirty_image_precond./reshape(max_psf, [1, 1, nChannels]);
dirty_image_precond = dirty_image_precond/Anorm; % Phi applied twice, and squared preconditioned operator norm for the normalisation
[~,S0,~] = svd(reshape(dirty_image_precond, [N, nChannels]),'econ');
nuclear_norm = sum(abs(diag(S0)));
[~,S0,~] = svd(reshape(x0, [N, nChannels]),'econ');
nuclear_norm_x0 = sum(abs(diag(S0)));

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel);
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)

% evaluate l21 norm of the dirty image -> regularization parameter
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
d1 = sqrt(sum(Psit_full(dirty_image_precond).^2, 2));
d0 = sqrt(sum(Psit_full(x0).^2, 2));
l21_norm = sum(d1); 
l21_norm_x0 = sum(d0); 

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
B = B/squared_operator_norm; %! normalize noise by the squared norm of the operator

sig = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));
mu0 = 1/l21_norm_x0;
mu = 1/l21_norm;
mu_bar = 1/nuclear_norm;

end