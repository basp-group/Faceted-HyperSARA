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

% TODO: do it directly in parallel with faceted SARA?
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel);
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);

%! the number of elements reported below is only valid for the periodic
%('per') boundary condition (implicitly assumed to be used in HS_forward_sparsity)
% TODO: enable other different boundary conditions
l21_norm = sum(sqrt(sum(Psit_full(dirty_image).^2, 2))); 

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
B = zeros(N, nChannels);
dirac = zeros(Ny, Nx);
dirac(floor([Ny, Nx]/2) + 1) = 1;
AD = A(dirac);
max_psf = zeros(nChannels, 1);
for l = 1:nChannels
    temp = zeros(No, 1);
    z = zeros(No, 1);
    for b = 1:numel(y{l})
        noise = (randn(size(y{l}{b})) + 1i*randn(size(y{l}{b})))/sqrt(2);
        temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * noise;
        z(W{l}{b}) = z(W{l}{b}) + G{l}{b}' * (G{l}{b} * AD(W{l}{b}));  
    end
    B(:,l) = reshape(At(temp),[N,1]);
    max_psf(l) = max(reshape(At(z),[N,1]));     
end
B = B./reshape(max_psf, [1, nChannels]);
[~,S0,~] = svd(B,'econ');
sig = std(diag(S0));
sig_bar = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));

end
