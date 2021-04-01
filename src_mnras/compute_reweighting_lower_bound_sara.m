function [sig, max_psf, dirty_image] = ...
    compute_reweighting_lower_bound_sara(y, W, G, A, At, Ny, Nx, oy, ox, Psit)

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
sig = std(abs(Psit(B)));

end
