function [B, max_psf] = create_dirty_noise_dr(y, A, At, H, T, W, Nx, Ny, No)

%! to be documented
%! create "dirty" noise matrix

nChannels = numel(y);
max_psf = zeros(nChannels, 1);
B = zeros(Nx*Ny, nChannels);
dirac = zeros(Ny, Nx);
dirac(floor([Ny, Nx]/2) + 1) = 1;
N = Nx*Ny;
AD = A(dirac);

for l = 1:nChannels
    b_l = zeros(No, 1);
    psf_l = zeros(No, 1);
    for b = 1:numel(y{l})
        noise = (randn(size(y{l}{b})) + 1i*randn(size(y{l}{b})))/sqrt(2);
        b_l(W{l}{b}) = b_l(W{l}{b}) + H{l}{b}' * (T{l}{b} .* noise);
        psf_l(W{l}{b}) = psf_l(W{l}{b}) + H{l}{b}' * (T{l}{b} .* (T{l}{b} .* (H{l}{b} * AD(W{l}{b}))));  
    end
    B(:,l) = reshape(At(b_l),[N,1]);
    max_psf(l) = max(reshape(At(psf_l),[N,1]));     
end

end
