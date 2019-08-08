function res = compute_residual_images_dr(x, y, T, A, At, H, W)
%compute_residual_images_dr: compute the residual images for each channel 
% of interest. (no data blocking)
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x             wideband image cube [N(1), N(2), L]
% > y             visibilities (blocked) {L}{nblocks}
% > T             pseudo singular values from the 
%                 reduction operator {L}{nblocks} (Sigma)
% > A             measurement operator @[1]
% > At            adjoint measurement operator @[1]
% > H             holographic matrices G'*G {L}
% > W             Fourier masking operators (selection Fourier plane) {L}{nblocks}
%
% Output:
%
% < res           residual image cube [N(1), N(2), L]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

[Ny, Nx, ci] = size(x);
res = zeros(size(x));

for i = 1 : ci
    Fx = FT2(real(At(H{i} * A(real(x(:,:,i))))));
    Fx = Fx(:);
    g2 = zeros(numel(Fx),1);
    for j = 1 : length(T{i})
        res_f = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
        u2 = T{i}{j} .* res_f;
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
    end
    res(:,:,i) = real(At(H{i}*A(real(IFT2(reshape(g2, Ny, Nx))))));
end

end
