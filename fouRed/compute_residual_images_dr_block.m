function res = compute_residual_images_dr_block(x, y, T, A, At, H, W)
%compute_residual_images_dr: compute the residual images for each channel 
% of interest.
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
% > W             masking operators (selection Fourier plane) {L}{nblocks}
%
% Output:
%
% < residual_image  residual image cube [N(1), N(2), L]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

[Ny, Nx, ci] = size(x);
res = zeros(size(x));

for i = 1 : ci
    g2 = zeros(Ny*Nx,1);
    H_mat = sparse(size(H{i}{1}));
    for j = 1:length(H{i})
        Fx = FT2(real(At(H{i}{j} * A(real(x(:,:,i))))));
        Fx = Fx(:);
        res_f = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
        u2 = T{i}{j} .* res_f;
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
        H_mat = H{i}{j} + H_mat;
    end
    res(:,:,i) = real(At(H_mat*A(real(IFT2(reshape(g2, Ny, Nx)))))); % adapt
end

end
