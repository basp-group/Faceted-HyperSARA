function res = compute_residual_images_dr_block(x, y, T, A, At, H, W)
%compute_residual_images_dr_block: compute the residual images for each 
% channel of interest.
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
% [../../2019]
%-------------------------------------------------------------------------%
%%

FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

[Ny, Nx, ci] = size(x);
res = zeros(size(x));

for i = 1 : ci
    g2 = zeros(Ny*Nx,1);
    Fx0 = A(real(x(:,:,i)));
    for j = 1:length(H{i})
        dimH = sqrt(numel(H{i}{j}));
        Fx = FT2(real(At(reshape(H{i}{j}, dimH, dimH) * Fx0)));
        Fx = Fx(:);
        res_f = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
        u2 = T{i}{j} .* res_f;
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
   
        res(:,:,i) = res(:,:,i) + real(At(reshape(H{i}{j}, dimH, dimH)*A(real(IFT2(reshape(g2, Ny, Nx))))));
    end
end

end
