function res = compute_residual_images_dr_block_new(x, y, T, A, At, H, W, Wm, reduction_version)
%compute_residual_images_dr_block: compute the residual images for each 
% channel of interest.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x                       wideband image cube [N(1), N(2), L]
% > y                       visibilities (blocked) {L}{nblocks}
% > T                       pseudo singular values from the 
%                           reduction operator {L}{nblocks} (Sigma)
% > A                       measurement operator @[1]
% > At                      adjoint measurement operator @[1]
% > H                       holographic matrices G'*G {L}
% > W                       mask of non-zero columns of G {L}{nblocks}
%                           [] for the case W is absent
% > Wm                      Fourier masking operators (selection Fourier plane) {L}{nblocks}
% > reduction_version       option 1: embedding operator F\Phi^t
%                           option 2: embedding operator G^t
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

% Variable flag for the case where W is present
flagW = 0;
if ~isempty(W)
    flagW = 1;
end

% number of over-sampled pixels
if flagW
    No = size(W{1}{1}, 1);
else
    No = size(H{1}{1}, 2);
end

R = length(H{1});
[Ny, Nx, ci] = size(x);
res = zeros(size(x));

for i = 1 : ci
    Fx = A(x(:,:,i));
    g2 = zeros(No,1);
    for j = 1 : R
        if reduction_version == 1
            if flagW
                tmp = FT2(real(At(H{i}{j} * Fx(W{i}{j}))));
            else
                tmp = FT2(real(At(H{i}{j} * Fx)));
            end
            tmp = tmp(:);
            res_f{i}{j} = y{i}{j} - T{i}{j} .* tmp(Wm{i}{j});
            tmp = zeros(size(Wm{i}{j}));
            tmp(Wm{i}{j}) = T{i}{j} .* res_f{i}{j};
            if flagW
                g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            else
                g2 = g2 + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            end
        elseif reduction_version == 2
            if flagW
                res_f{i}{j} = y{i}{j} - T{i}{j} .* (H{i}{j} * Fx(W{i}{j}));
                g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j}' * (T{i}{j} .* res_f{i}{j});
            else
                res_f{i}{j} = y{i}{j} - T{i}{j} .* (H{i}{j} * Fx);
                g2 = g2 + H{i}{j}' * (T{i}{j} .* res_f{i}{j});
            end
        end
    end
    res(:,:,i) = real(At(g2));
end

end
