function res = compute_residual_images_dr(x, y, T, A, At, H, W)
FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

[Ny, Nx, ci] = size(x);
res = zeros(size(x));

for i = 1 : ci
    Fx = FT2(real(At(H{i} * A(real(x(:,:,i))))));    % adapt
    Fx = Fx(:);
    g2 = zeros(numel(Fx),1);
    for j = 1 : length(T{i})
        res_f = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
        u2 = T{i}{j} .* res_f;
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
    end
    res(:,:,i) = real(At(H{i}*A(real(IFT2(reshape(g2, Ny, Nx)))))); % adapt
end

end
