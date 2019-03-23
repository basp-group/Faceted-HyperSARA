function y = operatorRPhit(x, Ipsf, Sigma, Mask, imsize)
% Adjoint of the new reduced measurement operator: Ipsf * F^T * S * Sigma
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));
Ny = imsize(1);
Nx = imsize(2);
x1 = zeros(Ny * Nx, 1);
x1(Mask) = Sigma .* x(:);
x1 = reshape(x1, Ny, Nx);
y = Ipsf(IFT2(x1));
