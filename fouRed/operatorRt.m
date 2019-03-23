function y = operatorRt(x, imsize, Phi, Sigma, Mask)

IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

Ny = imsize(1);
Nx = imsize(2);

y = zeros(Ny*Nx, 1);
x1 = x(:) .* Sigma;
y(Mask) = x1;
y = Phi(real(IFT2(reshape(y, Ny, Nx))));
