function y = operatorRPhi(x, Ipsf, Sigma, Mask, imsize)
% New reduced measurement operator: Sigma * S * F * Ipsf
% Real -> Complex
FT2 = @(x) fftshift(fft2(ifftshift(x)));
Ny = imsize(1);
Nx = imsize(2);
tmp = FT2(Ipsf(reshape(x, Ny, Nx)));
y = Sigma .* tmp(Mask);
