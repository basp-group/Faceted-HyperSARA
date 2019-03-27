function y = operatorR(x, Phi_t, Sigma, Mask)
% the embeddding operator R = \sigma * S * F * Phi^T
% Complex -> Complex
FT2 = @(x) fftshift(fft2(ifftshift(x)));
im = FT2(Phi_t(x));
im = im(:);
y = Sigma .* im(Mask);
