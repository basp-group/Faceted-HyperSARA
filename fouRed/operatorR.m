function y = operatorR(x, Phi_t, Sigma, Mask)
% the embeddding operator R = \sigma * S * F * Phi^T = \sigma * S * F * Dr^T * Z^T * F^-1 * G^T

FT2 = @(x) fftshift(fft2(ifftshift(x)));
im = FT2(real(Phi_t(x)));
im = im(:);
y = Sigma .* im(Mask);
