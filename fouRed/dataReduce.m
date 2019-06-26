function y = dataReduce(x, Gt, At, No, W, Sigma, Mask)
% the embeddding operator R = \sigma * S * F * Phi^T
% Complex -> Complex
FT2 = @(x) fftshift(fft2(ifftshift(x)));
x1 = Gt * x;
if exist('No', 'var') && exist('W', 'var')
    tmp = zeros(No(1) * No(2), 1);
    tmp(W) = x1;
    x1 = tmp;
end
im = FT2(real(At(x1)));
im = im(:);
y = Sigma .* im(Mask);
