function y = dataReduce(x, Gt, W, At, Sigma, Mask)
% the embeddding operator R = \sigma * S * F * Phi^T
% Complex -> Complex
%
% Author: Ming Jiang, E-mail: ming.jiang@epfl.ch
%
FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));
x1 = Gt * x;
if ~isempty(W)
    tmp = zeros(size(W));
    tmp(W) = x1;
    x1 = tmp;
end
im = FT2(real(At(x1)));
im = im(:);
y = Sigma .* im(Mask);
