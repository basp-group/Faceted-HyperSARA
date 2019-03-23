function [rPhi, rPhit, rsPhi, rsPhit] = oper_fourierReduction(Ipsf, Sigma, Mask, imsize)
% Final reduction operators Phi_sing and (Phi_sing)^T
    rPhi = @(x) fftshift(fft2(ifftshift(Ipsf(x))));  % R Phi = F * Phi^T * Phi = F Ipsf: image -> vect
    rPhit = @(x) Ipsf(fftshift(ifft2(ifftshift(reshape(full(x), imsize)))));  % Phi^T R^T = Phi^T * Phi * F^T = Ipsf * F^T: vect -> image
    
    rsPhi = @(x) operator_rsPhi(x, rPhi, Sigma, Mask, imsize);              % \Sigma * S * C: vect/image -> vect
    rsPhit = @(x) operator_rsPhit(x, rPhit, Sigma, Mask, imsize);           % C^T * S' * \Sigma: vect/image -> vect
    
end

function y = operator_rsPhi(x, rPhi, Sigma, Mask, imsize)
    tmp = rPhi(reshape(x,imsize));
    y = Sigma.*tmp(Mask);
end

function y = operator_rsPhit(x, rPhit, Sigma, Mask, imsize)
    tmp = zeros(prod(imsize), 1);
    tmp(Mask) = Sigma.*x(:);
    y = rPhit(tmp);
end
