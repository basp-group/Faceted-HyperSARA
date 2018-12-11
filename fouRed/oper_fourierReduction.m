function [rPhi, rPhit, rsPhi, rsPhit] = oper_fourierReduction(R, Rt, Sigma, Mask, Gw, A, At, imsize)
% Final reduction operators Phi_sing and (Phi_sing)^T
    rPhi = @(x) R(Gw * A(x));  % R Phi = F * Phi^T * Phi: image -> vect
    rPhit = @(x) real(At(Gw' * Rt(x)));  % Phi^T R^T = Phi^T * Phi * F^T: vect -> image
    
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
