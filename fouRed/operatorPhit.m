function y = operatorPhit(x, Gt, At, No, W)
% This function implements the operator Phi*T = A^T * G^T
% Complex -> Real
x1 = Gt * x;
if exist('No', 'var') && exist('W', 'var')
    tmp = zeros(No(1) * No(2), 1);
    tmp(W) = x1;
    x1 = tmp;
end
y = real(At(x1));
