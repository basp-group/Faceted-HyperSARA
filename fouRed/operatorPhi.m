function y = operatorPhi(x, G, A, W)
% This function implements the operator Phi = G * A
% Real -> Complex
x1 = A(real(x));
if exist('W', 'var')
    x1 = x1(W);
end
y = G * x1;