function y = operatorPhi(x, G, A, W)
% This function implements the operator Phi = G * A
x1 = A(x);
if exist('W', 'var')
    x1 = x1(W);
end
y = G * x1;