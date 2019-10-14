function y = dataReduce_degrid(x, Gt,Sigma, Mask)
% the embeddding operator R = \sigma * S * Gt
% Complex -> Complex
x1 = Gt * x;
y = Sigma .* x1(Mask);
