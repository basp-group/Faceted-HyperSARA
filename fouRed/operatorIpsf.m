function y = operatorIpsf(x, A, At, H, No, W)
Ny = No(1);
Nx = No(2);

x1 = A(real(x));
if exist('W', 'var') && exist('No', 'var')
    x1 = x1(W);
end
x2 = H * x1;
if exist('W', 'var') && exist('No', 'var')
    x3 = zeros(Ny * Nx, 1);
    x3(W) = x2;
    x2 = x3;
end
y = real(At(x2));
end