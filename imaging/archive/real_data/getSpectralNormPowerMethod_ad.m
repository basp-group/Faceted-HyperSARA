function val = getSpectralNormPowerMethod_ad(A, At, im_size, epsilon, nmax, vargin)
% Computes the maximum eigen value of the compund
% operator AtA
%
% % epsilon = 10^(-8) ;
% % nmax = 800;
x = randn(im_size);
x_norm = x / norm(x(:));

p = 1 + 10^(-6);
pnew = 1;

n = 1;

cond = abs(pnew - p) / pnew;

% Iterations

while cond >= epsilon && n < nmax

    xnew = At(A(x_norm)); % input normalised

    cond = norm(xnew(:) - x(:))  / norm(x(:));

    x = xnew;
    x_norm = x / norm(x(:));

    n = n + 1;
end

val = norm(xnew(:)); % output not normalised

if nargin > 4
    if vargin
        display(['norm it n=' num2str(n) ', rel_obj=' num2str(cond)]);
    end
end
end
