function y = HS_adjoint_sparsity(x, Psit, b)

% Parameters
[n,m,c] = size(x);
y = zeros(n*m*b,c);

%
for i = 1 : c
    %dwtmode('per');
    y(:,i) = Psit(x(:,:,i));
end

end
