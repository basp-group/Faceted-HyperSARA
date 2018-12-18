function x = HS_adjoint_operator(y, Gw, At,N, M)

% Parameters
c = size(y, 2);
x = zeros(N, M, c);

%
for ind = 1:c
    x(:, :, ind) = At(Gw{ind}' * y(:,ind));
end

end
