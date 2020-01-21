function x = HS_adjoint_operator(y, Gw, At,N, M)

% Parameters
c = length(y);
x = zeros(N, M, c);

%
for ind = 1:c
    x(:,:,ind) = At(Gw{ind}' * y{ind});
end

end
