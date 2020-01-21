function x = HS_adjoint_operator_precond(y, Gw, At, aW, N, M)

% Parameters
c = length(y);
x = zeros(N, M, c);

%
for ind = 1:c
    x(:,:,ind) = At(Gw{ind}' * (sqrt(cell2mat(aW{ind})) .* y{ind}));
end

end
