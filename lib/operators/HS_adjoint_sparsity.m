function x_wave = HS_adjoint_sparsity(x, Psit, b)

    [~, ~, c] = size(x);
    x_wave = zeros(b, c);

    for i = 1:c
        x_wave(:, i) = Psit(x(:, :, i));
    end

end
