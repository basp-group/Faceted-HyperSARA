function x = HS_forward_sparsity(x_wave, Psi, n, m)

    % Parameters
    [~, c] = size(x_wave);
    x = zeros(n, m, c);

    %
    for i = 1:c
        %     dwtmode('per');
        x(:, :, i) = Psi(x_wave(:, i));
    end

end
