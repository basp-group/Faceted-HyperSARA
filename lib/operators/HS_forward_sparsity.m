function x = HS_forward_sparsity(x_wave, Psi, n, m)
    % Apply the adjoint of the SARA dictionary to a wideband image (adjoint of
    % :mat:func:`lib.operators.HS_adjoint_sparsity`).
    %
    % Parameters
    % ----------
    % x_wave : double[:, :]
    %     Wavelet coefficients in the SARA domain [b, L].
    % Psi : anonymous function
    %     Function handle representing the serial implementation of the adjoint
    %     SARA dictionary.
    % n : int
    %     Spatial dimension of the output wideband image (y axis).
    % m : int
    %     Spatial dimension of the output wideband image (x axis).
    %
    % Returns
    % -------
    % x : double[:, :, :]
    %     Wideband image [n, m, L]
    %

    [~, c] = size(x_wave);
    x = zeros(n, m, c);

    for i = 1:c
        % dwtmode('per');
        x(:, :, i) = Psi(x_wave(:, i));
    end

end
