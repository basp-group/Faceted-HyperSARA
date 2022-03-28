function x_wave = HS_adjoint_sparsity(x, Psit, b)
    % Apply the SARA dictionary to a wideband image (adjoint of
    % :mat:func:`lib.operators.HS_forward_sparsity`).
    %
    % Parameters
    % ----------
    % x : double[:, :, :]
    %     Wideband image [N(1), N(2), L]
    % Psit : anonymous function
    %     Function handle representing the serial implementation of the
    %     SARA dictionary.
    % b : int
    %     Number of coefficients obtained after application of the SARA
    %     dictionary to a monochromatic image.
    %
    % Returns
    % -------
    % x_wave : double[:, :]
    %     Wavelet dictionaries obtained for each channel [b, L].
    %

    [~, ~, c] = size(x);
    x_wave = zeros(b, c);

    for i = 1:c
        x_wave(:, i) = Psit(x(:, :, i));
    end

end
