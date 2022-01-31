function y = HS_forward_operator(x, Gw, A)
    % Apply the forward wideband measurement operator (w/o data blocking,
    % adjoint of :mat:func:`lib.operators.HS_adjoint_operator`).
    %
    % Parameters
    % ----------
    % x : double[:, :, :]
    %     Wideband image.
    % Gw : cell of sparse complex[:, :]
    %     Degridding matrix (per channel).
    % A : anonymous function
    %     Weighted FFT involved in the NUFFT.
    %
    % Returns
    % -------
    % y : cell of complex[:]
    %     Visibilities.
    %

    [~, ~, c] = size(x);
    y = cell(c, 1);

    for ind = 1:c
        y{ind} =  Gw{ind} * A(x(:, :, ind));
    end

end
