function x = HS_adjoint_operator(y, Gw, At, N, M)
    % Apply adjoint of the wideband measurement operator (w/o data blocking,
    % adjoint of :mat:func:`lib.operators.HS_forward_operator`).
    %
    % Parameters
    % ----------
    % y : cell of complex[:]
    %     Input visibilities.
    % Gw : cell of sparse complex[:, :]
    %     Degridding matrix (per channel).
    % At : anonymous function
    %     Weighted iFFT involved in the adjoint NUFFT.
    % N : int
    %     Spatial dimension of the wideband image (y axis).
    % M : int
    %     Spatial dimension of the wideband image (x axis).
    %
    % Returns
    % -------
    % x : double[:, :, :]
    %     Output wideband image.
    %
    c = length(y);
    x = zeros(N, M, c);

    for ind = 1:c
        x(:, :, ind) = At(Gw{ind}' * y{ind});
    end

end
