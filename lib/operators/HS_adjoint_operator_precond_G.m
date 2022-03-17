function x = HS_adjoint_operator_precond_G(y, G, W, At, aW, N, M, flag_dr, Sigma)
    % Apply adjoint of the preconditioned wideband measurement operator
    % (with or w/o data dimensionality reduction, adjoint of
    % :mat:func:`lib.operators.HS_forward_operator_precond_G`).
    %
    % Parameters
    % ----------
    % y : cell of cell of complex[:]
    %     Input visibilities.
    % G : cell of cell of sparse complex[:, :]
    %     Degridding matrix (per channel per block).
    % W : cell of cell of int[:]
    %     Selection vector to extract data blocks from the full Fourier plane.
    % At : anonymous function
    %     Weighted iFFT involved in the adjoint NUFFT.
    % aW : cell of cell of double
    %     Diagonal preconditioning matrices (encoded as vector, for each
    %     channel and data block within a channel).
    % N : int
    %     Spatial dimension of the wideband image (y axis).
    % M : int
    %     Spatial dimension of the wideband image (x axis).
    % flag_dr : bool
    %     Flag indicating whether data dimensionality reduction is considered.
    % Sigma : cell of cell of double[:]
    %     Weighting matrix involved in data dimensionality reduction.
    %
    % Returns
    % -------
    % x : double[:, :, :]
    %     Output wideband image.
    %

    if ~exist('flag_dr', 'var')
        flag_dr = 0;
        Sigma = [];
    elseif ~exist('Sigma', 'var')
        Sigma = [];
    end

    c = length(y);
    x = zeros(N, M, c);
    % No = size(G{1}{1}, 2);
    No = size(W{1}{1}, 1);

    for i = 1:c
        g2 = zeros(No, 1);
        for j = 1:length(G{i})
            if flag_dr
                if istril(G{i}{j})
                    weighted_y = (sqrt(aW{i}{j}) .* Sigma{i}{j}) .* y{i}{j};
                    g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * weighted_y  +  G{i}{j} * weighted_y;

                else
                    g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * ((sqrt(aW{i}{j}) .* Sigma{i}{j}) .* y{i}{j}); %
                end
            else
                g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (sqrt(aW{i}{j}) .* y{i}{j});
            end
        end
        x(:, :, i) = At(g2);
    end

end
