function y = HS_forward_operator_precond_G(x, G, W, A, aW, flag_dr, Sigma)
    % Apply the forward preconditioned wideband measurement operator
    % (with or w/o data dimensionality reduction, adjoint of
    % :mat:func:`lib.operators.HS_adjoint_operator_precond_G`).
    %
    % Parameters
    % ----------
    % x : double[:, :, :]
    %     Wideband image.
    % G : cell of cell of sparse complex[:, :]
    %     Degridding matrix (per channel per block).
    % W : cell of cell of int[:]
    %     Selection vector to map data blocks to the full Fourier plane.
    % A : anonymous function
    %     Weighted FFT involved in the NUFFT.
    % aW : cell of cell of double
    %     Diagonal preconditioning matrices (encoded as vector, for each
    %     channel and data block within a channel).
    % flag_dr : bool
    %     Flag indicating whether data dimensionality reduction is considered.
    % Sigma : cell of cell of double[:]
    %     Dimensionality reduction matrix.
    %
    % Returns
    % -------
    % y : cell of cell of complex[:]
    %     Output visibilities.
    %

    if ~exist('flag_dr', 'var')
        flag_dr = 0;
        Sigma = [];
    elseif ~exist('Sigma', 'var')
        Sigma = [];
    end

    [~, ~, c] = size(x);
    y = cell(c, 1);

    for i = 1:c
        Fx = A(x(:, :, i));
        y{i}  = cell(size(G{i}));
        for j = 1:length(G{i})
            if flag_dr
                if  istril(G{i}{j})
                    y{i}{j} = (sqrt(aW{i}{j}) .* (Sigma{i}{j})) .* (G{i}{j} * Fx(W{i}{j}) + G{i}{j}' * Fx(W{i}{j}));
                else
                    y{i}{j} = (sqrt(aW{i}{j}) .* (Sigma{i}{j})) .* (G{i}{j} * Fx(W{i}{j}));
                end

            else
                y{i}{j} = sqrt(aW{i}{j}) .* (G{i}{j} * Fx(W{i}{j}));
            end
        end
    end

end
