function y = HS_forward_operator_G(x, G, W, A, flag_visibility_gridding, Sigma)
    % Apply the forward wideband measurement operator (with or w/o
    % data dimensionality reduction, adjoint of
    % :mat:func:`lib.operators.HS_adjoint_operator_G`).
    %
    % Parameters
    % ----------
    % x : double[:, :, :]
    %     Wideband image.
    % G : cell of cell of sparse complex[:, :]
    %     Degridding matrix  or holographic matrix if visibility gridding is active (per channel per block).
    % W : cell of cell of int[:]
    %     Selection vector to map data blocks to the full Fourier plane.
    % A : anonymous function
    %     Weighted FFT involved in the NUFFT.
    % flag_visibility_gridding : bool
    %     Flag indicating whether data dimensionality reduction  via visibility gridding is considered.
    % Sigma : cell of cell of double[:]
    %     Weighting matrix involved in data dimensionality reduction via visibility gridding.
    %
    % Returns
    % -------
    % y : cell of cell of complex[:]
    %     Output visibilities.
    %

    if ~exist('flag_visibility_gridding', 'var')
        flag_visibility_gridding = 0;
        Sigma = [];
    elseif ~exist('Sigma', 'var')
        Sigma = [];
    end
    [~, ~, c] = size(x);
    y = cell(c, 1);

    for i = 1:c
        Fx = A(x(:, :, i));
        y{i} = cell(length(G{i}), 1);
        for j = 1:length(G{i})
            if flag_visibility_gridding
                if istril(G{i}{j})
                    FxSlice = Fx(W{i}{j});
                    y{i}{j} = Sigma{i}{j} .* (G{i}{j} * FxSlice + (FxSlice' * G{i}{j})');
                    FxSlice = [];
                else
                    y{i}{j} = Sigma{i}{j} .* (G{i}{j} * Fx(W{i}{j}));
                end
            else
                y{i}{j} = (G{i}{j} * Fx(W{i}{j}));
            end
        end
    end

end
