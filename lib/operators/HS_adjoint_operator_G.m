function x = HS_adjoint_operator_G(y, G, W, At, N, M, flag_visibility_gridding, Lambda)
    % Apply adjoint of the wideband measurement operator (with or w/o
    % data dimensionality reduction, adjoint of
    % :mat:func:`lib.operators.HS_forward_operator_G`).
    %
    % Parameters
    % ----------
    % y : cell of cell of complex[:]
    %     Input visibilities.
    % G : cell of cell of sparse complex[:, :]
    %     Degridding matrix  or holographic matrix if visibility gridding is active (per channel per block).
    % W : cell of cell of int[:]
    %     Selection vector to extract data blocks from the full Fourier plane.
    % At : anonymous function
    %     Weighted iFFT involved in the adjoint NUFFT.
    % N : int
    %     Spatial dimension of the wideband image (y axis).
    % M : int
    %     Spatial dimension of the wideband image (x axis).
    % flag_visibility_gridding : bool
    %     Flag indicating whether data dimensionality reduction via visibility gridding is considered.
    % Lambda : cell of cell of double[:]
    %     Weighting matrix involved in data dimensionality reduction via visibility gridding.
    %
    % Returns
    % -------
    % x : double[:, :, :]
    %     Output wideband image.
    %

    if ~exist('flag_visibility_gridding', 'var')
        flag_visibility_gridding = 0;
        Lambda = [];
    elseif ~exist('Lambda', 'var')
        Lambda = [];
    end
    c = length(y);
    x = zeros(N, M, c);
    % No = size(G{1}{1}, 2);
    No = size(W{1}{1}, 1);

    for i = 1:c
        g2 = zeros(No, 1);
        for j = 1:length(G{i})
            if flag_visibility_gridding
                if istril(G{i}{j})
                    weighted_data = (Lambda{i}{j} .* y{i}{j});
                    g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * weighted_data  +  G{i}{j} * weighted_data;
                else
                    g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (Lambda{i}{j} .* y{i}{j});
                end
            else
                g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (y{i}{j});
            end
        end
        x(:, :, i) = At(g2);
    end

end
