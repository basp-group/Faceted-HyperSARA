function [A, At, G, W] = op_p_nufft_calib(p, D, N, Nn, No, Ns, S, ww, param)
    % Create the nonuniform gridding matrix and fft operators to be used for
    % parallel processing
    %
    % Parameters
    % ----------
    % p : cell of double[:, 2]
    %     Nonuniformly distributed frequency location points for each
    %     cell member which will be treated in parallel.
    % D : complex[:, :]
    %     DDE kernels (in the Fourier domain).
    % N : int[2]
    %     Size of the reconstructed image.
    % Nn : int[2]
    %     Size of the kernels (number of neighbors considered on each direction).
    % No : int[2]
    %     Oversampled fft from which to recover the non uniform fft via
    %     kernel convolution.
    % Ns : int[2]
    %     FFT shift.
    % S : int[1]
    %     Size of the DDE support (square support).
    % ww : cell of double[:, 2]
    %     :`math`:`w` component for each cell member which will be treated in
    %     parallel.
    % param : struct
    %     Parameters to build the degridding matrices, composed of the following
    %     fields.
    % param.use_nufft_blocks : bool
    %     Flag to specify NUFFT blocks are used (data blocking).
    % param.gen_only_fft_op : bool
    %     Flag to only generate the weighted FFT and inverse FFT involved in the
    %     the forward and the adjoint NUFFT operators, respectively.
    %
    % Returns
    % -------
    % A : function handle
    %     Function handle for direct operator.
    % At : function handle
    %     Function handle for adjoint operator.
    % G : cell of sparse complex[:, :]
    %     Convolution kernel matrix (small) associated with each patch in the
    %     Fourier plane.
    % W : cell of int[:]
    %     Mask of the values that contribute to the convolution.
    % Gw : sparse complex[:, :]
    %     Global convolution kernel matrix.
    %

    if ~exist('param', 'var')
        param = struct();
    end
    if ~isfield(param, 'use_nufft_blocks')
        param.use_nufft_blocks = 1;
    end
    if ~isfield(param, 'gen_only_fft_op')
        param.gen_only_fft_op = 0;
    end
    if ~exist('ww', 'var')
        ww = cell(length(p), 1);
        for q = 1:length(p)
            ww{q} = ones(length(p{q}(:, 1)), 1);
        end
    end

    R = size(p, 1);
    if param.gen_only_fft_op
        [A, At, ~, ~] = op_nufft([0, 0], N, Nn, No, Ns);
        G = [];
        W = [];
        Gw = [];
    else
        if ~param.use_nufft_blocks
            %% compute the overall gridding matrix and its associated kernels
            % [A, At, Gw, ~] = op_nufft(cell2mat(p), N, Nn, No, Ns);
            [A, At, Gw, ~] = op_nufft_calib(p, D, N, Nn, No, Ns, S);

            %% compute small gridding matrices associated with each parallel block
            G = cell(R, 1);
            W = cell(R, 1);

            % block start position
            fprintf('\nComputing block matrices ...\n');
            b_st = 1;
            for q = 1:R
                tstart = tic;
                % current block length
                % the matrix Gw is structured identical to the structure of p thus we
                % grab it block by block
                b_l = length(p{q});

                % get a block out of the large G and trim it
                Gw(b_st:b_st + b_l - 1, :) = spdiags(ww{q}, 0, b_l, b_l) * Gw(b_st:b_st + b_l - 1, :);
                Gb = Gw(b_st:b_st + b_l - 1, :);

                %% now trim the zero rows and store a mask in W

                % preallocate W for speed
                %             W{q} = false(No(1)*No(2), 1);

                % use the absolute values to speed up the search
                %             Gb_a = abs(Gb);

                % check if eack line is entirely zero
                %             W{q} = Gb_a' * ones(b_l, 1) ~= 0;
                W{q} = any(Gb, 1).';

                % store only what we need from G
                %             G{q} = Gb(:, W{q});
                G{q} = Gb;

                % iterate among the blocks
                b_st = b_st + b_l;
                tend = toc(tstart);
                fprintf('Block matrix %d: %ds \n', q, ceil(tend));
            end
        else

            %% compute small gridding matrices associated with each parallel block

            % Gw = spalloc(length(cell2mat(p)), No(1)*No(2), 16 * length(cell2mat(p)));
            G = cell(R, 1);
            W = cell(R, 1);

            b_st = 1;
            % block start position
            fprintf('\nComputing block matrices ...\n');
            for q = 1:R

                tstart = tic;
                b_l = length(p{q});

                %% compute the small gridding matrix and its associated kernels
                % [~, ~, Gb, ~] = op_nufft([p{q, 1} p{q, 2}], N, Nn, No, Ns);
                [~, ~, Gb, ~] = op_nufft_calib([p{q, 1} p{q, 2}], D, N, Nn, No, Ns, S);

                %% now trim the zero rows and store a mask in W

                % preallocate W for speed
                W{q} = false(No(1) * No(2), 1);

                Gb = spdiags(ww{q}, 0, b_l, b_l) * Gb;

                % use the absolute values to speed up the search
                %             Gb_a = abs(Gb);

                % check if eack line is entirely zero
                %             W{q} = Gb_a' * ones(size(Gb, 1), 1) ~= 0;
                W{q} = any(Gb, 1).';

                % store only what we need from G
                %             G{q} = Gb(:, W{q});
                G{q} = Gb;

                %% fill the whole Gw
                % Gw(b_st:b_st+b_l-1, :) = Gb;

                b_st = b_st + b_l;
                tend = toc(tstart);
                fprintf('Block matrix %d: %ds \n', q, ceil(tend));
            end

            [A, At, ~, ~] = op_nufft([0, 0], N, Nn, No, Ns);
        end
    end

end
