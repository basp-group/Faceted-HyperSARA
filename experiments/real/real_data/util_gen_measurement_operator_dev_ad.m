function [A, At, G, W, aW] = util_gen_measurement_operator_dev_ad(u, v, w, nWw, ...
    param_precond, param_blocking, nchans, Nx, Ny, param_nufft, param_wproj, param_preproc)
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Parameters
% ----------
% u : array (vector)
%     `u` coverage.
% v : array (vector)
%     `v` coverage.
% param_precond : struct
%     Structure to configure the preconditioning matrices.
% param_blocking : struct
%     Structure to configure data blocking.
% fc : array (vector)
%     Frequencies at which the operator needs to be defined.
% Nx : int
%     Image dimension (x-axis).
% Ny : int
%     Image dimension (y-axis).
% Kx : int
%     Dimension interpolation kernel (x-axis).
% Ky : int
%     Dimension interpolation kernel (y-axis).
% ox : int
%     Fourier oversampling factor(x-axis).
% oy : int
%     Fourier oversampling factor(y-axis).
% kernel : string
%     Type of interpolation kernel selected ('kaiser' or 'minmax:tuned').
%
% Returns
% -------
% A : lambda function
%     Lambda function to compute the rescaled 2D Fourier transform involved
%     in the emasurement operator.
% At : lambda function
%     Lambda function to compute the adjoint of ``A``.
% G : cell
%     Cell containing the trimmed-down interpolation kernels for each
%     channel, and each data block within a channel.
% W : cell
%     Cell containing the selection vector for each channel, and
%     data block within a channel.
% aWw : cell
%     Cell containing the preconditioning vectors for each channel, and
%     data block within a channel.
%%
param_nufft.N = [Ny Nx];
param_nufft.Nn = [param_nufft.Ky param_nufft.Kx];
param_nufft.No = [param_nufft.oy * Ny param_nufft.ox * Nx];
param_nufft.Ns = [Ny / 2 Nx / 2];

G = cell(nchans, 1);
W = cell(nchans, 1);
aW = cell(nchans, 1);

param.gen_only_fft_op = 1;
[A, At, ~, ~] = op_p_nufft_wproj_dde([0 0], 1, param_nufft, param_wproj, param);

for i = 1:nchans
    % set the blocks structure
    if ~isempty(param_blocking)
        % compute uniform weights (sampling density) for the preconditioning
        aWw = util_gen_preconditioning_matrix(u{i}, v{i}, param_precond);
        % set the weighting matrix, these are the natural weights for real data
        %     nWw = ones(length(u{i}), 1);
        [u{i}, v{i}, ~, ~, aW{i}, nW] = util_gen_block_structure(u{i}, v{i}, aWw{i}, nWw{i}, param_blocking);
    else
    aW{i} = cell(numel(u{i}), 1);
        for j = 1:numel(u{i})
            aW{i}{j} = util_gen_preconditioning_matrix(u{i}{j}, v{i}{j}, param_precond);
        end
        nW =  nWw{i};
    end
    % measurement operator initialization
    % [A, At, G{i}, W{i}] = op_p_nufft_irt([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW, kernel);

    param_nufft.nW = nW; % [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW %(p, N, Nn, No, Ns, ww, param)
    if ~(param_preproc.done)
        [~, ~, G{i}, W{i}] = op_p_nufft_wproj_dde([v{i} u{i}], w{i}, param_nufft, param_wproj);
    else
    load(param_preproc.G_filename(param_preproc.subcube, param_preproc.ch(i)), 'Gw');
        G{i} = cell(numel(u{i}), 1);
    W{i} = cell(numel(u{i}), 1);
        for j = 1:numel(u{i})

        W{i}{j} = (sum(abs(Gw{j}), 1) > 0).'; % , 1).';
            G{i}{j} = Gw{j}(:, W{i}{j});
            Gw{j} = [];
        end
    clear Gw nW;
    end
end

end
