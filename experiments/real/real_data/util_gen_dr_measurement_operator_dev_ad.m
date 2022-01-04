function [A, At, H, W, aW, Sigma, data, noise] = util_gen_dr_measurement_operator_dev_ad(y, u, v, w, nWw, ...
    param_precond, param_blocking, nchans, Nx, Ny, param_nufft, param_wproj, preproc_dr_residuals, ddes)
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
noise = [];
param_nufft.N = [Ny Nx];
param_nufft.Nn = [param_nufft.Ky param_nufft.Kx];
param_nufft.No = [param_nufft.oy * Ny param_nufft.ox * Nx];
param_nufft.Ns = [Ny / 2 Nx / 2];
res_data =  cell(nchans, 1);
H = cell(nchans, 1);
W = cell(nchans, 1);
aW = cell(nchans, 1);
data =  cell(nchans, 1);
Sigma =  cell(nchans, 1);
param.gen_only_fft_op = 1;
[A, At, ~, ~] = op_p_nufft_wproj_dde([0 0], 1, param_nufft, param_wproj, param);

for i = 1:nchans

    % set the blocks structure
    if ~isempty(param_blocking)
        % compute uniform weights (sampling density) for the preconditioning
        aWw = util_gen_preconditioning_matrix(u{i}, v{i}, param_precond);
        % set the weighting matrix, these are the natural weights for real data
        [u{i}, v{i}, ~, ~, aW{i}, nW] = util_gen_block_structure(u{i}, v{i}, aWw, nWw{i}, param_blocking);
        aWw = [];
    else
        % no -precond for DR ..
        for j = 1:numel(u{i})
            aW{i}{j, 1} = 1;
        end
        nW =  nWw{i};
    end
    % measurement operator initialization
    H{i}{1} = sparse(prod(param_nufft.No), prod(param_nufft.No));
    W_ch = sparse(prod(param_nufft.No), 1);
    for j = 1:numel(u{i})
        fprintf('\nData Block %d \n', j);
        param_nufft.nW = nW(j);
        if isempty(ddes)
            [~, ~, G, Wcurr] = op_p_nufft_wproj_dde([v{i}(j) u{i}(j)], w{i}(j), param_nufft, param_wproj, []);
        else
            [~, ~, G, Wcurr] = op_p_nufft_wproj_dde([v{i}(j) u{i}(j)], w{i}(j), param_nufft, param_wproj, ddes{i}(j));
        end
        Gcurr = sparse(prod(param_nufft.No), size(G{1}, 1));
        Gcurr(Wcurr{1}, :) = G{1}';
        clear G;
        if  ~isempty(preproc_dr_residuals)
            try % check is residual  data are available to get updatesd l2bounds
                if j == 1; res_data{i}{1} = Gcurr * preproc_dr_residuals{i}{j};
                else; res_data{i}{1} =  res_data{i}{1}  + Gcurr * preproc_dr_residuals{i}{j};
                end
            end
        end

        % gridded data & active Fourier modes & H

        if j == 1; data{i}{1} = Gcurr * y{i}{j};
        else; data{i}{1} =  data{i}{1}  + Gcurr * y{i}{j};
        end
        H{i}{1} = H{i}{1} + tril(Gcurr * Gcurr'); % keep lower half only of H
        clear Gcurr;
        W_ch = sparse((W_ch + Wcurr{1}) > 0);
    end
    W{i}{1} = W_ch; clear  W_ch;

    % apply Sigma to data
    diagH = diag(H{i}{1});
    data{i}{1}(diagH > 0) = data{i}{1}(diagH > 0) ./ sqrt(diagH(diagH > 0));

    % get noise stats
    try
        res_data{i}{1}(diagH > 0) =  res_data{i}{1}(diagH > 0) ./ sqrt(diagH(diagH > 0));
        noise.l2bounds{i}{1} = norm(nonzeros(res_data{i}{1}));
        noise.sigma{i} = std(nonzeros([real(res_data{i}{1}); imag(res_data{i}{1})]));
        dummy =    std([real(res_data{i}{1}); imag(res_data{i}{1})]);
        fprintf('\nResidual data vector is provided, nnz %d', nnz(res_data{i}{1}));
        fprintf('\nnoise estimates are computed. sigma: %f (unmasked res %f ) , l2 norm  %f', ...
            noise.sigma{i}, dummy, noise.l2bounds{i}{1});
        clear dummy;
    end
    res_data{i}{1} = [];

    W{i}{1} = (diagH > 0);
    H{i}{1} = H{i}{1}(W{i}{1}, W{i}{1});
    data{i}{1} = data{i}{1}(W{i}{1});
    Sigma{i}{1} = 1 ./ sqrt(diagH(W{i}{1}));

    if istril(H{i}{1})
        fprintf('\nH is lower triangular (for mem reasons)..updating diag elements of H\n');
        for idiag = 1:size(H{i}{1}, 1)
            H{i}{1}(idiag, idiag) = 0.5 * H{i}{1}(idiag, idiag);
        end
    end
    aW{i}{1} = 1; % no precond ...

end

end
