function [A, At, G, W, aW] = util_gen_measurement_operator(u, v, w, nW, ...
     nchans, Nx, Ny, param_nufft, param_wproj,param_precond ,ddes)
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Parameters
% ----------
% u : array (vector)
%     `u` coordinate.
% v : array (vector)
%     `v` coordinate.
% w : array (vector)
%     `w` coordinate.
% param_precond : struct
%     Structure to configure the preconditioning matrices.
% param_blocking : struct
%     Structure to configure data blocking.
% param_nufft : struct
%     Structure to configure NUFFT.
% param_wproj : struct
%     Structure to configure w-projection.
% nchans : Number of channels
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
% ddes: array
%     Array of antenna gains to be incorporated in the de-gridding matrix.
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
% aW : cell
%     Cell containing the preconditioning vectors for each channel, and
%     data block within a channel.
%%
if ~exist('ddes','var')
    ddes = [];
end
param_nufft.N = [Ny Nx];
param_nufft.Nn = [param_nufft.Ky param_nufft.Kx];
param_nufft.No = [param_nufft.oy * Ny param_nufft.ox * Nx];
param_nufft.Ns = [Ny / 2 Nx / 2];

G = cell(nchans, 1);
W = cell(nchans, 1);
aW = cell(nchans, 1);

% get fft operators
[A, At, ~, ~] = op_p_nufft_wproj_dde(param_nufft);

for i = 1:nchans
    % set the blocks structure
    aW{i} = cell(numel(u{i}), 1);
    for j = 1:numel(u{i})
        % compute uniform weights (sampling density) for the preconditioning
        aW{i}{j} = util_gen_preconditioning_matrix(u{i}{j}, v{i}{j}, param_precond);
    end
    % measurement operator initialization
    if isempty(ddes)
        [~, ~, G{i}, W{i}] = op_p_nufft_wproj_dde(param_nufft,[v{i} u{i}], w{i}, nW{i}, param_wproj);
    else
        [~, ~, G{i}, W{i}] = op_p_nufft_wproj_dde(param_nufft,[v{i} u{i}], w{i}, nW{i}, param_wproj,ddes{i});
    end
end

end
