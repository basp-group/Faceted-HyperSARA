function [A, At, G, W, aW] = util_gen_measurement_operator(u, v, ...
    param_precond, param_blocking, fc, fmax, Nx, Ny, Kx, Ky, ox, oy, kernel)
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

nchans = numel(fc);
G = cell(nchans, 1);
W = cell(nchans, 1);
aW = cell(nchans, 1);

for i = 1:nchans

    % need to normalize by the maximum over all the available frequencies
    uw = (fc(i)/fmax) * u;
    vw = (fc(i)/fmax) * v;
    
    % compute uniform weights (sampling density) for the preconditioning
    aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw), 1);
    
    % set the blocks structure
    [u1, v1, ~, ~, aW{i}, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_blocking);
    
    % measurement operator initialization
    [A, At, G{i}, W{i}] = op_p_nufft_irt([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW, kernel);
end

end