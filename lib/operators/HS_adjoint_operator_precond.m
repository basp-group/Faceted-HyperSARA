function x = HS_adjoint_operator_precond(y, Gw, At, aW, N, M)
% Apply adjoint of the preconditioned wideband measurement operator 
% (w/o data blocking, adjoint of 
% :mat:func:`lib.operators.HS_forward_operator_precond`).
% 
% Parameters
% ----------
% y : cell of complex[:]
%     Input visibilities.
% Gw : cell of sparse complex[:, :]
%     Degridding matrix (per channel).
% At : anonymous function
%     Weighted iFFT involved in the adjoint NUFFT.
% aW : cell of double[:]
%     Diagonal preconditioning matrices (encoded as vector, each cell entry
%     associated with a different channel).
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
        x(:, :, ind) = At(Gw{ind}' * (sqrt(cell2mat(aW{ind})) .* y{ind}));
    end

end
