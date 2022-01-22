function y = HS_forward_operator_precond(x, Gw, A, aW)
% Apply the forward preconditioned wideband measurement operator 
% (w/o data blocking, adjoint of 
% :mat:func:`lib.operators.HS_adjoint_operator_precond`).
% 
% Parameters
% ----------
% x : double[:, :, :]
%     Wideband image.
% Gw : cell of sparse complex[:, :]
%     Degridding matrix (per channel).
% A : anonymous function
%     Weighted FFT involved in the NUFFT.
% aW : cell of double[:]
%     Diagonal preconditioning matrices (encoded as vector, each cell entry
%     associated with a different channel).
% 
% Returns
% -------
% y : cell of complex[:]
%     Output visibilities.
%

    [~, ~, c] = size(x);
    y = cell(c, 1);

    for ind = 1:c
        y{ind} =  sqrt(cell2mat(aW{ind})) .* (Gw{ind} * A(x(:, :, ind)));
    end

end
