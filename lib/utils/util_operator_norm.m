function [Anorm, squared_operator_norm, rel_var, squared_precond_operator_norm, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, tol, max_iter)

    nchans = numel(G);
    squared_precond_operator_norm = zeros(nchans, 1);
    rel_var_precond = zeros(nchans, 1);

    for l = 1:nchans
        F = afclean(@(x) HS_forward_operator_precond_G(x, G(l), W(l), A, aW(l)));
        Ft = afclean(@(y) HS_adjoint_operator_precond_G(y, G(l), W(l), At, aW(l), Ny, Nx));
        [squared_precond_operator_norm(l), rel_var_precond(l)] = op_norm(F, Ft, [Ny, Nx], tol, max_iter, 0);
    end

    % ! operator is block diagonal: it norm is thus given by the max singular value
    % ! make sure the pdfb convergence criterion is strictly satisfied
    % ! by taking a smaller step-size (take precision of the estimation
    % ! of the operator norm into account)
    Anorm = max(squared_precond_operator_norm .* (1 + rel_var_precond));

    squared_operator_norm = zeros(nchans, 1);
    rel_var = zeros(nchans, 1);
    for l = 1:nchans
        F = afclean(@(x) HS_forward_operator_G(x, G(l), W(l), A));
        Ft = afclean(@(y) HS_adjoint_operator_G(y, G(l), W(l), At, Ny, Nx));
        [squared_operator_norm(l), ~] = op_norm(F, Ft, [Ny, Nx], tol, max_iter, 0);
    end

end
