function [v0, g0] = run_par_nuclear_spmd(v0, xhat, weights0, beta0)

    [M,N,c] = size(xhat);
    xhatm = reshape(xhat,numel(xhat)/c,c);
    [U0,S0,V0] = svd(v0 + xhatm,'econ');
    v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
    g0 = reshape(v0,M,N,c);

end