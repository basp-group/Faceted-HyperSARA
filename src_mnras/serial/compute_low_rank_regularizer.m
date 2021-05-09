function gamma0 = compute_low_rank_regularizer(x, sig_bar, alph_bar)
% TODO: to be updated

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [...]
%-------------------------------------------------------------------------%

[M,  N, c] = size(x);
[~, D, ~] = svd(reshape(x, [M*N, c]),'econ');
gamma0 = alph_bar / (sig_bar * sum(log(abs(diag(D))/sig_bar + 1)));

end
