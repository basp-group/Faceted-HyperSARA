function gamma0 = compute_low_rank_regularizer(x, sig_bar, alph_bar, regtype)
% TODO: to be updated

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [...]
%-------------------------------------------------------------------------%

[M,  N, c] = size(x);
[~, D, ~] = svd(reshape(x, [M*N, c]),'econ');

switch regtype
    case "log"
        gamma0 = alph_bar / (sig_bar * sum(log(abs(diag(D))/sig_bar + 1)));
    case "inv"
        gamma0 = alph_bar / sum(abs(diag(D)));
    otherwise
        gamma0 = alph_bar / sum(abs(diag(D)));
end

end
