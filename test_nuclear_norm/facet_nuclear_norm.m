function f = facet_nuclear_norm(x)

[M, N, L] = size(x);
[~,S,~] = svd(reshape(x, [N*M, L]),'econ');
f = sum(abs(diag(S)));

end