function weights = update_weights_nuclear_serial(x, reweight_alpha)

[M,  N, c] = size(x);
[~,D,~] = svd(reshape(x, [M*N, c]),'econ');
d = abs(diag(D));
weights = reweight_alpha ./ (reweight_alpha + d);

end