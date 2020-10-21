function f = nuclear_norm(x)
    [~,D,~] = svd(reshape(x, [numel(x)/size(x, 3), size(x, 3)]),'econ');
    f = sum(abs(diag(D)));
end