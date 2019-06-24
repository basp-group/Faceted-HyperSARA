function weights = update_weights_l21_serial(x, Psit, weights, reweight_alpha)

w = zeros(size(weights, 1), size(x, 3));
for l = 1 : size(x, 3)
    w(:, l) = Psit(x(:, :, l));
end
d = sqrt(sum(abs((w)).^2,2));
weights = reweight_alpha ./ (reweight_alpha + d);

end
