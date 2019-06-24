function res = compute_residual_images_dr(x, y, T, A, At, W)

ci = size(x, 3);
res = zeros(size(x));

for i = 1 : ci
    Fx = A{i}(x(:,:,i));
    g2 = zeros(numel(Fx),1);
    for j = 1 : length(T{i})
        res_f = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
        u2 = T{i}{j} .* res_f;
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
    end
    res(:,:,i) = real(At{i}(g2));
end

end
