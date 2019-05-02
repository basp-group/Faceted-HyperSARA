function res = compute_residual_images(x, y, G, A, At, W)

ci = size(x, 3);
res = zeros(size(x));

for i = 1 : ci
    Fx = A(x(:,:,i));
    g2 = zeros(numel(Fx),1);
    for j = 1 : length(G{i})
        res_f = y{i}{j} - G{i}{j} * Fx(W{i}{j});
        u2 = G{i}{j}' * res_f;
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
    end
    res(:,:,i) = real(At(g2));
end

end