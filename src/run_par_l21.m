function [v1, g1] = run_par_l21(v1, Psit, Psi, xhat, weights1, beta1, sigma1)

r1 = zeros(size(v1));
g1 = zeros(size(xhat));

for l = 1:size(xhat, 3)
    r1(:, l) = v1(:, l) + Psit(xhat(:,:,l));
end
l2 = sqrt(sum(abs(r1).^2,2));
l2_soft = max(l2 - beta1*weights1, 0)./ (l2+eps);
v1 = r1 - (l2_soft .* r1);

for l = 1:size(xhat, 3)
    g1(:,:,l) = sigma1*Psi(v1(:,l));
end

end