function [v2, Ftx, proj, norm_res, norm_epsilon] = update_data_fidelity_dr(v2, y, xhat, proj, A, At, T, W, pU, epsilon, ...
                     elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps, sigma22)

% 1. the structure of this function needs to be completely changed in the case
% of multiple, distributed data blocks for each ingle frequency (much more completated)
% 2. possibly hange the structure of the data given here...
                 
Ftx = zeros(size(xhat));
N = numel(xhat(:,:,1));
ci = size(xhat, 3);
norm_res = 0;
norm_epsilon = 0;
for i = 1 : ci
    Fx = A{i}(xhat(:,:,i));
    g2 = zeros(N, 1);
    for j = 1 : length(T{i})
        r2 = T{i}{j} .* Fx(W{i}{j});
        proj{i}{j} = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
            r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
        v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2 - pU{i}{j} .* proj{i}{j};        
        u2 = T{i}{j} .* v2{i}{j};
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
        
        norm_res = norm_res + power(norm(r2 - y{i}{j}), 2);
        norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
    end
    Ftx(:,:,i) = real(At{i}(g2));
end
Ftx = sigma22*Ftx;

end