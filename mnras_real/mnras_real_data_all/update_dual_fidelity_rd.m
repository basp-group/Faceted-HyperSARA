function [v2, Ftx, proj, norm_res, global_norm_res_c, norm_epsilon_c, global_norm_res_a, norm_epsilon_a] = update_dual_fidelity_rd(v2, y, xhat, proj, A, At, G, W, pU, epsilon, ...
                     elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps, sigma22)
                 
Ftx = zeros(size(xhat));
ci = size(xhat, 3);
norm_res = cell(ci, 1);
norm_epsilon_c = 0;
global_norm_res_c = 0;
norm_epsilon_a = 0;
global_norm_res_a = 0;
for i = 1 : ci
    Fx = A(xhat(:,:,i));
    g2 = zeros(size(Fx));
    norm_res{i} = cell(length(G{i}), 1);
    for j = 1 : length(G{i})
        r2 = G{i}{j} * Fx(W{i}{j});
        proj{i}{j} = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
            r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
        v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2 - pU{i}{j} .* proj{i}{j};        
        u2 = G{i}{j}' * v2{i}{j};
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
        
        norm_res{i}{j} = norm(r2 - y{i}{j}, 2);

        if j == 1  %%%% C block
            global_norm_res_c = global_norm_res_c + norm_res{i}{j}^2;
            norm_epsilon_c = norm_epsilon_c + power(epsilon{i}{j}, 2);
        else %%%% A block
            global_norm_res_a = global_norm_res_a + norm_res{i}{j}^2;
            norm_epsilon_a = norm_epsilon_a + power(epsilon{i}{j}, 2);
        end
    end
    Ftx(:,:,i) = real(At(g2));
end
Ftx = sigma22*Ftx;

end
