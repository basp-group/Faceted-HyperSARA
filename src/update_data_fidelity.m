function [v2, Ftx, proj, norm_res, global_norm_res, norm_epsilon] = ...
    update_data_fidelity(v2, y, xhat, proj, A, At, G, W, pU, epsilon, ...
                         elipse_proj_max_iter, elipse_proj_min_iter, ...
                         elipse_proj_eps, sigma22)
%update_data_fidelity: update step of the data fidelity term (with 
% preconditioning -> projection onto an ellipsoid).
%-------------------------------------------------------------------------%
%%
% Reference:
%
% [1] Dabbech, A. and Onose, A. and Abdulaziz, A. and Perley, R. A. and 
% Smirnov, O. M. and Wiaux, Y.Cygnus A Super-Resolved via Convex 
% Optimization from VLA Data, Mon. Not. R. Astron. Soc., 476(3), 
% pp. 2853-2866
%-------------------------------------------------------------------------%
%%
% Input:
%
% > v2                              data fidelity dual variable {L}{nblocks}[M, 1]
% > y                               blocks of visibilities {L}{nblocks}[M, 1]
% > xhat                            auxiliary image variable [N(1), N(2), L]
% > proj                            ... {L}{nblocks}[M, 1]
% > A                               measurement operator @[1]
% > At                              adjoint measurement operator @[1]
% > G                               interpolation matrix
% > W                               masking operator {L}{nblocks}
% > pU                              preconditioning matrices {L}{nblocks}
% > epsilon                         l2-ball radius {L}{nblocks}
% > elipse_proj_max_iter            maximum number of iterations (projection)
%                                   [1]
% > elipse_proj_min_iter            minimum number of iterations (projection)
%                                   [1]
% > elipse_proj_eps                 stopping criterion for the projection [1]
% > sigma22                         update step for the dual variable (tau*sigma2)
%
% Output:
%
% < v2                              data fidelity dual variable {L}{nblocks}[M, 1]
% < Ftx                             auxiliary variable for the update of
%                                   the primal variable [N(1), N(2)]
% < proj                            ... {L}{nblocks}[M, 1]
% < norm_res                        norm of the residual {L}{nblocks}[1]    
% < global_norm_res                 square global norm of the residual [1]
% < norm_epsilon                    square global norm of epsilon [1]
%
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%                         
             
Ftx = zeros(size(xhat));
nChannels = size(xhat, 3);
norm_res = cell(nChannels, 1);
norm_epsilon = 0;
global_norm_res = 0;
for i = 1 : nChannels
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
        global_norm_res = global_norm_res + norm_res{i}{j}^2;
        norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
    end
    Ftx(:,:,i) = real(At(g2));
end
Ftx = sigma22*Ftx;

end
