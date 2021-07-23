function [v2, Ftx, Fx_old, proj, norm_res, global_norm_res, norm_epsilon] = ...
    update_dual_fidelity(v2, y, x, Fx_old, proj, A, At, G, W, pU, epsilon, ...
                         elipse_proj_max_iter, elipse_proj_min_iter, ...
                         elipse_proj_eps, sigma22)
% Update the data fidelity term in the preconditioned primal-dual algorithm.
%
% Update the data fidelity terms owned by each worked involved in the group
% of data nodes (with preconditioning -> projection onto an ellipsoid).
%
% Args:
%     v2 (cell): data fidelity dual variable {L}{nblocks}[M, 1].
%     y (cell): blocks of visibilities {L}{nblocks}[M, 1].
%     x (array): primal variable [N(1), N(2), L].
%     Fx_old (array): ... [N(1), N(2), L].
%     proj (cell): value of the projection at the previous global 
%                    iteration, taken as a starting {L}{nblocks}[M, 1].
%     A (lambda): measurement operator @[1].
%     At (lambda): adjoint measurement operator @[1].
%     G (cell): blocked interpolation matrix {L}{nblocks}.
%     W (cell): blocked masking operator {L}{nblocks}.
%     pU (cell): preconditioning matrices {L}{nblocks}.
%     epsilon (cell): l2-ball constraints {L}{nblocks}[1].
%     elipse_proj_max_iter (int): maximum number of iterations 
%                                 (projection onto ellipsoid).
%     elipse_proj_min_iter (int): minimum number of iterations 
%                                 (projection onto ellipsoid).
%     elipse_proj_eps (int): stopping criterion for the projection.
%     sigma22 (array): step-size for the update of the dual variable 
%                      (tau*sigma2)
%
% Returns:
%     v2 (cell): data fidelity dual variable {L}{nblocks}[M, 1].
%     Ftx (array): auxiliary variable for the update of the primal 
%       variable [N(1), N(2)].
%     Fx_old: [N(1), N(2), L]
%     proj (cell): result of the projection step  {L}{nblocks}[M, 1].
%     norm_res (cell): norm of the residual {L}{nblocks}[1].    
%     global_norm_res (double): square global norm of the residual.
%     norm_epsilon (double): square global norm of epsilon.

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
% Code: P.-A. Thouvenin.
% Last revised: [29/04/2021]
% TODO: update name of variable
%-------------------------------------------------------------------------%
%%                         
             
Ftx = zeros(size(x));
nChannels = size(x, 3);
norm_res = cell(nChannels, 1);
norm_epsilon = 0;
global_norm_res = 0;
for i = 1 : nChannels
    Fx = A(x(:,:,i));
    g2 = zeros(size(Fx));
    norm_res{i} = cell(length(G{i}), 1);
    for j = 1 : length(G{i})
        r2 = G{i}{j} * (2*Fx(W{i}{j}) - Fx_old(W{i}{j}, i));
        proj{i}{j} = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
            r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, ...
            elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
        v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2 - pU{i}{j} .* proj{i}{j};        
        u2 = G{i}{j}' * v2{i}{j};
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
        
        norm_res{i}{j} = norm(G{i}{j}*Fx(W{i}{j}) - y{i}{j}, 2);
        global_norm_res = global_norm_res + norm_res{i}{j}^2;
        norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
    end
    Fx_old(:,i) = Fx; 
    Ftx(:,:,i) = sigma22(i)*real(At(g2));
end

end
