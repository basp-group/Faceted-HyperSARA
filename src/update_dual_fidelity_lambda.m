function [v2, Ftx, Fx_old, proj, norm_res, global_norm_res, norm_epsilon] = ...
    update_dual_fidelity_lambda(v2, y, x, Fx_old, proj, A, At, G, W, pU, epsilon, ...
                         elipse_proj_max_iter, elipse_proj_min_iter, ...
                         elipse_proj_eps, sigma22, apply_G, apply_Gdag)
% Update the data fidelity term in the preconditioned primal-dual algorithm.
%
% Update the data fidelity terms owned by each worked involved in the group
% of data nodes (with preconditioning -> projection onto an ellipsoid).
%
% Parameters
% ----------
% v2 : cell
%     Data fidelity dual variable {L}{nblocks}[M, 1].
% y : cell
%     Blocks of visibilities {L}{nblocks}[M, 1].
% x : array (3d)
%     Primal variable [N(1), N(2), L].
% Fx_old : array (3d)
%     ... [N(1), N(2), L].
% proj : cell
%     Value of the projection at the previous global iteration, taken as a
%     starting point {L}{nblocks}[M, 1].
% A : lambda
%     Measurement operator @[1].
% At : lambda
%     Adjoint measurement operator @[1].
% G : cell
%     Blocked interpolation matrix {L}{nblocks}.
% W : cell
%     Blocked masking operator {L}{nblocks}.
% pU : cell
%     preconditioning matrices {L}{nblocks}.
% epsilon : cell
%     :math:`\ell_2`-ball constraints {L}{nblocks}[1].
% elipse_proj_max_iter : int
%     Maximum number of iterations (projection onto the ellipsoid) 
% elipse_proj_min_iter : int
%     Minimum number of iterations (projection onto the ellipsoid) 
% elipse_proj_eps : double
%     Stopping criterion for the projection.
% sigma22 : array (1d)
%     Step-size for the update of the dual variable (tau*sigma2).
%
% Returns
% -------
% v2 : cell
%     Data fidelity dual variable {L}{nblocks}[M, 1].
% Ftx (array):
%     Auxiliary variable for the update of the primal variable [N(1), N(2)].
% Fx_old
%     [N(1), N(2), L]
% proj (cell)
%     Result of the projection step  {L}{nblocks}[M, 1].
% norm_res (cell)
%     Norm of the residual {L}{nblocks}[1].    
% global_norm_res (double)
%     Square global norm of the residual.
% norm_epsilon (double)
%     Square global norm of epsilon.
%    

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
        r2 = apply_G(2*Fx(W{i}{j}) - Fx_old(W{i}{j}, i), G{i}{j});
        proj{i}{j} = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
            r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, ...
            elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
        v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2 - pU{i}{j} .* proj{i}{j};        
        u2 = apply_Gdag(v2{i}{j}, G{i}{j});
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
        
        norm_res{i}{j} = norm(apply_G(Fx(W{i}{j}), G{i}{j}) - y{i}{j}, 2);
        global_norm_res = global_norm_res + norm_res{i}{j}^2;
        norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
    end
    Fx_old(:,i) = Fx; 
    Ftx(:,:,i) = sigma22(i)*real(At(g2));
end

end
