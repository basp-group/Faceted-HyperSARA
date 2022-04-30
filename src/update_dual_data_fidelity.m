function [v2, Ftx, Fx_old, proj, norm_res, global_norm_res, norm_epsilon] = ...
    update_dual_data_fidelity(v2, y, x, Fx_old, proj, A, At, G, W, pU, ...
    epsilon, elipse_proj_max_iter, ...
    elipse_proj_min_iter, ...
    elipse_proj_eps, sigma22, ...
    flag_dimensionality_reduction, Lambda)
% Update the data fidelity term in the preconditioned primal-dual algorithm.
%
% Update the data fidelity terms owned by each worked involved in the group
% of data nodes (with preconditioning :math:`\rightarrow` projection onto
% an ellipsoid).
%
% Parameters
% ----------
% v2 : cell
%     Data fidelity dual variable ``{L}{nblocks}[M, 1]``.
% y : cell
%     Blocks of visibilities ``{L}{nblocks}[M, 1]``.
% x : double[:, :, :]
%     Primal variable ``[N(1), N(2), L]``.
% Fx_old : complex[:, :, :]
%     Scaled Fourier transform of :math:`x` from the previous iterations
%     ``[N(1), N(2), L]``.
% proj : cell
%     Value of the projection at the previous global iteration, taken as a
%     starting point ``{L}{nblocks}[M, 1]``.
% A : anonymous function
%     Measurement operator @[1].
% At : anonymous function
%     Adjoint measurement operator @[1].
% G : cell
%     Blocked interpolation matrix {L}{nblocks}.
% W : cell
%     Blocked masking operator {L}{nblocks}.
% pU : cell
%     Preconditioning matrices {L}{nblocks}.
% epsilon : cell
%     :math:`\ell_2`-ball constraints {L}{nblocks}[1].
% elipse_proj_max_iter : int
%     Maximum number of iterations (projection onto the ellipsoid)
% elipse_proj_min_iter : int
%     Minimum number of iterations (projection onto the ellipsoid)
% elipse_proj_eps : double
%     Stopping criterion for the projection.
% sigma22 : double[:]
%     Step-size for the update of the dual variable (tau*sigma2).
% flag_dimensionality_reduction : bool
%     Flag to activate DR functionality.
% Lambda : cell
%     Dimensionality reduction weights {L}{nblocks}.
%
% Returns
% -------
% v2 : cell
%     Data fidelity dual variable ``{L}{nblocks}[M, 1]``.
% Ftx : double[:, :]
%     Auxiliary variable for the update of the primal variable ``[N(1), N(2)]``.
% Fx_old
%     Scaled Fourier transform of the updated image ``[N(1), N(2), L]``.
% proj : cell
%     Result of the projection step  ``{L}{nblocks}[M, 1]``.
% norm_res : cell
%     Norm of the residual ``{L}{nblocks}[1]``.
% global_norm_res : double
%     Square global norm of the residual.
% norm_epsilon : double
%     Square global norm of epsilon.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [29/04/2021]
% TODO: update name of variable
% -------------------------------------------------------------------------%
%%

Ftx = zeros(size(x));
n_channels = size(x, 3);
norm_res = cell(n_channels, 1);
norm_epsilon = 0;
global_norm_res = 0;
sc = @(z, radius) z * (1 - min(radius / norm(z(:)), 1));

if flag_dimensionality_reduction
    for i = 1:n_channels
        Fx = A(x(:, :, i));
        g2 = zeros(size(Fx, 1), size(Fx, 2));
        norm_res{i} = cell(length(G{i}), 1);
        for j = 1:length(G{i})
            % no-precond.
            FxPrev = (2 * Fx(W{i}{j}) - Fx_old(W{i}{j}, i));
            if istril(G{i}{j})
                v2{i}{j} = sc(v2{i}{j} + Lambda{i}{j} .* (G{i}{j} * FxPrev + (FxPrev' * G{i}{j})') - y{i}{j}, epsilon{i}{j});
            else
                v2{i}{j} = sc(v2{i}{j} + Lambda{i}{j} .* (G{i}{j} * FxPrev) - y{i}{j}, epsilon{i}{j});
            end; clear FxPrev;

% % % %             %preconditioning
% % % %             proj{i}{j} = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
% % % %                 r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, ...
% % % %                 elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
% % % %             v2{i}{j} = v2{i}{j} + pU{i}{j} .* (r2 - proj{i}{j});
% % % %             clear r2;

            if istril(G{i}{j})
                u2 = Lambda{i}{j} .* v2{i}{j};
                g2(W{i}{j}) = g2(W{i}{j}) +  (u2' * G{i}{j})' + G{i}{j} * u2;
                clear u2;
            else;  g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (Lambda{i}{j} .* v2{i}{j});
            end

            if istril(G{i}{j})
                FxSlice = Fx(W{i}{j});
                norm_res{i}{j} = sqrt(sum(abs(Lambda{i}{j} .* (G{i}{j} * FxSlice + (FxSlice' * G{i}{j})') - y{i}{j}).^2));
                clear FxSlice;
            else; norm_res{i}{j} = sqrt(sum(abs(Lambda{i}{j} .* (G{i}{j} * Fx(W{i}{j})) - y{i}{j}).^2));
            end

            global_norm_res = global_norm_res + norm_res{i}{j}^2;
            norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
        end
        Fx_old(:, i) = Fx; clear Fx;
        Ftx(:, :, i) = sigma22(i) * real(At(g2)); clear g2;
    end
else
    for i = 1:n_channels
        Fx = A(x(:, :, i));
        g2 = zeros(size(Fx));
        norm_res{i} = cell(length(G{i}), 1);
        for j = 1:length(G{i})
            r2 = G{i}{j} * (2 * Fx(W{i}{j}) - Fx_old(W{i}{j}, i));
            proj{i}{j} = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
                r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, ...
                elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
            v2{i}{j} = v2{i}{j} + pU{i}{j} .* (r2 -  proj{i}{j});
            clear r2;
            g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * v2{i}{j};

            norm_res{i}{j} = sqrt(sum(abs(G{i}{j} * Fx(W{i}{j}) - y{i}{j}).^2));
            global_norm_res = global_norm_res + norm_res{i}{j}^2;
            norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
        end
        Fx_old(:, i) = Fx; clear Fx;
        Ftx(:, :, i) = sigma22(i) * real(At(g2));
    end
end

end
