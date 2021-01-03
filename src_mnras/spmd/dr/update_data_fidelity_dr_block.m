function [v2, Ftx, proj, norm_res, global_norm_res, norm_epsilon] = update_data_fidelity_dr_block(v2, y, xhat, proj, A, At, H, T, W, pU, epsilon, ...
                     elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps, sigma22)
%update_data_fidelity_dr_block: update step of the data fidelity term (with 
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
% > proj                            value of the projection at the previous
%                                   global iteration, taken as a starting 
%                                   point for the new projection step {L}{nblocks}[M, 1]
% > A                               measurement operator @[1]
% > At                              adjoint measurement operator @[1]
% > H                               holographic matric G'*G {L}
% > T                               pseudo singular values from the 
%                                   reduction operator {L}{nblocks} (Sigma)
% > W                               Fourier masking operator {L}{nblocks}
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
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%% 

FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

Ftx = zeros(size(xhat));
[Ny, Nx] = size(xhat(:,:,1));
N = Ny * Nx;
ci = size(xhat, 3);
norm_res = cell(ci, 1);
global_norm_res = 0;
norm_epsilon = 0;
for i = 1 : ci
    Fx0 = A(real(xhat(:,:,i)));
    H_mat = sparse(size(H{i}{1}, 1), size(H{i}{1}, 2));
    g2 = zeros(N, 1);
    norm_res{i} = cell(length(T{i}), 1);
    for j = 1:length(H{i})
        Fx = FT2(real(At(H{i}{j} * Fx0)));
        Fx = Fx(:);
        r2 = T{i}{j} .* Fx(W{i}{j});
        proj{i}{j} = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
            r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
        v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2 - pU{i}{j} .* proj{i}{j};        
        u2 = T{i}{j} .* v2{i}{j};
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
        H_mat = H{i}{j} + H_mat;

        norm_res{i}{j} = norm(r2 - y{i}{j}, 2);
        global_norm_res = global_norm_res + norm_res{i}{j}^2;
        norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
    end
    Ftx(:,:,i) = real(At(H_mat*A(real(IFT2(reshape(g2, Ny, Nx))))));
end
Ftx = sigma22*Ftx;

end
