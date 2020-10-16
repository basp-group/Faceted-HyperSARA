function [v2, Ftx, proj, norm_res, norm_residual_check, norm_epsilon_check] = update_data_fidelity_dr_block_new_basic(v2, y, xhat, proj, A, At, H, W, T, Wm, pU, epsilon, ...
                     elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps, sigma22, precondition, reduction_version, realdatablocks)
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
% > W                               mask of non-zero columns of G {L}{nblocks}
%                                   [] for the case W is absent
% > T                               pseudo singular values from the 
%                                   reduction operator {L}{nblocks} (Sigma)
% > Wm                              Fourier masking operator {L}{nblocks}
% > pU                              preconditioning matrices {L}{nblocks}
% > epsilon                         l2-ball radius {L}{nblocks}
% > elipse_proj_max_iter            maximum number of iterations (projection)
%                                   [1]
% > elipse_proj_min_iter            minimum number of iterations (projection)
%                                   [1]
% > elipse_proj_eps                 stopping criterion for the projection [1]
% > sigma22                         update step for the dual variable (tau*sigma2)
% > precondition                    preconditioning [1]
% > reduction_version               option 1: embedding operator F\Phi^t
%                                   option 2: embedding operator G^t
% > realdatablocks                  only for VLA realdata, currently valid
%                                   for either 2 blocks or 9 blocks
%
% Output:
%
% < v2                              data fidelity dual variable {L}{nblocks}[M, 1]
% < Ftx                             auxiliary variable for the update of
%                                   the primal variable [N(1), N(2)]
% < proj                            ... {L}{nblocks}[M, 1]
% < norm_res                        norm of the residual {L}{nblocks}[1]    
% < norm_residual_check_[c,a]       ! only for real data, global square norm of
%                                   the [C,A]-config residual [1]
% < norm_epsilon_check_[c,a]        ! only for real data, global square norm of
%                                   the [C,A]-config epsilon [1]
%
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%% 

FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

% Variable flag for the case where W is present
flagW = 0;
if ~isempty(W)
    flagW = 1;
end

% number of over-sampled pixels
if flagW
    No = size(W{1}{1}, 1);
else
    No = size(H{1}{1}, 2);
end

R = length(H{1});

Ftx = zeros(size(xhat));
[Ny, Nx] = size(xhat(:,:,1));
ci = size(xhat, 3);
norm_res = cell(ci, 1);
% only for real data
norm_residual_check = 0;
norm_epsilon_check = 0;
for i = 1 : ci
    Fx = A(xhat(:,:,i));
    g2 = zeros(No,1);
    for j = 1 : R
        if reduction_version == 1
            if flagW
                tmp = FT2(real(At(H{i}{j} * Fx(W{i}{j}))));
            else
                tmp = FT2(real(At(H{i}{j} * Fx)));
            end
            tmp = tmp(:);
            r2 = T{i}{j} .* tmp(Wm{i}{j});
        elseif reduction_version == 2             
            if flagW
                r2 = T{i}{j} .* (H{i}{j} * Fx(W{i}{j}));
            else                
                r2 = T{i}{j} .* (H{i}{j} * Fx);
            end
        end

        if precondition
            [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, ...
                r2, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
            v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2 - pU{i}{j} .* proj{i}{j};
        else
            v2{i}{j} = v2{i}{j} + r2{i}{j} - y{i}{j} - sc(v2{i}{j} + r2{i}{j} - y{i}{j}, epsilon{i}{j});
        end
        if reduction_version == 1
            tmp = zeros(size(Wm{i}{j}));
            tmp(Wm{i}{j}) = T{i}{j} .* v2{i}{j};
            if flagW
                g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            else
                g2 = g2 + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            end
        elseif reduction_version == 2
            if flagW
                g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j}' * (T{i}{j} .* v2{i}{j});
            else
                g2 = g2 + H{i}{j}' * (T{i}{j} .* v2{i}{j});
            end
        end

        % norm of residual
        norm2 = sum(power(abs(r2 - y{i}{j}), 2));     % faster than built-in norm
        norm_res{i}{j} = sqrt(norm2);

        % Only for real data
        norm_residual_check = norm_residual_check + norm2;
        norm_epsilon_check = norm_epsilon_check + epsilon{i}{j}^2;
    end
    Ftx(:,:,i) = real(At(g2));
end
Ftx = sigma22*Ftx;

end
