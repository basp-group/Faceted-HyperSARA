function [u, y, l21_] = run_par_l21_spmd(u, x_overlap, weights, beta1, Iq, dims_q, I_overlap_q, dims_overlap_q, offset, status_q, nlevel, wavelet, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q)

zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding (to be checked again...)
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;
y = zeros(size(x_overlap));

l21_ = 0;
for l = 1 : size(x_, 3)
    u_old = u(:, l);
    w = sdwt2_sara(x_(:, :, l), Iq, dims_q, offset, status_q, nlevel, wavelet, Ncoefs_q);
    w = u(:, l) +  w;
    l2 = sqrt(sum(abs(w).^2,2));
    l2_soft = max(l2 - beta1*weights(:, l), 0)./(l2+eps);
    u(:, l) = w - l2_soft.*w;
    
    y(:, :, l) = isdwt2_sara(u(:, l)-u_old, Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
    
    % local L21 norm of current solution
    l21_ = l21_ + norm(l2(:),1);
end
