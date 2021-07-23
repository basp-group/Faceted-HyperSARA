function y = update_l21_spmd_debug(x_overlap, Iq, dims_q, I_overlap_q, dims_overlap_q, offset, status_q, nlevel, wavelet, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q)
% Compute the direct and inverse discrete wavelet transform for a
% single facet (for debugging purpose).
%%

zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding (to be checked again...)
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;
y = zeros(size(x_overlap));

for l = 1 : size(x_, 3)
    w = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
    
    y(:, :, l) = isdwt2_sara(w, Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
end
