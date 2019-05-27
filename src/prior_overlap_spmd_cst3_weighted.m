function [l21_norm, nuclear_norm] = prior_overlap_spmd_cst3_weighted(x_overlap, Iq, ...
    offset, status_q, nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, ...
    offsetLq, offsetRq, crop_l21, crop_nuclear, w)

% nuclear norm
c = size(x_overlap, 3);
xhatm = w.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
xhatm = reshape(xhatm,numel(xhatm)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear_norm = norm(diag(S0),1);

% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);

% compute local l21-norm
l21_norm = 0;
for l = 1 : c
    wl = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
    l2 = sqrt(sum(abs(wl).^2,2));
    l21_norm = l21_norm + norm(l2(:),1);
end

end
