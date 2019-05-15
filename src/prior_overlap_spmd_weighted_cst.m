function [l21_norm, nuclear_norm] = prior_overlap_spmd_weighted_cst(x_overlap, Iq, ...
    dims_q, offset, status_q, nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, ...
    offsetLq, offsetRq, w, crop)

% remark: add the weights in the formula? to be determined... (for now, keep as is)

% nuclear norm
c = size(x_overlap, 3);
xhatm = reshape(x_overlap(crop(1)+1:end, crop(2)+1:end).*w,prod(dims_overlap_ref_q),c);
[~,S0,~] = svd(xhatm,'econ');
nuclear_norm = norm(diag(S0),1);

% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;

% compute local l21-norm
l21_norm = 0;
for l = 1 : c
    wl = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
    l2 = sqrt(sum(abs(wl).^2,2));
    l21_norm = l21_norm + norm(l2(:),1);
end

end
