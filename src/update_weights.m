function [weights1, weights0] = update_weights(x_overlap, size_v1, overlap, ...
    I, dims, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha)

% Motivation:
% to be done inside a function (for loops are slow in spmd if not
% encapsulated in a function)

% nuclear norm
sol = reshape(x_overlap(overlap(1)+1:end, overlap(2)+1:end, :), [prod(dims), size(x_overlap, 3)]);
[~,S00,~] = svd(sol,'econ');
d_val0 = abs(diag(S00));
weights0 = reweight_alpha ./ (reweight_alpha + d_val0);

% l21 norm
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap;
w = zeros(size_v1);
for l = 1 : size(x_, 3)
    w(:, l) = sdwt2_sara(x_(:, :, l), I, dims, offset, status, nlevel, wavelet, Ncoefs);
end
d_val1 = sqrt(sum(abs((w)).^2,2));
weights1 = reweight_alpha ./ (reweight_alpha + d_val1);

end
