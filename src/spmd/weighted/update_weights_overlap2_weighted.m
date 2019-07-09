function [weights1, weights0] = update_weights_overlap2_weighted(x_overlap, size_v1, ...
    I, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, crop_l21, crop_nuclear, w)
%update_weights_overlap2_weighted: update the weights involved in the 
% reweighting procedure applied to the faceted l21 and nuclear norms. This
% version includes the spatial weighting correction for the nuclear norm.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x_overlap               overlapping image facet [M, N, L]
% > size_v1                 size of the dual variable associated with the 
%                           facet l21-norm [1, 2]
% > I                       starting index of the non-overlapping facet 
%                           [1, 2]
% > offset                  offset to be used from one dictionary to
%                           another (different overlap needed for each 
%                           dictionary -> cropping) {nDictionaries}
% > status                  status of the current facet (vert. and hrz) ...
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs                  size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref        dimension of the facet [1, 2] 
% > offsetL                 amount of zero-pading from the "left" [1, 2]
% > offsetR                 amout of zero-padding from the "right" [1, 2]
% > reweight_alpha          reweighting parameter [1]
% > crop_l21                relative cropping necessary for the facet l21- 
%                           norm [1, 2]
% > crop_nuclear            relative cropping necessary for the facet
%                           nuclear norm [1, 2]
% > w                       spatial weights applied to the facet nuclear
%                           norm (same size as x_overlap after cropping by 
%                           crop_nuclear)
%
% Note: the l21 and nuclear norms do not act on the same facet, hence the
% cropping step
%
% Output:
%
% < weights1                weights for the reweighting of the facet l21-
%                           norm [size_v1(1), size_v1(2)]
% < weights0                weights for the reweighting of the nuclear norm
%                           [min(M*N, L), 1]
%-------------------------------------------------------------------------%
%%
% Motivation: for loops are slow in spmd if not encapsulated in a function.
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

% nuclear norm
sol = w.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
sol = reshape(sol, [numel(sol)/size(sol, 3), size(x_overlap, 3)]);
[~,S00,~] = svd(sol,'econ');
d0 = abs(diag(S00));
weights0 = reweight_alpha ./ (reweight_alpha + d0);

% l21 norm
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
z = zeros(size_v1);
for l = 1 : size(x_, 3)
    z(:, l) = sdwt2_sara(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
end
d1 = sqrt(sum(z.^2,2));
weights1 = reweight_alpha ./ (reweight_alpha + d1);

end
