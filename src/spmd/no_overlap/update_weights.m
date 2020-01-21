function [weights1, weights0] = update_weights(x_overlap, size_v1, overlap, ...
    I, dims, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha)
%update_weights: update the weights involved in the reweighting procedure 
% applied to the faceted l21 and nuclear norms.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x_overlap               overlapping image facet [M, N, L]
% > size_v1                 size of the dual variable associated with the 
%                           facet l21-norm [1, 2]
% > I                       starting index of the non-overlapping facet 
%                           [1, 2]
% > dims                    size of the non-overlapping facet [1, 2] 
% > offset                  offset to be used from one dictionary to
%                           another (different overlap needed for each 
%                           dictionary -> cropping) {nDictionaries}
% > status                  status of the current facet (last or first 
%                           facet along vert. / hrz. direction)
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs                  size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref        size of the facet [1, 2] 
% > offsetL                 amount of zero-pading from the "left" [1, 2]
% > offsetR                 amout of zero-padding from the "right" [1, 2]
% > reweight_alpha          reweighting parameter [1]
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
% Last revision: [08/08/2019]
%-------------------------------------------------------------------------%
%%

% nuclear norm
sol = reshape(x_overlap(overlap(1)+1:end, overlap(2)+1:end, :), [prod(dims), size(x_overlap, 3)]);
[~,S00,~] = svd(sol,'econ');
d0 = abs(diag(S00));
weights0 = reweight_alpha ./ (reweight_alpha + d0);

% l21 norm
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap;
w = zeros(size_v1);
for l = 1 : size(x_, 3)
    w(:, l) = sdwt2_sara(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
end
d1 = sqrt(sum(abs((w)).^2,2));
weights1 = reweight_alpha ./ (reweight_alpha + d1);

end
