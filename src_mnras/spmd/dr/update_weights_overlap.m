function [weights1, weights0] = update_weights_overlap(x_overlap, size_v1, ...
    I, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, sig, sig_bar)
%update_weights_overlap: update the weights involved in the reweighting
% procedure applied to the faceted l21 and nuclear norms.
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
% > status                  status of the current facet (last or first 
%                           facet along vert. / hrz. direction)
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs                  size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref        dimension of the facet [1, 2] 
% > offsetL                 amount of zero-pading from the "left" [1, 2]
% > offsetR                 amout of zero-padding from the "right" [1, 2]
% > reweight_alpha          reweighting parameter [1]
% > sig                     noise level (wavelet space) [1]
% > sig_bar                 noise level (singular value space) [1]
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
% Last revised: [18/09/2020]
%-------------------------------------------------------------------------%
%%

% Update weights for the nuclear norm (one facet)
sol = reshape(x_overlap, [prod(dims_overlap_ref), size(x_overlap, 3)]);
[~, S00, ~] = svd(sol, 'econ');
d0 = abs(diag(S00));
upsilon_bar = sig_bar*reweight_alpha;
weights0 = upsilon_bar ./ (upsilon_bar + d0);

% Update weights for the l21-norm (one facet)
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap;
z = zeros(size_v1);
for l = 1 : size(x_, 3)
    z(:, l) = sdwt2_sara(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
end
d1 = sqrt(sum(abs((z)).^2, 2));
upsilon = sig*reweight_alpha;
weights1 = upsilon ./ (upsilon + d1);

end
