function [weights1, weights0] = update_weights_cst_overlap_weighted(x_overlap, size_v1, ...
    I, dims, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, w, crop, sig, sig_bar)
%update_weights_cst_overlap_weighted: update the weights involved in the 
% reweighting procedure applied to the faceted l21 and nuclear norms. This 
% version includes a spatial correction for the faceted nuclear norm 
% (tapering window). Assumption: overlap for the faceted wavelet transforms
% larger than for the nuclear norm.
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
%                           facet along vert. or hrz. direction)
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs                  size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref        dimension of the facet [1, 2] 
% > offsetL                 amount of zero-pading from the "left" [1, 2]
% > offsetR                 amout of zero-padding from the "right" [1, 2]
% > reweight_alpha          reweighting parameter [1]
% > w                       spatial weights applied to the facet nuclear
%                           norm [size(x_overlap)] 
% > crop                    relative cropping necessary for the nuclear 
%                           norm (smaller overlap than the l21) [1, 2]
% > sig                     noise level (wavelet space) [1]
% > sig_bar                 noise level (singular value space) [1]
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
% Last revised: [18/09/2020]
%-------------------------------------------------------------------------%
%%

% nuclear norm
sol = x_overlap(crop(1)+1:end, crop(2)+1:end, :).*w;
sol = reshape(sol, [numel(sol)/size(x_overlap, 3), size(x_overlap, 3)]);
[~,S00,~] = svd(sol,'econ');
d0 = abs(diag(S00));
upsilon_bar = sig_bar*reweight_alpha;
weights0 = upsilon_bar ./ (upsilon_bar + d0);

% l21 norm
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap;
u = zeros(size_v1);
for l = 1 : size(x_, 3)
    u(:, l) = sdwt2_sara(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
end
d1 = sqrt(sum(abs((u)).^2,2));
upsilon = sig*reweight_alpha;
weights1 = upsilon ./ (upsilon + d1);

end
