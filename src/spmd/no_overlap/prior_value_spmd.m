function [l21_norm, nuclear_norm] = prior_value_spmd(x_overlap, crop, Iq, ...
    offset, status_q, nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, ...
    offsetLq, offsetRq, size_v1)
%prior_value_spmd: compute the value of the faceted prior (l21 + nuclear 
% norm) for a single facet.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x_overlap               overlapping image facet [M, N, L]
% > crop                    relative cropping necessary for the nuclear 
%                           norm (no overlap for this term) [1, 2]
% > Iq                      starting index of the non-overlapping facet 
%                           [1, 2]
% > offset                  offset to be used from one dictionary to
%                           another (different overlap needed for each 
%                           dictionary -> cropping) {nDictionaries}
% > status_q                status of the current facet (vert. and hrz) ...
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs_q                size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref_q      dimension of the facet [1, 2] 
% > offsetLq                amount of zero-pading from the "left" [1, 2]
% > offsetRq                amout of zero-padding from the "right" [1, 2]
% > size_v1                 size of the dual variable associated with the 
%                           facet l21-norm [1, 2]
%
% Output:
%
% < l21_norm                facet l21-norm [1]
% < nuclear_norm            facet nuclear norm [1]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%% compute facet nuclear norm
c = size(x_overlap, 3);
xhatm = x_overlap(1+crop(1):end, 1+crop(2):end, :);
xhatm = reshape(xhatm,numel(xhatm)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear_norm = norm(diag(S0),1);

%% compute facet l21-norm
% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;

% faceted SARA
z = zeros(size_v1);
for l = 1 : c
    z(:,l) = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
end
l21_norm = sum(sqrt(sum(abs(z).^2,2)), 1);

end
