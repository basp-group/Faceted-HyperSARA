function [v1, g] = update_l21_spmd(v1, x_overlap, weights, beta1, Iq, dims_q,...
                                  I_overlap_q, dims_overlap_q, offset, ...
                                  status_q, nlevel, wavelet, Ncoefs_q, ...
                                  temLIdxs_q, temRIdxs_q, offsetLq, ...
                                  offsetRq, dims_overlap_ref_q)
%update_l21_spmd: update the dual variable realted to the facet l21-norm
% prior.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > v1                      dual variable associated with the facet l21-
%                           norm [s, L]
% > x_overlap               overlapping image facet [M, N, L]
% > weights                 weights to balance the effect of redundant
%                           pixels due to the overlap between facets [M, N]
% > beta1                   ratio between regularization and convergence 
%                           parameter (gamma1 / sigma1) [1]
% > Iq                      starting index of the non-overlapping base 
%                           facet [1, 2]
% > dims_q                  dimensions of the non-overlapping base facet
%                           [1, 2]
% > I_overlap_q             starting index of the facet [1, 2]
% > dims_overlap_q          dimensions of the facet [1, 2]
% > offset                  offset to be used from one dictionary to
%                           another (different overlap needed for each 
%                           dictionary -> cropping) {nDictionaries}
% > status_q                status of the current facet (last or first 
%                           facet along vert. or hrz. direction)
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs_q                size of the wavelet decompositions at each
%                           scale
% > temLIdxs_q              amount of cropping from the "left" [1, 2]
% > temRIdxs_q              amount of cropping from the "right" [1, 2]
% > offsetLq                amount of zero-pading from the "left" [1, 2]
% > offsetRq                amout of zero-padding from the "right" [1, 2]
% > dims_overlap_ref_q      dimension of the facet [1, 2]    
%
% Output:
%
% < v1                      dual variable associated with the l21-norm 
%                           prior [s, L]
% < g                       auxiliary variable for the update of the primal
%                           variable [M, N, L]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revision: [08/08/2019]
%-------------------------------------------------------------------------%
%%

zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;
g = zeros(size(x_overlap));

for l = 1 : size(x_, 3)
    w = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
    
    w = v1(:, l) +  w;
    l2 = sqrt(sum(abs(w).^2,2));
    l2_soft = max(l2 - beta1*weights, 0)./(l2+eps);
    v1(:, l) = w - l2_soft.*w;
    
    g(:, :, l) = isdwt2_sara(v1(:, l), Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
end

end
