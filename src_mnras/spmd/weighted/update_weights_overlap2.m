function [weights1, weights0] = update_weights_overlap2(x_overlap, size_v1, ...
    I, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, crop_l21, crop_nuclear, w, sig, sig_bar)
% Update the weights of the per facet priors.
%
% Update the weights involved in the reweighting procedure applied to the 
% faceted l21 and nuclear norms. This version includes a spatial weighting 
% correction for the nuclear norm (tapering window).
%
% Args:
%     x_overlap (array): overlapping image facet [M, N, L].
%     size_v1 (int): size of the dual variable associated with the facet 
%                    l21-norm [1, 2].
%     I (array): starting index of the non-overlapping facet [1, 2].
%     offset (array): offset to be used from one dictionary to
%                          another (different overlap needed for each 
%                          dictionary -> cropping) {nDictionaries}.
%     status (array): status of the current facet (last or first 
%                          facet along vert. / hrz. direction).
%     nlevel (int): depth of the wavelet decompositions.
%     wavelet (cell): name of the wavelet dictionaries.
%     Ncoefs (array): size of the wavelet decompositions at each scale.
%     dims_overlap_ref (array): dimension of the facet [1, 2].
%     offsetL (array): amount of zero-pading from the "left" [1, 2].
%     offsetR (array): amount of zero-padding from the "right" [1, 2].
%     reweight_alpha (double): reweighting parameter.
%     crop_l21 (array): relative cropping necessary for the facet 
%                            l21-norm [1, 2].
%     crop_nuclear (array): relative cropping necessary for the facet
%                                nuclear norm [1, 2].
%     w (array): spatial weights applied to the facet nuclear norm 
%                     (same size as x_overlap after cropping by 
%                      crop_nuclear)
%     sig (double): noise level (wavelet space) [1]
%     sig_bar (double): noise level (singular value space) [1]
%
% Returns:
%     weights1 (array): weights for the reweighting of the facet l21-
%                            norm [size_v1(1), size_v1(2)].
%     weights0 (array): weights for the reweighting of the nuclear 
%                            norm [min(M*N, L), 1].

%-------------------------------------------------------------------------%
%%
% Motivation: for loops are slow in spmd if not encapsulated in a function.
% Note: the l21 and nuclear norms do not act on the same facet, hence the
% cropping step.
%%
% Code: P.-A. Thouvenin.
% Last revised: [18/09/2020]
%-------------------------------------------------------------------------%
%%

% nuclear norm
sol = w.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
sol = reshape(sol, [numel(sol)/size(sol, 3), size(x_overlap, 3)]);
[~,S00,~] = svd(sol,'econ');
d0 = abs(diag(S00));
upsilon_bar = sig_bar*reweight_alpha;
weights0 = upsilon_bar ./ (upsilon_bar + d0);

% l21 norm
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
z = zeros(size_v1);
for l = 1 : size(x_, 3)
    z(:, l) = sdwt2_sara_faceting(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
end
d1 = sqrt(sum(z.^2,2));
upsilon = sig*reweight_alpha;
weights1 = upsilon ./ (upsilon + d1);

% parfor k = 1:P
%     d_val = abs(Psit{k}(xsol));
%     res_val = Psit{k}(res_im/peak_psf);
%     noise_val = std(res_val);
%     weights{k} = noise_val * max(reweight_alpha, 1) ./ (noise_val * max(reweight_alpha, 1) + d_val);
% %                 weights{k} = max(reweight_alpha, noise_val) ./ (max(reweight_alpha, noise_val) + d_val);
% %                 kstd = 3*1.4826*mad(res_val);
% %                 weights{k} = 1 ./ (1e-10 + d_val/(reweight_alpha * gamma));
% %                 weights{k}(d_val > max(d_val) * reweight_abs_of_max) = 0;
% end
% 
% reweight_alpha = reweight_alpha_ff * reweight_alpha;


end
