function [SPsitLx, Ij, dims_PsitLx, Ncoefs] = sdwt2_sara(x_overlap, I, dims, offset, status, J, wavelet)
% Forward operator to compute the contribution of the SARA prior relative
% to a single facet.
%%
%-------------------------------------------------------------------------%
% Input:
% > x_overlap  ... 
% > I          starting index of the facet (without overlap) [1, 2] 
% > dims       facet dimension (without overlap) [1, 2]
% > offset     offset to be considered for each dictionary w.r.t. the largest overlapping
%              facet x_overlap [M, 1]
% > status     indicates whether the facet is the first (or the last)
%              along each dimension [1, 2]
% > J          number of scales considered
% > wavelet    name of the wavelets considered {1, M}
%
% Output:
% < SPsitLx
% < Ij
% < dims_PsitLx
% < Ncoefs
%-------------------------------------------------------------------------%
%%
% Debug: 
% [23/10/18] ok.
%-------------------------------------------------------------------------%
%%
M = numel(wavelet);
SPsitLx = []; %cell(M, 1);
Ij = zeros(M*(J+1), 2);
dims_PsitLx = zeros(M*(J+1), 2);
Ncoefs = zeros(M*(J+1), 2);
dim = numel(I);

for m = 1:M   
    temLIdxs = zeros(2, 1);
    
    % forward transform
    if ~strcmp(wavelet{m}, 'self')
        % define the portion of x_overlap to be exploited for the dictionary considered
        for i = 1:dim
            if I(i) > 0
                temLIdxs(i) = offset(m); % no offset if I(q, i) == 0
            end
        end
        x_overlap_tmp = x_overlap(temLIdxs(1)+1:end, temLIdxs(2)+1:end);
        [lo, hi] = wfilters(wavelet{m}, 'd'); % decomposition filters
        [SPsitLx_m, Ij((m-1)*(J+1)+1:m*(J+1),:), dims_PsitLx((m-1)*(J+1)+1:m*(J+1),:), Ncoefs((m-1)*(J+1)+1:m*(J+1),:)] = sdwt2(x_overlap_tmp, I, dims, status, lo, hi, J, M);
    else
        % define the portion of x_overlap to be exploited for the dictionary considered
        for i = 1:dim
            if I(i) > 0
                temLIdxs(i) = offset(m) + mod(I(i), 2^J); % no offset if I(q, i) == 0 (see if this was the missing part)
            end
        end
        x_overlap_tmp = x_overlap(temLIdxs(1)+1:end, temLIdxs(2)+1:end);
        SPsitLx_m = x_overlap_tmp(:)/sqrt(M);
        dims_PsitLx((m-1)*(J+1)+1:m*(J+1),:) = ones(J+1,1)*size(x_overlap_tmp);
    end
    SPsitLx = [SPsitLx; SPsitLx_m];
end

end