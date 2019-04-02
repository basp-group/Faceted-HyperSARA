function SPsitLx = sdwt2_sara(x_overlap, I, dims, offset, status, J, wavelet, Ncoefs)
% Forward operator to compute the contribution of the SARA prior relative
% to a single facet.
%%
%-------------------------------------------------------------------------%
% Input:
%
% > x_overlap  redundant image facet 
% > I          starting index of the facet (without overlap) [1, 2] 
% > dims       facet dimension (without overlap) [1, 2]
% > offset     offset to be considered for each dictionary w.r.t. the 
%              largest overlapping facet x_overlap [M, 1]
% > status     indicates whether the facet is the first (or the last)
%              along each dimension [1, 2]
% > J          number of decomposition levels considered
% > wavelet    name of the wavelets considered {1, M}
% > Ncoefs     number of wavelet coefficients for each scale, corresponding
%              to each dictionary involved in the transform [M*(J+1), 1]
%
% Output:
%
% < SPsitLx    resulting wavelet coefficients [write details on the order 
%              of the coefficients inside the vector]
%-------------------------------------------------------------------------%
%%
% Debug: 
% [23/10/18] ok.
% [18/03/19] further code acceleration.
%-------------------------------------------------------------------------%
%%
% Auxiliary variables
M = numel(wavelet);
dim = numel(I);
id = 0;
id_Ncoefs = 0;

% Compute total number of wavelet coefficients
p = prod(Ncoefs, 2);
dirac_present = any(ismember(wavelet, 'self'));
if dirac_present
    s = 3*sum(p(1:end)) - 2*sum(p(J+1:J+1:end)) + prod(dims);
else
    s = 3*sum(p) - 2*sum(p(J+1:J+1:end));
end
SPsitLx = zeros(s, 1);

% Renormalize wavelet coefficients (use of several dictionaries)
x_overlap = x_overlap/sqrt(M); % fewer coefs at this stage, faster then

for m = 1:M   
    crop_offset = zeros(2, 1); % cropping due to the differences in the 
    % length of the wavelet filters (different wavelet transforms)
    
    % forward transform
    if ~strcmp(wavelet{m}, 'self')
        % define the portion of x_overlap to be exploited for the dictionary considered
        for i = 1:dim
            if I(i) > 0
                crop_offset(i) = offset(m); % no offset if I(q, i) == 0
            end
        end
        x_overlap_tmp = x_overlap(crop_offset(1)+1:end, crop_offset(2)+1:end);
        [lo, hi] = wfilters(wavelet{m}, 'd'); % decomposition filters
        sm = 3*sum(p((m-1)*(J+1)+1:m*(J+1)-1)) + p(m*(J+1)); 
        SPsitLx_m = sdwt2(x_overlap_tmp, dims, status, lo, hi, J, Ncoefs(id_Ncoefs+1:id_Ncoefs+(J+1),:));
        SPsitLx(id+1:id+sm) = SPsitLx_m;
        id = id + sm;
        id_Ncoefs = id_Ncoefs + (J+1);
    else
        % define the portion of x_overlap to be exploited for the
        % dictionary considered (remove the overlap from the facet)
        for i = 1:dim
            if I(i) > 0
                crop_offset(i) = offset(m) + mod(I(i), 2^J); % no offset if I(q, i) == 0 (see if this was the missing part)
            end
        end
        x_overlap_tmp = x_overlap(crop_offset(1)+1:end, crop_offset(2)+1:end);
        SPsitLx_m = x_overlap_tmp(:);
        sm = numel(x_overlap_tmp);
        SPsitLx(id+1:id+sm) = SPsitLx_m;
        id = id + sm;
        id_Ncoefs = id_Ncoefs + 1;
    end
end

end