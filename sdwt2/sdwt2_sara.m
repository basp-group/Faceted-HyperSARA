function SPsitLx = sdwt2_sara(x_overlap, I, offset, status, J, wavelet, Ncoefs)
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
szS = size(x_overlap);
M = numel(wavelet);
id = 0;
%id_Ncoefs = 0;

% Compute total number of wavelet coefficients
p = prod(Ncoefs, 2);
s = 3*sum(p(1:end)) - 2*sum(p(J+1:J+1:end)) - 2*p(end); % number of coeffs with the Dirac basis
SPsitLx = zeros(s, 1);

% Renormalize wavelet coefficients (use of several dictionaries)
x_overlap = x_overlap/sqrt(M); % fewer coefs at this stage

% Offsets
LId   = zeros(2, M);
for i = 1:2
    if I(i) > 0
        LId(i,:) = offset.';% no offset if I(q, i) == 0
        LId(i,M) = LId(i,M) + mod(I(i), 2^J); % for Dirac basis
    end
end
LId = LId + 1;

for m = 1:M-1
    %crop_offset = zeros(2, 1); % cropping due to the differences in the
    % length of the wavelet filters (different wavelet transforms)
    
    % forward transform (not Dirac dictionary)
    % define the portion of x_overlap to be exploited for the dictionary considered
    %x_overlap_tmp = x_overlap(crop_offset(1)+1:end, crop_offset(2)+1:end);
    %[lo, hi] = wfilters(wavelet{m}, 'd'); % decomposition filters
    sm = 3*sum(p((m-1)*(J+1)+1:m*(J+1)-1)) + p(m*(J+1));    
    SPsitLx(id+1:id+sm) = sdwt2(x_overlap(LId(1,m):szS(1),LId(2,m):szS(2)), status, wavelet{m}, J, Ncoefs((m-1)*(J+1)+1:m*(J+1),:));
    id = id + sm;
end

SPsitLx(id+1:end) = col(x_overlap(LId(1,M):end,LId(2,M):end));

end