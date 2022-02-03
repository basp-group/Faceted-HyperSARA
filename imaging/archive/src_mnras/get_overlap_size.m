function overlapSize = get_overlap_size(N, Q, overlapFraction)
% Convert fraction of overlap into the corresponding number of pixels.
%
% Compute the number of pixels corresponding to a prescribed overlap 
% fraction.
%
% Parameters
% ----------
% N : array, int
%     Size of the full image [1, 2].
% Q : array, int
%     Number of facets along each dimension [1, 2].
% overlapFraction : array, double
%     Overlap fraction between two consecutive facets, given along each
%     dimension (fraction expressed with respect to the final facet size).
%
% Raises
% ------
% AssertionError
%     All entries in ``overlapFraction`` need to be strictly lower than 0.5
%     (only overlap between consecutive facets is supported).
%
% Returns
% -------
% overlapSize : array, int
%     Number of pixels contained in the overlap region.
%

if any(overlapFraction > 1)
    error('get_overlap_size:InputValueError', 'All entries in overlapFraction need to be < 0.5')
end

overlapSize = floor((Q>1).*(N./Q).*(overlapFraction./(1-overlapFraction)));

end