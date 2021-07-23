function overlapSize = get_overlap_size(N, Q, overlapFraction)


if any(overlapFraction > 1)
    error('get_overlap_size:InputValueError', 'All entries in overlapFraction need to be < 0.5')
end

overlapSize = floor((Q>1).*(N./Q).*(overlapFraction./(1-overlapFraction)));

end