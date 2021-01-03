function rg = interleaved_facets(nchunks, N)
% Tessellates 1:N into interleaved subsets (subsampling).
%
% Tessellates 1:N into subsets of interleaved indices, each containing
% approximately the same number of indices (downsampling of 1:N)
%
% Args:
%     nchunks (int): number of output segments.
%     N (int): total number of indices.
%
% Returns:
%     rg (array_like): first/last index of each segment [nchunks, 2].
% 

%-------------------------------------------------------------------------%
%%
if nchunks > N
    error('Number of facets Q=%i greater than the dimension L=%i', ...
    nchunks, N);
end

rg = cell(nchunks, 1);
for q = 1:nchunks
    rg{q} = q:nchunks:N;
end

end
