function rg = interleaved_facets(nchunks, N)
% Tessellates 1:L into subsets of interleaved indices, each containing
% approximately the same number of indices (downsampling of 1:N)
%-------------------------------------------------------------------------%
%%
% Input: 
% > nchunks  number of output segments [1]
% > N        total number of indices [1]
%
% Output:
% < rg       first/last index of each segment [nchunks, 2]
%-------------------------------------------------------------------------%
%%
if nchunks > N
    error('Number of facets Q=%i greater than the dimension L=%i', nchunks,
           N);
end

rg = cell(nchunks, 1);
for q = 1:nchunks
    rg{q} = q:nchunks:N;
end

end
