function rg = domain_decomposition_overlap(nchunks, N, d)
% Tessellates 1:N into overlapping subsets, each containing
% approximately the same number of indices.
%-------------------------------------------------------------------------%
%%
% Input: 
% > nchunks  number of output segments [1]
% > N        total number of indices [1]
% < d        number of indices in the overlap
%
% Output:
% < rg       first/last index of each segment [nchunks, 2]
%-------------------------------------------------------------------------%
%%
splits = round(linspace(0, N, nchunks + 1));

r = rem(d, 2);
rg = zeros(nchunks, 2);
for q = 1:nchunks
    if q > 1
        id_start = (splits(q) + 1) - floor(d/2) - r;
    else
        id_start = (splits(q) + 1);
    end

    if q < nchunks
        id_end = splits(q + 1) + floor(d/2);
    else
        id_end = splits(q + 1);
    end

    rg(q, :) = [id_start, id_end];
end

end
