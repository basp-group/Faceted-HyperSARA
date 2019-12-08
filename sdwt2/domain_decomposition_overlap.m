function rg = domain_decomposition_overlap(nchunks, N, d)
% Tessellates 1:N into overlapping subsets.
%
% Tessellates 1:N into overlapping subsets, each containing
% approximately the same number of indices (overlap by ca. d/2 on both 
% sides).
%
% Args:
%     nchunks (int): number of output segments.
%     N (int): total number of indices.
%     d (int): number of indices in the overlap.
%
% Returns:
%     rg (array_like): first/last index of each segment [nchunks, 2].
% 

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
