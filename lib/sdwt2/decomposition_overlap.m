function rg = decomposition_overlap(nchunks, N, d)
% Tessellates 1:N into overlapping subsets.
%
% Tessellates 1:N into overlapping subsets, each containing
% approximately the same number of indices (overlap from the left).
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

rg = zeros(nchunks, 2);
rg(1, :) = [splits(1)+1, splits(2)];
for q = 2:nchunks
   rg(q, :) = [(splits(q) + 1) - d, splits(q + 1)]; 
end

end
