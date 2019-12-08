function rg = domain_decomposition(nchunks, N)
% Tessellates 1:N into non-overlapping subsets.
%
% Tessellates 1:N into non-overlapping subsets, each containing
% approximately the same number of indices (overlap from the left).
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
splits = round(linspace(0, N, nchunks + 1));
% rg = zeros(nchunks, 2);
% for q = 1:nchunks
%    rg(q, :) = [(splits(q) + 1), splits(q + 1)]; 
% end
rg = [splits(1:end-1).'+1, splits(2:end).'];

end
