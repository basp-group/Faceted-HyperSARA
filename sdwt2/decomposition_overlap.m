function rg = decomposition_overlap(nchunks, N, d)
% Tessellates 1:N into overlapping subsets, each containing
% approximately the same number of indices (overlap from the left).
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

rg = zeros(nchunks, 2);
rg(1, :) = [splits(1)+1, splits(2)];
for q = 2:nchunks
   rg(q, :) = [(splits(q) + 1) - d, splits(q + 1)]; 
end

end
