function rg = domain_decomposition(nchunks, N)
% Tessellates 1:N into non-overlapping subsets, each containing
% approximately the same number of indices.
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
splits = round(linspace(0, N, nchunks + 1));
rg = zeros(nchunks, 2);

for q = 1:nchunks
   rg(q, :) = [(splits(q) + 1), splits(q + 1)]; 
end

end
