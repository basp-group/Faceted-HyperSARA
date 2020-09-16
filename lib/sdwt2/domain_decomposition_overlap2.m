function rg = domain_decomposition_overlap2(nchunks, N, d)
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
% sanity check
if (d > floor(N/nchunks))
   error('Overlap greater than the dimension of a single non-redundant facet');
end

% only extend the facet toward the left
splits = round(linspace(0, N, nchunks + 1));

% rg = zeros(nchunks, 2);
% for q = 1:nchunks
%     if q > 1
%         id_start = (splits(q) + 1) - d;
%     else
%         id_start = (splits(q) + 1);
%     end
% 
%     id_end = splits(q + 1);
% 
%     rg(q, :) = [id_start, id_end];
% end

rg = [vertcat(splits(1)+1, splits(2:end-1).'+1-d), splits(2:end).'];

end
