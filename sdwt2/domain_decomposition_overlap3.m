function rg = domain_decomposition_overlap3(nchunks, N, p)
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
if (p > 0.5)
   error('Overlap fraction greater than 50% of a given facet'); 
   % need to change the communication scheme in this case?
end

% only extend the facet towards the left
% splits = round(linspace(0, N, nchunks + 1));
% rg = zeros(nchunks, 2);
% for q = 1:nchunks
%     if q > 1
%         id_start = (splits(q) + 1) - floor(p*(splits(q + 1) - splits(q)));
%     else
%         id_start = (splits(q) + 1);
%     end
% 
%     id_end = splits(q + 1);
% 
%     rg(q, :) = [id_start, id_end];
% end

% splits = round(linspace(0, N, nchunks + 1));
% rg = zeros(nchunks, 2);
% for q = 1:nchunks
%    rg(q, :) = [(splits(q) + 1), splits(q + 1)]; 
% end
% rg(2:end,:) = rg(2:end,:) - floor(p.*(rg(1:end-1,2)-rg(1:end-1,1)+1)).*[1:nchunks-1].'; % to be checked
% 
% if (rg(end, 2) - rg(end, 1) + 1 >= floor(N/nchunks))
%     rg = [rg; [rg(end,1)-floor(p.*(rg(end,2)-rg(end,1)+1)),N]];
% end

splits = round(linspace(0, N, floor(nchunks/p) + 1));
rg = [splits(1:end-floor(1/p)).'+1, splits(floor(1/p)+1:end).'];
% rg = [splits(1:end-1).'+1, splits(2:end).']; % no redundancy

% see if an additional chunk is needed or not

end
