function y = HS_operatorGtPhi_t(x, H, W, At, Sigma, aW)
% Adjoint of the new reduced measurement operator: At * H * S * Sigma
% Complex -> Real
% ! Attention: H is a reduced holographic matrix by removing unnecessary
% rows
%
% Author: Ming Jiang, E-mail: ming.jiang@epfl.ch
%
if iscell(x)
    c = length(x);
else
    [~, c] = size(x);
end

% Variable flagW for the case where W is present
flagW = 0;
if ~isempty(W)
    flagW = 1;
end

if flagW
    No = size(W{1}{1}, 1);
else
    No = size(H{1}{1}, 2);
end

for ind = 1:c
    if iscell(H{ind})
        x1 = zeros(No, 1);
        for j = 1:length(H{ind})
            if exist('aW', 'var')
                xtmp = sqrt(aW{ind}{j}) .* x{ind}{j};
            else
                xtmp = x{ind}{j};
            end
            if flagW
                x1(W{ind}{j}) = x1(W{ind}{j}) + H{ind}{j}' * (Sigma{ind}{j} .* xtmp(:));
            else
                x1 = x1 + H{ind}{j}' * (Sigma{ind}{j} .* xtmp(:));
            end
        end
    else
        if exist('aW', 'var')
            xtmp = sqrt(aW{ind}) .* x{ind};
        else
            xtmp = x{ind};
        end
        if flagW
            x1 = zeros(size(W{ind}));
            x1(W{ind}) = H{ind}' * (Sigma{ind} .* xtmp(:));
        else
            x1 = H{ind}' * (Sigma{ind} .* xtmp(:));
        end
    end
    y(:,:,ind) = real(At(x1));
end

% spmd
%     if iscell(H)
%         x1 = zeros(No, 1);
%         for j = 1:length(H)
%             if exist('aW', 'var')
%                 xtmp = sqrt(aW{j}) .* x{j};
%             else
%                 xtmp = x{j};
%             end
%             if flagW
%                 x1(W{j}) = x1(W{j}) + H{j}' * (Sigma{j} .* xtmp(:));
%             else
%                 x1 = x1 + H{j}' * (Sigma{j} .* xtmp(:));
%             end
%         end
%     else
%         if exist('aW', 'var')
%             xtmp = sqrt(aW) .* x;
%         else
%             xtmp = x;
%         end
%         if flagW
%             x1 = zeros(size(W));
%             x1(W) = H' * (Sigma .* xtmp(:));
%         else
%             x1 = H' * (Sigma .* xtmp(:));
%         end
%     end
%     y = real(At(x1));
% end

end
