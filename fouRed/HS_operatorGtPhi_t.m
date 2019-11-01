function y = HS_operatorGtPhi_t(x, H, At, Sigma, Mask, Nd, Kd, aW)
% Adjoint of the new reduced measurement operator: At * H * S * Sigma
% Complex -> Real
% Parameters

% if iscell(x)
%     c = numel(x);
% else
%     [~, c] = size(x);
% end

Ny = Nd(1);
Nx = Nd(2);

y = zeros(Ny, Nx);
%
% for ind = 1:c
%     if iscell(H{ind})
%         % old blocking structure
%         for j = 1:length(H{ind})
%             if exist('aW', 'var')
%                 xtmp = sqrt(aW{ind}{j}) .* x{ind}{j};
%             else
%                 xtmp = x{ind}{j};
%             end
%             x1 = zeros(Kd(1) * Kd(2), 1);
%             x1(Mask{ind}{j}) = Sigma{ind}{j} .* xtmp(:);
%             dimH = sqrt(numel(H{ind}{j}));
%             y(:,:,ind) = y(:,:,ind) + real(At(reshape(H{ind}{j}, dimH, dimH) * x1));
%         end
%     else
%         if exist('aW', 'var')
%             xtmp = sqrt(cell2mat(aW{ind})) .* x(:,ind);
%         else
%             xtmp = x(:,ind);
%         end
%         x1 = zeros(Kd(1) * Kd(2), 1);
%         x1(cell2mat(Mask{ind})) = cell2mat(Sigma{ind}) .* xtmp(:);
%         dimH = sqrt(numel(H{ind}));
%         y(:,:,ind) = real(At(reshape(H{ind}, dimH, dimH) * x1));
%     end
% end

if iscell(H)
    % old blocking structure
    for j = 1:length(H)
        if exist('aW', 'var')
            xtmp = sqrt(aW{j}) .* x{j};
        else
            xtmp = x{j};
        end
        x1 = zeros(Kd(1) * Kd(2), 1);
        x1(Mask{j}) = Sigma{j} .* xtmp(:);
        y = y + real(At(H{j} * x1));
    end
else
    if exist('aW', 'var')
        xtmp = sqrt(cell2mat(aW)) .* x;
    else
        xtmp = x;
    end
    x1 = zeros(Kd(1) * Kd(2), 1);
    x1(cell2mat(Mask)) = cell2mat(Sigma) .* xtmp(:);
    y = real(At(H * x1));
end

end
