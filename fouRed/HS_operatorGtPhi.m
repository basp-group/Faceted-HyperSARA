function y = HS_operatorGtPhi(x, H, A, Sigma, Mask, aW)
% New reduced measurement operator: Sigma * S * H * A
% Real -> Complex

% Parameters
% c = size(x, 3);

% for ind = 1:c
%     x1 = A(real(x(:,:,ind)));
%     
%     if iscell(H{ind})
%         for j = 1:length(H{ind})
%             dimH = sqrt(numel(H{ind}{j}));
%             x2 = reshape(H{ind}{j}, dimH, dimH) * x1;
%             xtmp = Sigma{ind}{j} .* x2(Mask{ind}{j});
%             if exist('aW', 'var')
%                 y{ind}{j} =  sqrt(aW{ind}{j}) .* xtmp;
%             else
%                 y{ind}{j} =  xtmp;
%             end
%         end
%     else    
%         dimH = sqrt(numel(H{ind}));
%         x2 = reshape(H{ind}, dimH, dimH) * x1;
%         xtmp = cell2mat(Sigma{ind}) .* x2(cell2mat(Mask{ind}));
%         if exist('aW', 'var')
%             y(:,ind) =  sqrt(cell2mat(aW{ind})) .* xtmp;
%         else
%             y(:,ind) =  xtmp;
%         end
%     end
%     
% end

x1 = A(x);

if iscell(H)
    for j = 1:length(H)
        dimH = sqrt(numel(H{j}));
        x2 = reshape(H{j}, dimH, dimH) * x1;
        xtmp = Sigma{j} .* x2(Mask{j});
        if exist('aW', 'var')
            y{j} =  sqrt(aW{j}) .* xtmp;
        else
            y{j} =  xtmp;
        end
    end
else    
    dimH = sqrt(numel(H));
    x2 = reshape(H, dimH, dimH) * x1;
    xtmp = cell2mat(Sigma) .* x2(cell2mat(Mask));
    if exist('aW', 'var')
        y =  sqrt(cell2mat(aW)) .* xtmp;
    else
        y =  xtmp;
    end
end

end
