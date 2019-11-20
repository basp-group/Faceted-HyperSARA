function y = HS_operatorGtPhi(x, H, W, A, Sigma, aW)
% New reduced measurement operator: Sigma * S * H * A
% Real -> Complex
% ! Attention: H is a reduced holographic matrix by removing unnecessary
% rows

c = size(x, 3);

% Variable flagW for the case where W is present
flagW = 0;
if ~isempty(W)
    flagW = 1;
end

for ind = 1:c
    x1 = A(real(x(:,:,ind)));
    if iscell(H{ind})
        for j = 1:length(H{ind})
            if flagW
                x2 = H{ind}{j} * x1(W{ind}{j});
            else
                x2 = H{ind}{j} * x1;
            end
            xtmp = Sigma{ind}{j} .* x2;
            if exist('aW', 'var')
                y{ind}{j} =  sqrt(aW{ind}{j}) .* xtmp;
            else
                y{ind}{j} =  xtmp;
            end
        end
    else
        if flagW
            x2 = H{ind} * x1(W{ind});
        else
            x2 = H{ind} * x1;
        end
        xtmp = Sigma{ind} .* x2;
        if exist('aW', 'var')
            y{ind} =  sqrt(aW{ind}) .* xtmp;
        else
            y{ind} =  xtmp;
        end
    end
end

end
