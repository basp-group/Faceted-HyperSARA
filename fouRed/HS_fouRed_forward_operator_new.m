function y = HS_fouRed_forward_operator_new(x, A, At, H, Sigma, Mask, aW)
% Phi' = Sigma * F * Phi^t * Phi = Sigma * F * A^t * H * A

FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));

% Parameters
[Ny, Nx, c] = size(x);

%
for ind = 1:c
    x1 = A(real(x(:,:,ind)));
    
    if iscell(H{ind})
        for j = 1:length(H{ind})
            x2 = H{ind}{j} * x1;
            xtmp = FT2(real(At(x2)));
            xtmp = xtmp(:);
            xtmp = Sigma{ind}{j} .* xtmp(Mask{ind}{j});
            if exist('aW', 'var')
                y{ind}{j} =  sqrt(aW{ind}{j}) .* xtmp;
            else
                y{ind}{j} =  xtmp;
            end
        end
    else    
        x2 = H{ind} * x1;
        xtmp = FT2(real(At(x2)));
        xtmp = xtmp(:);
        xtmp = cell2mat(Sigma{ind}) .* xtmp(cell2mat(Mask{ind}));
        if exist('aW', 'var')
            y(:,ind) =  sqrt(cell2mat(aW{ind})) .* xtmp;
    %         y(:,ind) =  sqrt(cell2mat(aW{ind})) .* B{ind}(x(:,:,ind));
        else
            y(:,ind) =  xtmp;
    %         y(:,ind) =  B{ind}(x(:,:,ind));
        end
    end
    
end

end
