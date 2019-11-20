function y = HS_fouRed_forward_operator_new(x, A, At, H, W, Sigma, Mask, aW)
% Phi' = Sigma * F * Phi^t * Phi = Sigma * F * A^t * H * A
% Real -> Complex
% ! Attention: H is a self-adjoint complete holographic matrix, or the 
% number of rows is equal to oversampled image size.

FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));

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
        if flagW
            x2 = H{ind} * x1(W{ind});
        else
            x2 = H{ind} * x1;
        end
        xtmp = FT2(real(At(x2)));
        xtmp = xtmp(:);
        xtmp = Sigma{ind} .* xtmp(Mask{ind});
        if exist('aW', 'var')
            y{ind} =  sqrt(aW{ind}) .* xtmp;
        else
            y{ind} =  xtmp;
        end
    end
    
end

end
