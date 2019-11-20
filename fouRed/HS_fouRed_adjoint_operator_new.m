function y = HS_fouRed_adjoint_operator_new(x, A, At, H, W, Sigma, Mask, aW)
% Phi'^t = Phi^t * Phi * F^t * Sigma = A^t * H * A * F^t * Sigma
% Complex -> Real
% ! Attention: H is a self-adjoint complete holographic matrix, or the 
% number of rows is equal to oversampled image size.

IFT2 = @(x) fftshift(ifft2(ifftshift(x))) * sqrt(numel(x));

% Parameters
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
[Ny, Nx] = size(At(zeros(No, 1)));

for ind = 1:c
    if iscell(H{ind})
        x1 = zeros(No, 1);
        for j = 1:length(H{ind})
            if exist('aW', 'var')
                xtmp = sqrt(aW{ind}{j}) .* x{ind}{j};
            else
                xtmp = x{ind}{j};
            end
            xtmp = Sigma{ind}{j} .* xtmp(:);
            tmp = zeros(size(Mask{ind}{j}));
            tmp(Mask{ind}{j}) = xtmp;
            if flagW
                x1(W{ind}{j}) = x1(W{ind}{j}) + H{ind}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            else
                x1 = x1 + H{ind}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            end
        end
    else
        if exist('aW', 'var')
            xtmp = sqrt(aW{ind}) .* x{ind};
        else
            xtmp = x{ind};
        end
        tmp = zeros(size(Mask{ind}));
        tmp(Mask{ind}) = Sigma{ind} .* xtmp(:);
        if flagW
            x1 = zeros(No, 1);
            x1(W{ind}) = x1(W{ind}) + H{ind} * A(real(IFT2(reshape(tmp, Ny, Nx)))); 
        else
            x1 = H{ind} * A(real(IFT2(reshape(tmp, Ny, Nx))));
        end
    end
    y(:,:,ind) = real(At(x1));
end

end
