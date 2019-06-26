function y = HS_fouRed_adjoint_operator(x, A, At, H, Sigma, Mask, N, No, W, aW)
% Phi'^t = Phi^t * Phi * F^t * Sigma = A^t * H * A * F^t * Sigma

IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

% Parameters
[~, c] = size(x);
Ny = N(1);
Nx = N(2);
if exist('No', 'var')
    Noy = No(1);
    Nox = No(2);
end
%
for ind = 1:c
    if exist('aW', 'var')
        xtmp = sqrt(cell2mat(aW{ind})) .* x(:,ind);
    else
        xtmp = x(:,ind);
%         y(:,:,ind) = Bt{ind}(sqrt(cell2mat(aW{ind})) .* x(:,ind));
%     else
%         y(:,:,ind) = Bt{ind}(x(:,ind));
    end
    x1 = zeros(Ny * Nx, 1);
    x1(Mask{ind}) = Sigma{ind} .* xtmp(:);
    x1 = reshape(x1, Ny, Nx);
    x1 = A(IFT2(x1));
    if exist('No', 'var') && exist('W', 'var')
        tmp = H{ind} * x1(W{ind});
        x2 = zeros(Noy*Nox,1);
        x2(W{ind}) = tmp;
    else
        x2 = H{ind} * A(IFT2(x1));
    end
    y(:,:,ind) = At(x2);
end

end
