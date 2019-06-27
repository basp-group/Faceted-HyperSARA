function y = HS_fouRed_adjoint_operator_new(x, A, At, H, Sigma, Mask, N, aW)
% Phi'^t = Phi^t * Phi * F^t * Sigma = A^t * H * A * F^t * Sigma

IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

% Parameters
[~, c] = size(x);
Ny = N(1);
Nx = N(2);

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
    y(:,:,ind) = At(H{ind} * A(IFT2(x1)));
end

end
