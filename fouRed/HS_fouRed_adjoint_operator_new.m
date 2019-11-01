function y = HS_fouRed_adjoint_operator_new(x, A, At, H, Sigma, Mask, N, aW)
% Phi'^t = Phi^t * Phi * F^t * Sigma = A^t * H * A * F^t * Sigma

IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

% Parameters
[~, c] = size(x);
Ny = N(1);
Nx = N(2);

y = zeros(Ny, Nx, c);
%
for ind = 1:c
    if iscell(H{ind})
        % old blocking structure
        for j = 1:length(H{ind})
            if exist('aW', 'var')
                xtmp = sqrt(aW{ind}{j}) .* x{ind}{j};
            else
                xtmp = x{ind}{j};
            end
            xtmp = Sigma{ind}{j} .* xtmp(:);
            x1 = zeros(Ny * Nx, 1);
            x1(Mask{ind}{j}) = xtmp;
            x1 = reshape(x1, Ny, Nx);
            y(:,:,ind) = y(:,:,ind) + real(At(H{ind}{j} * A(real(IFT2(x1)))));
        end
    else
        if exist('aW', 'var')
            xtmp = sqrt(cell2mat(aW{ind})) .* x(:,ind);
        else
            xtmp = x(:,ind);
    %         y(:,:,ind) = Bt{ind}(sqrt(cell2mat(aW{ind})) .* x(:,ind));
    %     else
    %         y(:,:,ind) = Bt{ind}(x(:,ind));
        end
        x1 = zeros(Ny * Nx, 1);
        x1(cell2mat(Mask{ind})) = cell2mat(Sigma{ind}) .* xtmp(:);
        x1 = reshape(x1, Ny, Nx);
        y(:,:,ind) = real(At(H{ind} * A(real(IFT2(x1)))));
    end
end

end
