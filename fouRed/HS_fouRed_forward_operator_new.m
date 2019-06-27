function y = HS_fouRed_forward_operator_new(x, A, At, H, Sigma, Mask, aW)
% Phi' = Sigma * F * Phi^t * Phi = Sigma * F * A^t * H * A

FT2 = @(x) fftshift(fft2(ifftshift(x)));

% Parameters
[Ny, Nx, c] = size(x);

%
for ind = 1:c
    x1 = A(real(x(:,:,ind)));
    x2 = H{ind} * x1;
    xtmp = FT2(real(At(x2)));
    xtmp = xtmp(:);
    xtmp = Sigma{ind} .* xtmp(Mask{ind});
    if exist('aW', 'var')
        y(:,ind) =  sqrt(cell2mat(aW{ind})) .* xtmp;
%         y(:,ind) =  sqrt(cell2mat(aW{ind})) .* B{ind}(x(:,:,ind));
    else
        y(:,ind) =  xtmp;
%         y(:,ind) =  B{ind}(x(:,:,ind));
    end
end

end
