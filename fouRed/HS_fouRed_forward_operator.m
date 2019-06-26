function y = HS_fouRed_forward_operator(x, A, At, H, Sigma, Mask, No, W, aW)
% Phi' = Sigma * F * Phi^t * Phi = Sigma * F * A^t * H * A

FT2 = @(x) fftshift(fft2(ifftshift(x)));

% Parameters
[Ny, Nx, c] = size(x);
if exist('No', 'var')
    Noy = No(1);
    Nox = No(2);
end

%
for ind = 1:c
    x1 = A(real(x(:,:,ind)));
    if exist('No', 'var') && exist('W', 'var')
        tmp = H{ind} * x1(W{ind});
        x2 = zeros(Noy*Nox,1);
        x2(W{ind}) = tmp;
    else
        x2 = H{ind} * x1;    
    end
    xtmp = FT2(real(At(x2)));
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
