function psf_dr(datadir, ch, reduction_version, fouRed_gamma, fouRed_type, levelG, levelC)

addpath ../lib/operators/

if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end

Ny = 4096; Nx = 4096;

DRfilename = [datadir, '/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch), '.mat'];
fprintf('Read dimensionality reduction file: %s\n', DRfilename)
tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'A', 'At');
    
H = tmp.H{1,1};
W = tmp.W{1,1};
T = tmp.T{1,1};
Wm = tmp.Wm{1,1};
%     epsilon{i,1} = tmp.epsilon{1,1};
A = tmp.A;
At = tmp.At;

FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));
IFT2 = @(x) fftshift(ifft2(ifftshift(x))) * sqrt(numel(x));

% number of nodes
R = length(H);

% Variable flag for the case where W is present
flagW = 0;
if ~isempty(W)
    flagW = 1;
end

% number of over-sampled pixels
if flagW
    No = size(W{1}, 1);
else
    No = size(H{1}, 2);
end

dirac=zeros(Ny,Nx);
dirac(Ny/2,Nx/2)=1;
dirac_uv=A(dirac);
for q = 1:R        
    if reduction_version == 1
        if flagW
            tmp = FT2(real(At(H{q} * dirac_uv(W{q}))));
        else
            tmp = FT2(real(At(H{q} * dirac_uv)));
        end
        tmp = tmp(:);
        tmp2{q} = T{q} .* tmp(Wm{q});
    elseif reduction_version == 2
        if flagW
            tmp2{q} = T{q} .* (H{q} * dirac_uv(W{q}));
        else
            tmp2{q} = T{q} .* (H{q} * dirac_uv);
        end
    end
end

g_f2 = zeros(No, 1);
for q = 1 : R
    if reduction_version == 1   % H is self-adjoint in this case
        tmp = zeros(size(Wm{q}));
        tmp(Wm{q}) = T{q} .* tmp2{q};
        if flagW
            g_f2(W{q}) = g_f2(W{q}) + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
        else
            g_f2 = g_f2 + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
        end
    elseif reduction_version == 2
        if flagW
            g_f2(W{q}) = g_f2(W{q}) + H{q}' * (T{q} .* tmp2{q});
        else
            g_f2 = g_f2 + H{q}' * (T{q} .* tmp2{q});
        end
    end
end
psf = real(At(g_f2));
fitswrite(psf, ['./PSF_ch', num2str(ch), '.fits'])
fprintf("\nEffective channel: %i, Psf computation is finished\n", ch)