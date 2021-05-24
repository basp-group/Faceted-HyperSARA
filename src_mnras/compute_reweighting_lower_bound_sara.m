function [sig, max_psf, mu, l11_norm, l11_norm_x0, dirty_image] = ...
    compute_reweighting_lower_bound_sara(y, W, G, aW, A, At, Ny, Nx, oy, ox, ...
    wavelet_basis, filters_length, nlevel, sigma_noise, rw_type, x0, Anorm, squared_operator_norm,xapprox,noise_transfer,regtype)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% generate noise transferred to image domain
if strcmp(noise_transfer, "precond")
    [B, max_psf] = create_dirty_noise_precond(y, A, At, G, aW, W, Nx, Ny, No, sigma_noise, 1234);
    B = B/Anorm; %! normalize noise by the squared norm of the operator
else
    [B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No, sigma_noise, 1234);
    B = B/squared_operator_norm; %! normalize noise by the squared norm of the operator
end

% form dirty image 
if strcmp(xapprox, "precond")
    temp = zeros(No, 1);
    for b = 1:numel(G{1})
        temp(W{1}{b}) = temp(W{1}{b}) + G{1}{b}' * (aW{1}{b}.*y{1}{b});
    end
    dirty_image = At(temp);
    dirty_image = dirty_image/Anorm; % Phi applied twice, and square operator norm for the normalisation
else
    temp = zeros(No, 1);
    for b = 1:numel(G{1})
        temp(W{1}{b}) = temp(W{1}{b}) + G{1}{b}' * y{1}{b};
    end
    dirty_image = At(temp);
    dirty_image = dirty_image/squared_operator_norm;
end

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel); % suppose Dirac is the last entry
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
d0 = abs(Psit_full(x0));
d1 = abs(Psit_full(dirty_image));
l11_norm_x0 = sum(d0);
l11_norm = sum(d1);

% compute sig (estimate of the "noise level" in SARA space) involved in the
% reweighting scheme
switch regtype
    case "inv"
        sig = std(abs(Psit_full(reshape(B, [Ny, Nx]))));
        if strcmp(rw_type, "ground_truth")
            mu = 1/l11_norm_x0;
        else
            mu = 1/l11_norm;
        end
    case "log"
        sig = std(abs(Psit_full(reshape(B, [Ny, Nx]))));
        if strcmp(rw_type, "ground_truth")
            mu = 1/sum(sig*log(d0/sig + 1));
        else
            mu = 1/sum(sig*log(d1/sig + 1));
        end
    case "heuristic"
        sig = sqrt(N*(sigma_noise^2)/(squared_operator_norm*s));
        mu = sig;
    otherwise
        error("Unknown regularization type")
end
