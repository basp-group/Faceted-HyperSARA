function [sig, sig_bar, mu0, mu, mu_bar, mu_bar_full, l21_norm, nuclear_norm, l21_norm_x0, nuclear_norm_x0, dirty_image] = ...
    compute_reweighting_lower_bound_inverse_heuristic(y, W, G, aW, At, Ny, Nx, ...
    oy, ox, nChannels, wavelet_basis, filters_length, nlevel, sigma_noise, ...
    algo_version, Qx, Qy, overlap_size, window_type, x0, squared_precond_operator_norm, ...
    squared_operator_norm, xapprox)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% form dirty image (no normalization)
dirty_image = zeros([Ny, Nx, nChannels]);
if strcmp(xapprox, "precond")
    for l = 1:nChannels
        temp = zeros(No, 1);
        for b = 1:numel(G{l})
            temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * (aW{l}{b}.*y{l}{b});
        end
        dirty_image(:,:,l) = At(temp);
    end
    dirty_image = dirty_image./reshape(squared_precond_operator_norm, [1, 1, nChannels]); % Phi applied twice, and square operator norm for the normalisation
else
    for l = 1:nChannels
        temp = zeros(No, 1);
        for b = 1:numel(G{l})
            temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * y{l}{b};
        end
        dirty_image(:,:,l) = At(temp);
    end
    dirty_image = dirty_image./reshape(squared_operator_norm, [1, 1, nChannels]);
end


% evaluate nuclear norm of the dirty image -> regularization parameter
[~,S0,~] = svd(reshape(dirty_image, [N, nChannels]),'econ');
nuclear_norm = sum(abs(diag(S0)));
[~,S0,~] = svd(reshape(x0, [N, nChannels]),'econ');
nuclear_norm_x0 = sum(abs(diag(S0)));

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis_fhs(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel);
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)

% evaluate l21 norm of the dirty image -> regularization parameter
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
d1 = sqrt(sum(Psit_full(dirty_image).^2, 2));
d0 = sqrt(sum(Psit_full(x0).^2, 2));
l21_norm = sum(d1); 
l21_norm_x0 = sum(d0); 

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
mu0 = 1/l21_norm_x0;
mu = 1/l21_norm;
mu_bar_full = 1/nuclear_norm; % using the inverse
mu_bar = mu_bar_full;

% sig = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));
% sig = sqrt(N*sum((sigma_noise.^2)./squared_operator_norm)/s);
sig = sqrt(sum((sigma_noise.^2)./squared_operator_norm)/numel(filters_length));

% compute sig_bar
if strcmp(algo_version, 'hypersara')
    sig_bar = sqrt(Nx*Ny*sum(sigma_noise.^2./squared_operator_norm)/min(Nx*Ny, nChannels));
else
    Q = Qx*Qy;
    rg_y = split_range(Qy, Ny);
    rg_x = split_range(Qx, Nx);
    I = zeros(Q, 2);
    dims = zeros(Q, 2);
    for qx = 1:Qx
        for qy = 1:Qy
            q = (qx-1)*Qy+qy;
            I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
            dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
        end
    end
    clear rg_y rg_x;

    rg_yo = split_range(Qy, Ny, overlap_size(1));
    rg_xo = split_range(Qx, Nx, overlap_size(2));
    Io = zeros(Q, 2);
    dims_o = zeros(Q, 2);
    for qx = 1:Qx
        for qy = 1:Qy
            q = (qx-1)*Qy+qy;
            Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
            dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
        end
    end
    clear rg_yo rg_xo;

    sig_bar = zeros(Q, 1);
    for q = 1:Q
        [qy, qx] = ind2sub([Qy, Qx], q);
        w = generate_weights(qx, qy, Qx, Qy, window_type, dims(q,:), dims_o(q,:), overlap_size);
        
        Noq = prod(dims_o(q, :));
        sig_bar(q) = sqrt(sum(w(:).^2)*sum(sigma_noise.^2./squared_operator_norm)/min(Noq, nChannels));
    end
end
