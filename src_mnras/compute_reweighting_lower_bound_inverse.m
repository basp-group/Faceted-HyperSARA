function [sig, sig_bar, mu0, mu, mu_bar, mu_bar_full, max_psf, l21_norm, nuclear_norm, l21_norm_x0, nuclear_norm_x0, dirty_image_precond] = ...
    compute_reweighting_lower_bound_inverse(y, W, G, aW, A, At, Ny, Nx, oy, ox, ...
    nChannels, wavelet_basis, filters_length, nlevel, sigma_noise, rw_type, algo_version, Qx, Qy, overlap_size, window_type, x0, Anorm, squared_operator_norm)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% form dirty image (no normalization)
dirty_image_precond = zeros([Ny, Nx, nChannels]);
for l = 1:nChannels
    temp = zeros(No, 1);
    for b = 1:numel(G{l})
        temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * (aW{l}{b}.*y{l}{b}); % ... sqrt(aW{i}{j}) .* sqrt(aW{i}{j}) .* y{l}{b}
    end
    dirty_image_precond(:,:,l) = At(temp);
end

% generate noise transferred from data to image domain, compute maximum of 
% the psf in each channel
[B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No, sigma_noise, 1234);

% evaluate nuclear norm of the dirty image -> regularization parameter
% dirty_image_precond = dirty_image_precond./reshape(max_psf, [1, 1, nChannels]);
dirty_image_precond = dirty_image_precond/Anorm; % Phi applied twice, and squared preconditioned operator norm for the normalisation
[~,S0,~] = svd(reshape(dirty_image_precond, [N, nChannels]),'econ');
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
d1 = sqrt(sum(Psit_full(dirty_image_precond).^2, 2));
d0 = sqrt(sum(Psit_full(x0).^2, 2));
l21_norm = sum(d1); 
l21_norm_x0 = sum(d0); 

% compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
% SARA space) involved in the reweighting scheme
B = B/squared_operator_norm; %! normalize noise by the squared norm of the operator

sig = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));
mu0 = 1/l21_norm_x0;
mu = 1/l21_norm;
mu_bar_full = 1/nuclear_norm; % using the inverse

% compute sig_bar
if strcmp(algo_version, 'hypersara')
    switch rw_type
        
    case "ground_truth"
        [U0,~,V0] = svd(reshape(x0, [N, nChannels]),'econ');
        sig_bar = std(abs(diag(U0'*B*V0))); 

    case "dirty"
        [U0,~,V0] = svd(reshape(dirty_image_precond, [N, nChannels]),'econ');
        sig_bar = std(abs(diag(U0'*B*V0)));

    end
    mu_bar = mu_bar_full;
else
    Q = Qx*Qy;
    B = reshape(B, [Ny, Nx, nChannels]);

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
    mu_bar = 0;
    for q = 1:Q
        [qy, qx] = ind2sub([Qy, Qx], q);
        w = generate_weights(qx, qy, Qx, Qy, window_type, dims(q,:), dims_o(q,:), overlap_size);
        
        Noq = prod(dims_o(q, :));
        Bq = w.*B(Io(q, 1)+1:Io(q, 1)+dims_o(q, 1), Io(q, 2)+1:Io(q, 2)+dims_o(q, 2), :);
        bq = reshape(Bq, [Noq, nChannels]);

        switch rw_type
            
        case "ground_truth"
            % SVD space of X0
            X0q = w.*x0(Io(q, 1)+1:Io(q, 1)+dims_o(q, 1), Io(q, 2)+1:Io(q, 2)+dims_o(q, 2), :);
            [U0,S0,V0] = svd(reshape(X0q, [Noq, nChannels]),'econ');
            sig_bar(q) = std(abs(diag(U0'*bq*V0))); 
            
        case "dirty"
            % SVD space of Xdirty
            dirty_im_q = w.*dirty_image_precond(Io(q, 1)+1:Io(q, 1)+dims_o(q, 1), Io(q, 2)+1:Io(q, 2)+dims_o(q, 2), :);
            [U0,S0,V0] = svd(reshape(dirty_im_q, [Noq, nChannels]),'econ');
            sig_bar(q) = std(abs(diag(U0'*bq*V0)));
        end
        mu_bar = mu_bar + sum(abs(diag(S0)));   
    end
    mu_bar = 1/mu_bar;
end