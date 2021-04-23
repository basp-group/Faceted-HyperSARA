function [sig, sig_bar, mu0, mu, mu_bar, m_bar, max_psf, l21_norm, nuclear_norm, l21_norm_x0, nuclear_norm_x0, dirty_image] = ...
    compute_reweighting_lower_bound(y, W, G, A, At, Ny, Nx, oy, ox, ...
    nChannels, wavelet_basis, filters_length, nlevel, sigma_noise, rw_type, algo_version, Qx, Qy, overlap_size, window_type, x0, Anorm)

%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
No = N*oy*ox; % size of oversampled Fourier space

% form dirty image (no normalization)
dirty_image = zeros([Ny, Nx, nChannels]);
for l = 1:nChannels
    temp = zeros(No, 1);
    for b = 1:numel(G{l})
        temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * y{l}{b};
    end
    dirty_image(:,:,l) = At(temp);
end

% generate noise transferred from data to image domain, compute maximum of 
% the psf in each channel
[B, max_psf] = create_dirty_noise(y, A, At, G, W, Nx, Ny, No, sigma_noise, 1234);

% evaluate nuclear norm of the dirty image -> regularization parameter
% dirty_image = dirty_image./reshape(max_psf, [1, 1, nChannels]);
dirty_image = dirty_image/Anorm; % Phi applied twice, and square operator norm for the normalisation
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
% B = B./reshape(max_psf, [1, nChannels]); %! normalize noise by the psf
B = B/Anorm; %! normalize noise by the squared norm of the operator
E = N*sum(var(B,0,1));

%! wavelet prior
sig = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));
mu0 = 1/sum(sig*log(d0/sig + 1));
if strcmp(rw_type, "ground_truth")
    mu = mu0;
else
    mu = 1/sum(sig*log(d1/sig + 1));
end

%! heuristic for the wavelet transform
% sig_w = sqrt(E/s);
% E_op = N*sum(var(B.*reshape(max_psf./Anorm_channel, [1, nChannels]),0,1));
% sig_w_op = sqrt(E_op/s);
% fprintf("Wavelet (upsilon^2): noise transfer = %e, heuristic (PSF) = %e, heuristic (op. norm) = %e\n", sig^2, sig_w^2, sig_w_op^2)

% [~,S0,~] = svd(B,'econ');
% sig_bar = std(abs(diag(S0)));
% sig_bar = std(abs(diag(U'*B*V)));

 %! need to compute operator norm per channel
%  [~,S0,~] = svd(reshape(X0, [N, nChannels]),'econ');
%  nuclear_norm_X0 = sum(abs(diag(S0)));            
%  [~, Psitw] = op_sp_wlt_basis_fhs(wlt_basis, nlevel, Ny, Nx);
%  Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
%  l21_norm_X0 = sum(sqrt(sum(Psit_full(x0).^2, 2))); 
%  fprintf("\\mu (X0) -> %e, (dirty, PSF) -> %e\n", ...
%      nuclear_norm_X0/l21_norm_X0, nuclear_norm/l21_norm)

% compute sig_bar
if strcmp(algo_version, 'hypersara')
    switch rw_type
        
    case "ground_truth"
        [U0,S0,V0] = svd(reshape(x0, [N, nChannels]),'econ');
        sig_bar = std(abs(diag(U0'*B*V0))); 

    case "dirty"
        [U0,S0,V0] = svd(reshape(dirty_image, [N, nChannels]),'econ');
        sig_bar = std(abs(diag(U0'*B*V0)));

    case "heuristic"
        sig_bar = sqrt(N*E/min(nChannels, N));
        [~,S0,~] = svd(reshape(dirty_image, [N, nChannels]),'econ');
    end        
    m_bar = sig_bar*sum(log(abs(diag(S0))/sig_bar + 1));
    mu_bar = 1/m_bar;    
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
    m_bar = zeros(Q, 1);
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
            dirty_im_q = w.*dirty_image(Io(q, 1)+1:Io(q, 1)+dims_o(q, 1), Io(q, 2)+1:Io(q, 2)+dims_o(q, 2), :);
            [U0,S0,V0] = svd(reshape(dirty_im_q, [Noq, nChannels]),'econ');
            sig_bar(q) = std(abs(diag(U0'*bq*V0)));
            
        case "heuristic"
            sig_bar(q)= norm(w(:))^2*E/(N*min(nChannels, Noq));
             dirty_im_q = w.*dirty_image(Io(q, 1)+1:Io(q, 1)+dims_o(q, 1), Io(q, 2)+1:Io(q, 2)+dims_o(q, 2), :);
            [~,S0,~] = svd(reshape(dirty_im_q, [Noq, nChannels]),'econ');
        end
        m_bar(q) = sig_bar(q)*sum(log(abs(diag(S0))/sig_bar(q) + 1));
    end    
    mu_bar = 1/sum(m_bar);
end

end