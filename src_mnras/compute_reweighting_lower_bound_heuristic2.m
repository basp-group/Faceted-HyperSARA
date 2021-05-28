function [sig, sig_bar, mu, mu_bar] = ...
    compute_reweighting_lower_bound_heuristic2(Ny, Nx, ...
    nChannels, filters_length, nlevel, sigma_noise, ...
    algo_version, Qx, Qy, overlap_size, window_type, ...
    squared_operator_norm)
 
N = Ny*Nx;    % number of image pixels

% compute number of wavelet coefficients
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel);
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)

% compute sig and sig_bar
% sig = sqrt(N*sum((sigma_noise.^2)./squared_operator_norm)/s);
sig_w = sqrt(mean((sigma_noise.^2)./squared_operator_norm)/numel(filters_length));
mu_c = sqrt(2)*gamma((nChannels+1)/2)/gamma(nChannels/2);
sig_c = sqrt(nChannels - mu_c^2); 
sig = sig_w*(mu_c + 3*sig_c);

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

mu = sig;
mu_bar = sig_bar;
