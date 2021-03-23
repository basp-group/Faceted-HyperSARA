function [sig, sig_bar, max_psf, l21_norm, nuclear_norm, dirty_image] = ...
    compute_reweighting_lower_bound_dr(yTp, Wp, Tp, Hp, Ap, Atp, Ny, Nx, No, ...
    nChannels, wavelet_basis, filters_length, nlevel, Q, cell_c_chunks)
%! TO BE DOCUMENTED
%! make sure the rng always starts from the same value for reproducibility 

N = Ny*Nx;    % number of image pixels
n_data_workers = length(cell_c_chunks);

spmd
    if labindex > Q % data cores   
        %! this instructions would need to be fixed (the function does not take as an input a structure similar to the data y, which is a problem)    
        % local_dirty_image = HS_operatorGtPhi_t(yTp, Hp, Wp, Atp, Tp, []);
        local_dirty_image = create_dirty_image_dr(y, At, H, T, W, Nx, Ny, No);
        [local_B, local_max_psf] = create_dirty_noise_dr(yTp, Ap, Atp, Hp, Tp, Wp, Nx, Ny, No);
    end
end

% retrieve results from data workers
max_psf = zeros(nChannels, 1);
dirty_image = zeros([Ny, Nx, nChannels]);
B = zeros([N, nChannels]);
for k = 1:n_data_workers
    dirty_image(:,:,cell_c_chunks{k}) = local_dirty_image{Q+k};
    B(:,cell_c_chunks{k}) = local_B{Q+k};
    max_psf(cell_c_chunks{k}) = local_max_psf{Q+k};
end

% estimate mu: ratio between nuclear and l21 norm priors applied to
% the dirty image
[~,S0,~] = svd(reshape(dirty_image, [N, nChannels]),'econ');
nuclear_norm = sum(abs(diag(S0)));

% set-up global SARA dictionary
dwtmode('zpd')
[~, Psitw] = op_sp_wlt_basis(wavelet_basis, nlevel, Ny, Nx);
[~, s] = n_wavelet_coefficients(filters_length(1:end-1), [Ny, Nx], 'zpd', nlevel);
s = s+N; % total number of SARA coefficients (adding number of elements from Dirac basis)
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
l21_norm = sum(sqrt(sum(Psit_full(dirty_image).^2, 2))); 

% compute reweighting lower bound (estimate of the "noise level" in "SVD" 
% and SARA space) involved in the reweighting scheme
B = B./reshape(max_psf, [1, nChannels]);
[~,S0,~] = svd(B,'econ');
sig = std(diag(S0));
sig_bar = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));

end
