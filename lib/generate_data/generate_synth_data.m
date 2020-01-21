function [Y, V, W, Omega, u_ab, v_ab, y, x, sp_scale, U, D, N, K, na] = generate_synth_data(S, J, N, K, P, T, F, dl, input_snr, cov_type, A, im_choice)
% Generate synthetic RI data with smoothly spatial- and time-varying 
% DDE kernels.
%-------------------------------------------------------------------------%
%%
% Input:
% > S         : spatial dimension of the DDE kernels (square kernels [S, S])
% > J         : size of the gridding kernels (square kernels [J, J])
% > N         : image dimension [2, 1]
% > K         : size of the spatial Fourier space [2, 1]
% > P         : size of the DDEs' temporal support
% > T         : number of snapshots
% > F         : size of the temporal Fourier space
% > dl        : pixel size
% > input_snr : SNR value (in dB)
% > cov_type  : type of u-v coverage ('vlaa', 'meerkat')
% > A         : value of the DDE standard deviation (randomly-geerated ground truth)
% > im_choice : image selected ('M31', 'W28')
%
% Output:
% < Y     : data with redundancy [na, na, T]
%           (Nan on the diag., y_ab = conj(y_ab) for a > b)
% < V     : spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na, na, T]
% < W     : spatial NUFFT gridding kernels [J2, T]
% < Omega : u-v components with redundancy (similar to Y) [na, na, 2, T],
%           Omega(:, :, 1, t) -> v_ab(t)
%           Omega(:, :, 2, t) -> u_ab(t)
% < u_ab  : u components [M, T]
% < v_ab  : v components [M, T]
% < y     : data vector [M*T, 1]
% < x     : ground truth image
% < sp_scale : NUFFT scaling coefficients
% < U     : DDE kernels in the spatial and temporal frequency domain [S2, P]
% < D     : DDE kernels in the spatial frequency domain [S2, T]
% < N     : image size
% < K     : size of the spatial Fourier space
% < na    : number of antennas
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

% Load image
x = fitsread([im_choice,'.fits']); 
if any(size(x) > N)
    x = imresize(x, N);
end
x = (x+abs(x))./2;
x = x./max(x(:));

% Parameters
S2 = S*S;
J2 = J*J;

% U-V coverage & antenna positions
[u_ab, v_ab, na] = generate_uv_coverage(T, dl, cov_type);
M = na*(na-1)/2;

% Generate ground truth DDEs
cs = floor(S2/2) + 1; % center of the spatial support
p = floor(P/2) + 1;
U = (A*(randn(S2, na, P) + 1i*randn(S2, na, P))/P)*sqrt(F);
U(:, :, p) = 0;
U(cs, :, p) = sqrt(F);
D = computeD(U, F, [], [], T);

% Spatial gridding coefficients and associated frequencies
V = zeros(J2, na, na, T);
W = cell(T, 1);
Omega = zeros(na, na, 2, T);
indices = getIndicesLow(na); % indices for the lower part of Y

parfor t = 1:T % "for" is enough for small T             
    % spatial gridding coefficients
    st = compute_interp_coeffs([v_ab(:, t), u_ab(:, t)], N, [J, J], K, N/2); % [Ny/2, Nx/2], N = K/2 [Mt, J2] 
    W{t} = st.uu; % [J2, M] 
    v_temp = flipud(conj(st.uu));
    Z = zeros(na, na, J2);
    w_temp = zeros(na);
    for j = 1:J2
        w_temp(indices) = st.uu(j,:); % explained by the symmetry of the coefficients
        w_temp = w_temp.';
        w_temp(indices) = v_temp(j,:);
        Z(:, :, j) = w_temp;
    end
    V(:, :, :, t) = permute(Z, [3, 1, 2]); % [J2, na, na]
end
st = compute_interp_coeffs([v_ab(:, 1), u_ab(:, 1)], N, [J, J], K, N/2);
sp_scale = st.sn;

% Created duplicate frequency matrices [na, na, 2, T]
for t = 1:T
    w_temp = NaN(na);
    w_temp(indices) = v_ab(:, t);
    w_temp = w_temp.';
    w_temp(indices) = -v_ab(:, t);
    Omega(:,:,1,t) = w_temp; % -w_temp + w_temp.';
    w_temp(indices) = u_ab(:, t);
    w_temp = w_temp.';
    w_temp(indices) = -u_ab(:, t);
    Omega(:,:,2,t) = w_temp; % -w_temp + w_temp.'; % [na, na, 2] for each t
end
    
% Data generation
y = zeros(M, T);
Y = zeros(na, na, T);

for t = 1:T
    % Create operator G [for each time instant]
    G = createGnufft_T2(D(:, :, t), D(:, :, t), [v_ab(:, t), u_ab(:, t)], K, S, J, W{t});
    A_ = @(x) G*so_fft2(x, K, sp_scale);
    y0 = A_(x);
    sigma_noise = sqrt(sum(abs(y0).^2)*10^(-input_snr/10));
    noise = (randn(size(y0)) + 1i*randn(size(y0)))*sigma_noise/sqrt(2);
    y(:,t) = y0 + noise;
    % Create duplicate data Y [na, na]
    Y(:,:,t) = createY(y(:,t), na);
end

% Replace invalid entries by NaN
Y(Y == 0) = NaN; 
V(V == 0) = NaN;

end
