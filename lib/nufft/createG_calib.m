function [G, scale, ll, v, V] = createG_calib(D, om, N, K, J, S, nshift)
% Create the gridding matrix G including the DDE kernels.
%-------------------------------------------------------------------------%
% Input:
% > D      : DDEs kernels for a single time instant [S2, na]
% > om     : normalized u-v frequencies [M, 2] 
% > N      : image size [1, 2]
% > K      : size of the spatial Fourier space [1, 2] 
% > J      : size of the gridding kernels
% > S      : size of the DDE kernels (in the spatial Fourier domain)
% > nshift : shift in Fourier space [1, 2]
%
% Output:
% < G     : gridding matrix [M, prod(K)] (sparse matrix)
% < scale : scaling coefficients for the NUFFT [N]
% < ll    : position of the nonzero values of G
% < v     : convolutions between the DDE kernels [Q2, M]
% < V     : values contained in G [Q2, J2, M]
%-------------------------------------------------------------------------%
%% 
% [06/08/2019], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%
% Compute gridding coefficients  
% compute_interp_coeffs(om, Nd, Jd, Kd, n_shift)
st = compute_interp_coeffs(om, N, [J, J], K, nshift); % st.uu of size [J^2, M]
Q  = 2*S - 1;             % size of the kernels after convolution
Q2 = Q^2;
J2 = J^2;

Qprime = floor(Q/2);      % Q is always odd (one possibility only)
tmp1 = (-Qprime:Qprime).';
tmp1 = tmp1(:, ones(Q,1));

[~, na, ~] = size(D);    % [S2, na] number of antennas at time t
M_true  = na*(na-1)/2; % number of acquisitions -> check value...
M = size(om, 1);       % M = M_true if all the measurements are present, M < M_true if some of the data have been flagged: need to add flagging in this case
T = size(D, 3);
v = zeros(Q^2, M); % convolution values (stored in column for each pair)

%% Perform 2D convolutions and gridding using D1 and D2 (to be possibly performed in parallel)
for t = 1:T
    q = 0; % global counter
    for alpha = 1:na-1
        for beta = alpha+1:na % modify the double loop to exclusively select the appropriate elements, apply nonzeros on W
            % 2D convolutions
            q = q+1;
            v(:,q) = reshape(conv2(rot90(reshape(D(:,alpha,t),[S,S]),2),reshape(conj(D(:,beta,t)),[S,S])), [Q^2,1]); % only select the appropriate entries...
        end
    end
end

% % shortcut for testing purposes
% vt = reshape(conv2(rot90(reshape(D(:,1,1),[S,S]),2),reshape(conj(D(:,2,1)),[S,S])), [Q^2,1]);
% v = vt(:, ones(M, 1));

% Generate indices in the sparse G matrix
if rem(J,2) > 0 % odd
   c0 = round(om.*K/(2*pi)) - (J+1)/2; % [M, 2]
else
   c0 = floor(om.*K/(2*pi)) - J/2; 
end
kdy = bsxfun(@plus, (1:J).', c0(:,1).'); % [J M]
kdx = bsxfun(@plus, (1:J).', c0(:,2).'); % [J M]
ii = mod(bsxfun(@plus, tmp1(:), reshape(kdy, [1,J,M])), K(1)) + 1; % [Q2, J, n] % row indices of the elements within each area, 
                                                                   % whose leftmost element row indices are given above
jj = mod(bsxfun(@plus, reshape(tmp1.', [Q2,1]), reshape(kdx, [1,J,M])), K(2)) + 1; % [Q2, J, M] % column indices ...
ll = reshape(bsxfun(@plus, reshape((jj-1)*K(1), [Q2,1,J,M]), reshape(ii, [Q2,J,1,M])), [Q2,J2,M]);

% Duplicate values to have all the convolutions centered in the different elements
V = bsxfun(@times, reshape(v, [Q2, 1, M]), reshape(st.uu, [1, J2, M])); % [Q2, J2, M]   
V(isnan(V(:))) = []; % [J2*M, 1] there are zeros in W at the positions where the 
                 % measurements are missing, right size once the zeros are 
                 % filtered

% Generate row indices (within G) [remark: Jt^2 nz elements per row]
kk = bsxfun(@times, ones(J2*Q2,1), 1:M); % [J2*Q2, M]
% try
G = sparse(kk(:), ll(:), V(:), M, K(1)*K(2));  % [M, prod(K)]
scale = st.sn;
% catch 
%     keyboard
% end
end
