function [G, scale] = createG(om, J, N, K, nshift)
% Create the gridding matrix G to perform the (spatial) NUFFT.
%-------------------------------------------------------------------------%
% Input:
% > om     : normalized u-v frequencies [M, 2] 
% > J      : size of the gridding kernels
% > N      : image size [1, 2]
% > K      : size of the spatial Fourier space [1, 2] 
% > nshift : number of time instants
%
% Output:
% < G     : gridding matrix [M, prod(K)]
% < scale : scaling coefficients for the NUFFT [N]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%
% Compute gridding coefficients  
% compute_interp_coeffs(om, Nd, Jd, Kd, n_shift)
st = compute_interp_coeffs(om, N, [J, J], K, nshift); % st.uu of size [J^2, M]
M = size(om, 1);
J2 = J*J;

% Generate indices in the sparse G matrix
if rem(J,2) > 0 % odd
   c0 = round(om.*K/(2*pi)) - (J+1)/2; % [M, 2]
else
   c0 = floor(om.*K/(2*pi)) - J/2; 
end
% keyboard
kdy = mod(bsxfun(@plus, (1:J).', c0(:,1).'), K(1)) + 1; % [J M]
kdx = mod(bsxfun(@plus, (1:J).', c0(:,2).'), K(2)) + 1; % [J M] row indices of the elements within each area, whose leftmost element row indices are given above
x = reshape(bsxfun(@plus, reshape((kdx-1)*K(1), [1,J,M]), reshape(kdy, [J,1,M])), [J2,M]);

y = bsxfun(@times, ones(J2, 1), 1:M); % [J2 M]
% Create the sparse matrix Gt of size [T, F] -> save for the following 
% temporal interpolations
G = sparse(y(:), x(:), st.uu, M, prod(K)); % trim down Gt for application to U (remove unused columns...)
scale = st.sn;
end
