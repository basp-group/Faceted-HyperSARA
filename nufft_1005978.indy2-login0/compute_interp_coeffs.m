function st = compute_interp_coeffs(om, Nd, Jd, Kd, n_shift)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% default/recommended interpolator is minmax with KB scaling factors
M = size(om,1);
dd = length(Nd);

st.alpha = cell(dd,1);
st.beta = cell(dd,1);
for id = 1:dd
    [st.alpha{id}, st.beta{id}] = nufft_alpha_kb_fit(Nd(id), Jd(id), Kd(id)); % ok
end

st.tol = 0;

% scaling factors: "outer product" of 1D vectors
st.sn = 1;
for id=1:dd
    tmp = nufft_scale(Nd(id), Kd(id), st.alpha{id}, st.beta{id}); % ok
	st.sn = st.sn(:) * tmp';
end
if length(Nd) > 1
	st.sn = reshape(st.sn, Nd); % [(Nd)]
else
	st.sn = st.sn(:); % [(Nd)]
end

% keyboard
% [J? M] interpolation coefficient vectors.  will need kron of these later
ud = cell(dd,1);
for id=1:dd
	N = Nd(id);
	J = Jd(id);
	K = Kd(id);
    alpha = st.alpha{id};
    beta = st.beta{id};
    T = nufft_T(N, J, K, st.tol, alpha, beta, 'false'); % [J? J?] % should be ok
    [r, arg] = nufft_r(om(:,id), N, J, K, alpha, beta, 'false'); % [J? M]  % should be ok
    c = T*r;	clear T r

	gam = 2*pi/K;
	phase_scale = 1i * gam * (N-1)/2;

	phase = exp(phase_scale * arg); % [J? M] linear phase
	ud{id} = phase .* c; % [J? M]
end, clear c arg gam phase phase_scale koff N J K

% kk = kd{1}; % [J1 M] % [PA] values of the indices locations
uu = ud{1}; % [J1 M] % [PA] values of the u coefficients
Jprod = Jd(1);
for id = 2:dd
% 	kk = block_outer_sum(kk, kd{id}); % outer sum of indices
% 	kk = reshape(kk, Jprod, M);
% 	uu = block_outer_prod(uu, ud{id}); % outer product of coefficients
    uu = bsxfun(@times, reshape(uu, [Jprod, 1, M]), reshape(ud{id}, [1, Jd(id), M]));
    Jprod = prod(Jd(1:id));
    uu = reshape(uu, Jprod, M);
end % now kk and uu are [*Jd M]

% apply phase shift
% pre-do Hermitian transpose of interpolation coefficients
phase = exp(1i * (om * n_shift(:))).'; % [1 M] % [PA] use of n_shift to shift the frequencies given in om
% st.uu = conj(uu) .* phase(ones(1,prod(Jd)),:); % [*Jd M] % [PA] use bsxfun instead of duplicating entries... 
st.uu = bsxfun(@times, conj(uu), phase); % [*Jd M]

% mm = [1:M]; mm = mm(ones(prod(Jd),1),:); % [*Jd M]
% make sparse matrix, ensuring arguments are double for stupid matlab
% st.p = sparse(mm(:), double(kk(:)), double(uu(:)), M, prod(Kd));

end

