function f = compute_sara_prior_distributed(x, Psit, s)
% Compute the l21-norm of the wavelet cofficients associated with the SARA 
% dictionary.
%
% Compute :math:`\Vert \boldsymbol{\Psi}^\dag \mathbf{X} \Vert_{2,1}` for 
% an input wideband image :math:`\mathbf{X} \in \mathbb{R}^{N\times L}`, 
% with :math:`\boldsymbol{\Psi}^\dag` the SARA dictionary :cite:`Carrillo2012`.
%
% Parameters
% ----------
% x : array (3d)
%     Wideband image cube [N(1), N(2), L].
% Psit : lambda
%     SARA dictionary @[1]
% s : int
%     Number of wavelet coefficients (per channel) [1]
%
% Returns
% -------
% f : double
%     Value of :math:`\Vert \boldsymbol{\Psi}^\dag \mathbf{X} \Vert_{2,1}`.
%

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

r = zeros(s, size(x, 3));

for l = 1:size(x, 3)
    r(:, l) = Psit(x(:, :, l)); 
end

d_ = sum(abs(r).^2,2);
d = gplus(d_);
f = sum(sqrt(d));

end
    