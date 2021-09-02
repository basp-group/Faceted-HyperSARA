function f = nuclear_norm(x)
% Compute the nuclear norm of a wideband image.
%
% Compute :math:`\Vert \mathbf{X} \Vert_*` for :math:`\mathbf{X} \in 
% \mathbb{R}^{N \times L}`.
%
% Parameters
% ----------
% x : array 3d()
%     Wideband image cube [N(1), N(2), L].
%
% Returns
% -------
% f : double
%     Nuclear norm of the wideband image.
%

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

[~,D,~] = svd(reshape(x, [numel(x)/size(x, 3), size(x, 3)]),'econ');
f = sum(abs(diag(D)));

end
